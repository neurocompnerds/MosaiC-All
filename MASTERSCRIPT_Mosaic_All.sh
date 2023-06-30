#!/bin/bash
# MasterScript: Phase 1 of mosaic variant finding pipeline, which includes
# 1. 
# 2. GATK-HC: germline variant calling in each family
# 3. Mutect2 and FilterMutect2: Parents,Probands and Siblings (if available)
# 4. MosaicHunter: Parents,Probands and Siblings (if available)
# 5. MosaicForecast on Mutect2 variant callset, followed by Filter
# UniAdelaide-HPC friendly
# Date: 9th June 2023
# 

usage()
{
echo "#MasterScript: Phase 1 of mosaic variant finding pipeline, which includes
# 1. Coverage Analysis of every bam file
# 2. GATK-HC: germline variant calling in each family
# 3. Mutect2: Parents and Probands and Siblings (if available)
# 4. MosaicHunter: Parents and Probands and Siblings (if available)
# 5. MosaicForecast on Mutect2 variant callset
#
#
# Usage $0 -s $sampleID_list -o $Output_folder -c $Config_File | [ - h | --help ]
#
# Options
#-s REQUIRED sampleID.list (tab-delimited columns $BAMdir,$ProbandID,$Gender,$Mother,$Father)
#-o REQUIRED Output_directory (all variant calls_output in a single output directory)
#-c REQUIRED ConfigFile (i.e /hpcfs/groups/phoenix-hpc-neurogenetics/Nandini/Mosaic-All/Mosaic-S/Mosaic-All.config)
#
# -h or --help  Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
#
# Original: Nandini Sandran, 9/6/2023
# Modified: (Date; Name; Description)
# 

#TEST (delete after finalising the scripts)
#MAIN=/hpcfs/groups/phoenix-hpc-neurogenetics/Nandini
#Test01: bash $MAIN/Mosaic-All/Mosaic-S/MasterScript_Mosaic_All.sh -s $MAIN/Mosaic-All/SampleID -o $MAIN/Mosaic-All/Outputs -p $MAIN/Mutect2_ReCalling_batch1/PON/PON_Batch01_hs37dh_GAparents.vcf
#Test02: bash $MAIN/Mosaic-All/Mosaic-S/MasterScript_Mosaic_All.sh -s $MAIN/Mosaic-All/SampleID -o $MAIN/Mosaic-All/Outputs -c $MAIN/Mosaic-All/Mosaic-S/Mosaic-All.config

"
}

## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -s )                    shift
                                        SAMPLELIST=$1
                                        ;;
                -o )                    shift
                                        OUTDIR=$1
                                        ;;
                -c )                    shift
                                        CONFIG_FILE=$1
                                        ;;
                 * )                    usage
                                        exit 1
        esac
        shift
done


## Define Directories## (how to change this accordingly)
SCRIPTDIR="/hpcfs/groups/phoenix-hpc-neurogenetics/Nandini/Mosaic-All/Mosaic-S"

if [ ! -d "${OUTDIR}" ]; then
    mkdir -p ${OUTDIR}
	echo "## INFO: output directory created, you'll find all of the outputs and log files in here: ${OUTDIR}" >> $LOGDIR/$ProbandID.pipeline.log
fi

#Array from list of Samples (ignoring the header of the file)
mapfile -t SAMPLEID < <(tail -n +2 "$SAMPLELIST")

#modules
module purge
module load BCFtools/1.17-GCC-11.2.0

source $CONFIG_FILE

# Iteration for variant calling starts here
for SAMPLEID in "${SAMPLEID[@]}"; do

    #Defining variables from each row
		BamDIR=$(awk '{print $1}' <<< "$SAMPLEID ")
    		ProbandID=$(awk '{print $2}' <<< "$SAMPLEID ")
    		ProbandGender=$(awk '{print $3}' <<< "$SAMPLEID ")
    		MotherID=$(awk '{print $4}' <<< "$SAMPLEID ")
    		FatherID=$(awk '{print $5}' <<< "$SAMPLEID ")

		echo "Pipeline for $ProbandID,$MotherID,$FatherID in $BamDIR" >> $OUTDIR/$ProbandID.pipeline.log

    #1.MosaicHunter 
            	# Check if both MotherID and FatherID are present
            	if [[ -n "$MotherID" && -n "$FatherID" ]]; then
            	sbatch $SCRIPTDIR/MosaicHunter_WES_Trio.sh -s $ProbandID -b $BamDIR -d $OUTDIR -g $ProbandGender -f $FatherID -m $MotherID -c $CONFIG_FILE
            	else
            	sbatch $SCRIPTDIR/MosaicHunter_WES_Singlemode.sh -s $ProbandID -b $BamDIR -d $OUTDIR -g $ProbandGender -c $CONFIG_FILE 
            	fi
		
		# For Parents (if available)
		# Check if either MotherID or FatherID is present
				if [[ -n "$MotherID" || -n "$FatherID" ]]; then
  					if [[ -n "$MotherID" ]]; then
    					sbatch "$SCRIPTDIR/MosaicHunter_WES_Singlemode.sh" -s "$MotherID" -b "$BamDIR" -d "$OUTDIR" -g "F" -c "$CONFIG_FILE"
  					fi

  					if [[ -n "$FatherID" ]]; then
    					sbatch "$SCRIPTDIR/MosaicHunter_WES_Singlemode.sh" -s "$FatherID" -b "$BamDIR" -d "$OUTDIR" -g "M" -c "$CONFIG_FILE"
  					fi
				fi
			
    #2.Mutect2

		#Check if the PON contains the sample in the family
       		for samples in "$ProbandID" "$MotherID" "$FatherID"; do 

		# Store the result of the grep command in a variable
			normalSample=$(bcftools view $PON | grep "$samples")

		# Check if $SampleID is present in the result
			if [ -n "$normalSample" ]; then
    			echo "$samples is present. No Mutect2 will be performed. Provide another Panel Of Normal." >> $OUTDIR/$ProbandID.pipeline.log
			else
			
			Mutect2="sbatch $SCRIPTDIR/Mutect2.singlemode.sh -b $BamDIR -s $samples -c $CONFIG_FILE -o $OUTDIR"

			# Submit the Mutect2 job
			Mutect2JobID=$($Mutect2 | awk '{print $NF}')

			# Submit the FilterMutect2 job with a dependency on Mutect2
			sbatch --export=ALL --dependency=afterok:${Mutect2JobID} $SCRIPTDIR/Mutect2.FilterMutect2.sh -s $samples -v $OUTDIR -c $CONFIG_FILE
			fi
   			

   #3.MosaicForecast

    			MF1="sbatch --export=ALL --dependency=afterok:${Mutect2JobID} $SCRIPTDIR/MF1_ProcessInput.sh -s $samples -b $BamDIR -o $OUTDIR -c $CONFIG_FILE"
    			MF1_job_id=$($MF1 | awk '{print $NF}')

    			MF2="sbatch --export=ALL --dependency=afterok:${MF1_job_id} $SCRIPTDIR/MF2_Extractreadlevel-singularity.sh -b $BamDIR -s $samples -c $CONFIG_FILE -o $OUTDIR"
    			MF2_job_id=$($MF2 | awk '{print $NF}')

    			MF3="sbatch --export=ALL --dependency=afterok:${MF2_job_id} $SCRIPTDIR/MF3.GenotypePredictionsl-singularity.sh -s $samples -c $CONFIG_FILE -o $OUTDIR"
    			MF3_job_id=$($MF3 | awk '{print $NF}')

		done

 #4. Germline variant calling- GATKHC (Cannot run in HPC yet with this script)
 
		for samples in "$ProbandID" "$MotherID" "$FatherID"; do

			 	if [ $CONFIG="hs37d5" ]; then

				CONFIG_for_GATKHC=$SCRIPTDIR/BWA-GATKHC.hs37d5_phoenix.cfg

				GATKHCjob=`sbatch --array=0-23 $SCRIPTDIR/GATK.HC_Universal_phoenix.sh -S $samples $BamDIR -c $CONFIG_for_GATKHC`
				GATKHCjob=$(echo $GATKHCjob | awk '{print $NF}')

				sbatch --export=ALL --dependency=afterok:${GATKHCjob} $SCRIPTDIR/GATK.gatherVCFs_Universal_phoenix.sh -c $CONFIG_for_GATKHC -S $samples -o $OUTDIR	

				else

				echo "Please provide alternate config for GATK-HC, that looks like 
					$SCRIPTDIR/BWA-GATKHC.TEMPLATE_phoenix.cfg" >> $OUTDIR/$ProbandID.pipeline.log
				fi

		done
	
done


#Ensure part2 is executed only after part1 completed

wait


#2. To compile all variant calls in single file according variant calls
## need to run this command only after the above ones
for samples in "$ProbandID" "$MotherID" "$FatherID"; do
awk -v ID="$samples" '$0 !~ /^##/ {print ID "\t" "Mut" "\t" $0}' $OUTDIR/$samples.mutect2.singlemode.PASS.aaf.vcf >> Mutect2.calls
awk -v ID="$samples" '$0 !~ /^##/ {print ID "\t" "MF" "\t" $0}' $OUTDIR/$samples.mosaicforecast.genotype.predictions.refined.bed >> MosaicForecast.calls
awk -v ID="$samples" '$0 !~ /^##/ {print ID "\t" "MH" "\t" $0}' $OUTDIR/$samples.final.passed.tsv >> MosaicHunter.calls.txt
done
