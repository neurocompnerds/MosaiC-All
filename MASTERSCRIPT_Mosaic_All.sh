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
# Usage $0 -s \$sampleID_list -o \$Output_folder -c \$Config_File | [ - h | --help ]
#
# Options
#-s REQUIRED sampleID.list (one header row and then tab-delimited columns \$BAMdir,\$ProbandID,\$Gender,\$Mother,\$Father)
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
#Test03: bash $SCRIPTDIR/MASTERSCRIPT_Mosaic_All.sh -s /hpcfs/users/a1742674/MosaiC-All/SampleID -o /hpcfs/users/a1742674/MosaiC-All/Outputs -c $SCRIPTDIR/MosaiC-All.config
#Test04: bash $SCRIPTDIR/MASTERSCRIPT_Mosaic_All.sh -s /hpcfs/users/a1742674/MosaiC-All/SampleID_All -o /hpcfs/users/a1742674/MosaiC-All/Output_160823 -c $SCRIPTDIR/MosaiC-All.config
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
                -h | --help )           usage
                                        exit 0
                                        ;;
                 * )                    usage
                                        exit 1
        esac
        shift
done

## Define Directories, Variables and 
# If the script lacks any requirements then fail immediately or create one

if [ -z "$SAMPLELIST" ]; then
    usage
    echo "## ERROR: You need to provide a sample list 
	#-s REQUIRED sampleID.list (one header row and then tab-delimited columns \$BAMdir,\$ProbandID,\$Gender,\$Mother,\$Father)"
	exit 1
fi

if [ -z "$CONFIG_FILE" ]; then
    usage
    echo "## ERROR: You need to provide a config file 
	#-c REQUIRED ConfigFile (i.e /hpcfs/groups/phoenix-hpc-neurogenetics/Nandini/Mosaic-All/Mosaic-S/Mosaic-All.config)"
	exit 1
fi

source $CONFIG_FILE

if [ ! -d "${OUTDIR}" ]; then
    mkdir -p ${OUTDIR}
	echo "## INFO: output directory created, you'll find all of the outputs and log files in here: ${OUTDIR}" >> $LOGDIR/$ProbandID.pipeline.log
fi

LOGDIR="/hpcfs/users/${USER}/MosaiC-All/Log"

if [ ! -d "$LOGDIR" ]; then
    mkdir -p $LOGDIR
	echo "## INFO: Slurm log files will be placed in this location $LOGDIR"
fi

#Array from a list of Samples (ignoring the header of the file)
mapfile -t SAMPLEID < <(tail -n +2 "$SAMPLELIST")

#modules
module purge
module load BCFtools/1.17-GCC-11.2.0

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
       	# Submit MHjob for Proband either in triomode or singlemode
	# so, need to Check if both MotherID and FatherID are present
     		if [[ -n "$MotherID" && -n "$FatherID" ]]; then
        	MH_P="sbatch $SCRIPTDIR/MosaicHunter_WES_Trio.sh -s $ProbandID -b $BamDIR -d $OUTDIR -g $ProbandGender -f $FatherID -m $MotherID -c $CONFIG_FILE"
	 	MH_P_JobID=$($MH_P | awk '{print $NF}')
		sbatch --export=ALL --dependency=afterok:${MH_P_JobID} $SCRIPTDIR/CombineCalls.sh -s $ProbandID -v MH -f $OUTDIR/$samples.final.passed.tsv -o $OUTDIR
 		 	if [ "$MH_P_JobID" ]; then
   			Count_MH=$(cat $OUTDIR/$samples.final.passed.tsv | wc -l)
			echo "$ProbandID\tTriomode\t$Count_MH" | tr " " "\t" >> Counts.txt
			fi
        	else
        	MH_P="sbatch $SCRIPTDIR/MosaicHunter_WES_Singlemode.sh -s $ProbandID -b $BamDIR -d $OUTDIR -g $ProbandGender -c $CONFIG_FILE"
	 	MH_P_JobID=$($MH_P | awk '{print $NF}')
		sbatch --export=ALL --dependency=afterok:${MH_P_JobID} $SCRIPTDIR/CombineCalls.sh -s $ProbandID -v MH -f $OUTDIR/$samples.final.passed.tsv -o $OUTDIR
        	fi
		
	# Submit MHjob for Parents
	# Check if either MotherID or FatherID is present
		if [[ -n "$MotherID" ]]; then
    		sbatch "$SCRIPTDIR/MosaicHunter_WES_Singlemode.sh" -s "$MotherID" -b "$BamDIR" -d "$OUTDIR" -g "F" -c "$CONFIG_FILE"
  		fi

  		if [[ -n "$FatherID" ]]; then
    	    	sbatch "$SCRIPTDIR/MosaicHunter_WES_Singlemode.sh" -s "$FatherID" -b "$BamDIR" -d "$OUTDIR" -g "M" -c "$CONFIG_FILE"
  		fi
			
    #2.Mutect2 and MosaicForecast

		#Check if the PON contains the sample in the family
       	for samples in "$ProbandID" "$MotherID" "$FatherID"; do 

		# Store the result of the grep command in a variable
			normalSample=$(bcftools view $PON_A | grep "$samples")

		# Check if $SampleID is present in the result
			if [ -n "$normalSample" ]; then
    			echo "$samples is present in $PON_A. So, let's do Mutect2 on this sample using $PON_B. But, we still have to check $PON_B first"
       			normalSample_B=$(bcftools view $PON_B | grep "$samples")
	  			#check sample in PON_B
	   			if [ -n "$normalSample_B" ]; then
       					#execute Mutect2 using PON_B
       					Mutect2="sbatch $SCRIPTDIR/Mutect2.singlemode.sh -b $BamDIR -s $samples -c $CONFIG_FILE -o $OUTDIR -p $PON_B"
                			Mutect2JobID=$($Mutect2 | awk '{print $NF}')
		   			#execute FilterMutect2 which depends on Mutect2
		  			sbatch --export=ALL --dependency=afterok:${Mutect2JobID} $SCRIPTDIR/Mutect2.FilterMutect2.sh -s $samples -v $OUTDIR -c $CONFIG_FILE
       					#execute MosaicForecast (step1) which depends on Mutect2
       					MF1="sbatch --export=ALL --dependency=afterok:${Mutect2JobID} $SCRIPTDIR/MF1_ProcessInput.sh -s $samples -b $BamDIR -o $OUTDIR -c $CONFIG_FILE"
    					MF1_job_id=$($MF1 | awk '{print $NF}')
					#execute MosaicForecast (step2) which depends on MosaicForecast (step1)
    					MF2="sbatch --export=ALL --dependency=afterok:${MF1_job_id} $SCRIPTDIR/MF2_Extractreadlevel-singularity.sh -b $BamDIR -s $samples -c $CONFIG_FILE -o $OUTDIR"
    					MF2_job_id=$($MF2 | awk '{print $NF}')
					#execute MosaicForecast (step3) which depends on MosaicForecast (step2)
   					MF3="sbatch --export=ALL --dependency=afterok:${MF2_job_id} $SCRIPTDIR/MF3.GenotypePredictionsl-singularity.sh -s $samples -c $CONFIG_FILE -o $OUTDIR"
   					MF3_job_id=$($MF3 | awk '{print $NF}')
      				
		  		else
      					echo "$samples present in both PON, provide another one. No Mutect2 is performed for this sample: MUTECT2 FAIL"
	  			fi
			else
			    	# Submit the Mutect2 job using PON_A
				echo "$samples is not present in $PON_A. So, let's do Mutect2 on this sample using $PON_A"
				Mutect2="sbatch $SCRIPTDIR/Mutect2.singlemode.sh -b $BamDIR -s $samples -c $CONFIG_FILE -o $OUTDIR -p $PON_A"
                		Mutect2JobID=$($Mutect2 | awk '{print $NF}')
                		#execute FilterMutect2 which depends on Mutect2
		  		sbatch --export=ALL --dependency=afterok:${Mutect2JobID} $SCRIPTDIR/Mutect2.FilterMutect2.sh -s $samples -v $OUTDIR -c $CONFIG_FILE
       				#execute MosaicForecast (step1) which depends on Mutect2
       				MF1="sbatch --export=ALL --dependency=afterok:${Mutect2JobID} $SCRIPTDIR/MF1_ProcessInput.sh -s $samples -b $BamDIR -o $OUTDIR -c $CONFIG_FILE"
    				MF1_job_id=$($MF1 | awk '{print $NF}')
				#execute MosaicForecast (step2) which depends on MosaicForecast (step1)
    				MF2="sbatch --export=ALL --dependency=afterok:${MF1_job_id} $SCRIPTDIR/MF2_Extractreadlevel-singularity.sh -b $BamDIR -s $samples -c $CONFIG_FILE -o $OUTDIR"
    				MF2_job_id=$($MF2 | awk '{print $NF}')
				#execute MosaicForecast (step3) which depends on MosaicForecast (step2)
   				MF3="sbatch --export=ALL --dependency=afterok:${MF2_job_id} $SCRIPTDIR/MF3.GenotypePredictionsl-singularity.sh -s $samples -c $CONFIG_FILE -o $OUTDIR"
   				MF3_job_id=$($MF3 | awk '{print $NF}')
			    	sbatch --export=ALL --dependency=afterok:${Mutect2JobID} $SCRIPTDIR/Mutect2.FilterMutect2.sh -s $samples -v $OUTDIR -c $CONFIG_FILE
			fi
			
	done

    #4. Germline variant calling- GATKHC (Cannot run in HPC yet with this script)
 
	for samples in "$ProbandID" "$MotherID" "$FatherID"; do # The GATK.HC config is now defined in the Mosaic-All.config

		GATKHCjob=`sbatch --array=0-23 $SCRIPTDIR/GATK.HC_Universal_phoenix.sh -S $samples -o $BamDIR -c $SCRIPTDIR/$CONFIG_for_GATKHC`
		GATKHCjob=$(echo $GATKHCjob | awk '{print $NF}')
		sbatch --export=ALL --dependency=afterok:${GATKHCjob} $SCRIPTDIR/GATK.gatherVCFs_Universal_phoenix.sh -c $SCRIPTDIR/$CONFIG_for_GATKHC -S $samples -o $OUTDIR	

	done
	
done


