#!/bin/bash
# MasterScript: Variant calling step 
# 1. Coverage Analysis of every bam file
# 2. Mutect2 and FilterMutect2: Parents,Probands and Siblings (if available)
# 3. MosaicHunter: Parents,Probands and Siblings (if available)
# 4. MosaicForecast on Mutect2 variant callset, followed by Filter
# UniAdelaide-HPC friendly
# Date: 9th June 2023
# 

usage()
{
echo "#MasterScript: Variant calling steps, which includes
# 1. Coverage Analysis of every bam file
# 2. Mutect2: Parents and Probands and Siblings (if available)
# 3. MosaicHunter: Parents and Probands and Siblings (if available)
# 4. MosaicForecast on Mutect2 variant callset
#
#
# Usage $0 -s /path/to/sampleID.list -o /path/to/output_folder -c /path/to/config_file | [ - h | --help ]
#
# Options
#-s <file>      REQUIRED: A file e.g. sampleID.list (one header row and then tab-delimited columns \$BAMdir,\$ProbandID,\$Gender,\$Mother,\$Father)
#-o <directory> REQUIRED: Output directory (all variant calls for all samples will output into a single directory)
#-c <file>      REQUIRED: Configuration File (This file sets paths and defaults relevant to your system e.g. see: config/Mosaic-All.config)
#
# -h or --help  Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
#
# Original: Nandini Sandran, 9/6/2023
# Modified: (Date; Name; Description)
# See: https://github.com/ngs1810/MosaiC-All for history and new updates.
#
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
            *  )                    usage
                                    exit 1
    esac
    shift
done

## Define Directories, Variables and 
# If the script lacks any requirements then fail immediately or create one

if [ -z "$SAMPLELIST" ]; then
    usage
    echo "## ERROR: You need to provide a sample list" 
    echo "#-s REQUIRED sampleID.list (one header row and then tab-delimited columns \$BAMdir,\$ProbandID,\$Gender,\$Mother,\$Father)"
	exit 1
fi

if [ -z "$CONFIG_FILE" ]; then
    usage
    echo "## ERROR: You need to provide a config file. Check the config/Mosaic-All.config file for an example."
    echo "#-c <file>  REQUIRED: Configuration File (This file sets paths and defaults relevant to your system e.g. see: config/Mosaic-All.config)"
	exit 1
fi

source $CONFIG_FILE

if [ ! -d "$LOGDIR" ]; then
    mkdir -p $LOGDIR
	echo "## INFO: Slurm log files will be placed in this location $LOGDIR"
fi

if [ ! -d "${OUTDIR}" ]; then
    mkdir -p ${OUTDIR}
	echo "## INFO: output directory created, you'll find all of the outputs and log files in here: ${OUTDIR}" >> $OUTDIR/Mosaic-All.pipeline.log
fi

#Array from a list of Samples (ignoring the header of the file)
mapfile -t SAMPLEID < <(tail -n +2 "$SAMPLELIST")

#modules
module purge
module load BCFtools/1.17-GCC-11.2.0

# Iteration for variant calling starts here
for SAMPLEID in "${SAMPLEID[@]}"; do

    #Defining variables from each row
	BAMDIR=$(awk '{print $1}' <<< "$SAMPLEID ")
    ProbandID=$(awk '{print $2}' <<< "$SAMPLEID ")
   	Gender=$(awk '{print $3}' <<< "$SAMPLEID ")
   	MotherID=$(awk '{print $4}' <<< "$SAMPLEID ")
   	FatherID=$(awk '{print $5}' <<< "$SAMPLEID ")

	echo "Pipeline for $ProbandID, $MotherID, $FatherID in $BAMDIR" >> $OUTDIR/$ProbandID.pipeline.log

    #1.MosaicHunter 
    # Submit MHjob for Proband either in triomode or singlemode
	# so, need to Check if both MotherID and FatherID are present
    if [[ -n "$MotherID" && -n "$FatherID" ]]; then
        sbatch $SCRIPTDIR/scripts/MosaicHunter_WES_Trio.sh -s $ProbandID -b $BAMDIR -d $OUTDIR -g $Gender -f $FatherID -m $MotherID -c $CONFIG_FILE
    else
	echo "## MH INFO: No parents provided, so, the proband is assumed singleton, thus, singlemodeMH will be used" >> $OUTDIR/$SampleID.pipeline.log
        sbatch $SCRIPTDIR/scripts/MosaicHunter_WES_Singlemode.sh -s $ProbandID -b $BAMDIR -d $OUTDIR -g $Gender -c $CONFIG_FILE
    fi
		
	# Submit MHjob for Parents
	# Check if either MotherID or FatherID is present
	if [[ -n "$MotherID" ]]; then
    	sbatch "$SCRIPTDIR/scripts/MosaicHunter_WES_Singlemode.sh" -s "$MotherID" -b "$BAMDIR" -d "$OUTDIR" -g "F" -c "$CONFIG_FILE"
	else
	echo "## INFO: MotherID is empty." >> "$OUTDIR/$SampleID.pipeline.log"
  	fi

  	if [[ -n "$FatherID" ]]; then
       	sbatch "$SCRIPTDIR/scripts/MosaicHunter_WES_Singlemode.sh" -s "$FatherID" -b "$BAMDIR" -d "$OUTDIR" -g "M" -c "$CONFIG_FILE"
  	else
	echo "## INFO: FatherID is empty." >> "$OUTDIR/$SampleID.pipeline.log"
	fi
			
    #2.Mutect2 and MosaicForecast

	#Check if the PON contains the sample in the family
    	for SampleID in "$ProbandID" "$MotherID" "$FatherID"; do 

	if [[ -n "$SampleID" ]]; then

	# Store the result of the grep command in a variable
		normalSample=$(bcftools view $PON_A |  grep -oE -- "$SampleID")

		# Check if $SampleID is present in the result
		if [ -n "$normalSample" ]; then
    		echo "## Mutect2 WARN: $SampleID is present in $PON_A. Checking for this sample in $PON_B." >> $OUTDIR/$SampleID.pipeline.log
    		normalSample_B=$(bcftools view $PON_B | grep -oE -- "$SampleID")
	    	#check sample in PON_B
	   		if [ -z "$normalSample_B" ]; then
                	echo "## Mutect2 INFO: $SampleID is not present in $PON_B. So, let's do Mutect2 on this sample using $PON_B" >> $OUTDIR/$SampleID.pipeline.log
       			PON=$PON_B
		  	else
   			echo "## WARN: $SampleID present in both Panel of Normal VCFs, you will need to provide another one and set this in $CONFIG_FILE. Mutect2 was not performed for this sample." >> $OUTDIR/$SampleID.pipeline.log
	  		fi
		else
		# Submit the Mutect2 job using PON_A
		echo "## Mutect2 INFO: $SampleID is not present in $PON_A. So, let's do Mutect2 on this sample using $PON_A" >> $OUTDIR/$SampleID.pipeline.log
            	PON=$PON_A
		fi

	# Run Mutect2 and MosaicForecast using the selected PON
        Mutect2="sbatch $SCRIPTDIR/scripts/Mutect2.singlemode.sh -b $BAMDIR -s $SampleID -c $CONFIG_FILE -o $OUTDIR -p $PON"
        Mutect2JobID=$($Mutect2 | awk '{print $NF}')

        #execute FilterMutect2 which depends on Mutect2
        sbatch --export=ALL --dependency=afterok:${Mutect2JobID} $SCRIPTDIR/scripts/Mutect2.FilterMutect2.sh -s $SampleID -v $OUTDIR -c $CONFIG_FILE

        #execute MosaicForecast (step1) which depends on Mutect2
        MF1="sbatch --export=ALL --dependency=afterok:${Mutect2JobID} $SCRIPTDIR/scripts/MF1_ProcessInput.sh -s $SampleID -b $BAMDIR -o $OUTDIR -c $CONFIG_FILE"
        MF1_job_id=$($MF1 | awk '{print $NF}')

        #execute MosaicForecast (step2) which depends on MosaicForecast (step1)
        MF2="sbatch --export=ALL --dependency=afterok:${MF1_job_id} $SCRIPTDIR/scripts/MF2_Extractreadlevel-singularity.sh -b $BAMDIR -s $SampleID -c $CONFIG_FILE -o $OUTDIR"
        MF2_job_id=$($MF2 | awk '{print $NF}')

        #execute MosaicForecast (step3) which depends on MosaicForecast (step2)
        MF3="sbatch --export=ALL --dependency=afterok:${MF2_job_id} $SCRIPTDIR/scripts/MF3.GenotypePredictions-singularity.sh -s $SampleID -c $CONFIG_FILE -o $OUTDIR"
        MF3_job_id=$($MF3 | awk '{print $NF}')

	
	else
	echo "##Mutect2 and MF INFO: MotherID or FatherID is empty. Skipping the loop." >> "$OUTDIR/$SampleID.pipeline.log"
	fi
	
	done
	
done
