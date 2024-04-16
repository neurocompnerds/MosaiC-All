#!/bin/sh
#SBATCH -J MosaicHunter_Single.sh
#SBATCH -o /hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/Mosaic/MosaiC-All/TestRun/MH_Singleton-slurm-%j.out

#SBATCH -p skylake
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=8:00:00
#SBATCH --mem=100GB
#SBATCH --gres=tmpfs:40G
# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=%u@adelaide.edu.au

usage()
{
echo "
#
#README
#This script detects variants in exome data in single mode
#There are three major steps (Prefilter-> MH variant caller-> Process output files -> Log files)
#
#DESIGNED for HPC, which needs to be executed with MasterScript only
#
#CHANGEME please
#Line 3 Specify slurm-log directory
#
#SCRIPT EXECUTE
#sbatch $0 -s \$SampleID -b /path/to/bam-file/ -d /path/to/output -g M -c $SCRIPTDIR/config 
#
#-s REQUIRED SampleID                         (e.g A unique identifier for your sample that is also in the BAM file name)
#-b REQUIRED Directory of Bam files           (/path/to/bam_file/)
#-d REQUIRED Directory of Outputs             (/path/to/output/)
#-g REQUIRED Biological sex of of the Sample  (M or F)
#-c REQUIRED The script configuration file    (/path/to/config_file.cfg)
#
#OUTPUT lists
#1.\$OUTDIR/\${SampleID[$SLURM_ARRAY_TASK_ID]}                            -> Two folders containing outputs from MH execution that are removed later to save on storage
#2.\$OUTDIR/\${SampleID[$SLURM_ARRAY_TASK_ID]}.final.passed.tsv           -> Raw variants output file
#3.\$OUTDIR/\${SampleID[$SLURM_ARRAY_TASK_ID]}.forAnnovar.singlemode.vcf  -> List of variants with useful info only
#4.\$OUTDIR/\${SampleID[$SLURM_ARRAY_TASK_ID]}.log                        -> StdOutput.lof file from MH
#5.\$OUTDIR/\${SampleID[$SLURM_ARRAY_TASK_ID]}.summary.log                -> Logfile produced by MH
"
}

## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -s )                    shift
                                        SampleID=$1
                                        ;;
                -b )                    shift
                                        BAMDIR=$1
                                        ;;
                -d )                    shift
                                        OUTDIR=$1
                                        ;;
                -g )                    shift
                                        Gender=$1
                                        ;;
                -c )                    shift
                                        CONFIG_FILE=$1
                                        ;;
                * )                     usage
                                        exit 1
        esac
        shift
done

## Check for all required arguments
for ARG in $BAMDIR $OUTDIR $CONFIG_FILE $SampleID $Gender; do
    if [ -z "${ARG}" ]; then
        usage
        echo "## ERROR: Please check that you have supplied all required argunments for this script"
        exit 1
    fi
done


#module for HPCs
module purge
module use /apps/modules/all
module load Java/1.8.0_121
module load BLAT/3.5-foss-2016b

#specify variables and directories based on config files

ProbandBamFile=$(find "$BAMDIR" -type f -name "$SampleID*.bam")
source $CONFIG_FILE

#1.prefilter

java -jar $MHDIR/build/mosaichunter.jar -C $MHDIR/conf/$CONFIG/exome_parameters.properties \
-P reference_file=$REFGEN \
-P input_file=$ProbandBamFile \
-P heterozygous_filter.sex=$Gender \
-P output_dir=$OUTDIR/$SampleID.parameters.log

echo "Pre-filter completed for $SampleID"

#2.define Alpha and beta value and remove white spaces before

Al=$(cat $OUTDIR/$SampleID.parameters.log/stdout*.log | grep "alpha" | cut -d ":" -f2)
Alpha=$(echo "$Al" | sed 's/^ *//g')

Be=$(cat $OUTDIR/$SampleID.parameters.log/stdout*.log | grep "beta" | cut -d ":" -f2)
Beta=$(echo "$Be" | sed 's/^ *//g')

Dp=$(cat $OUTDIR/$SampleID.parameters.log/stdout*.log | grep "average depth" | cut -d ":" -f2)
Depth=$(echo "$Dp" | sed 's/^ *//g')

#3.execute mosaic variant calling

java -Djava.io.tmpdir=${TMPDIR} -jar $MHDIR/build/mosaichunter.jar -C $MHDIR/conf/$CONFIG/exome.properties \
-P reference_file=$REFGEN \
-P input_file=$ProbandBamFile \
-P mosaic_filter.sex=$Gender \
-P depth=$Depth \
-P mosaic_filter.alpha_param=$Alpha \
-P mosaic_filter.beta_param=$Beta \
-P mosaic_filter.dbsnp_file=$DBSNP \
-P repetitive_region_filter.bed_file=$REPEATS \
-P common_site_filter.bed_file=$COMMONERROR \
-P output_dir=$OUTDIR/$SampleID

echo "## TNFO: (MosaicHunter_WES_Singlemode.sh) Somatic variant calling completed for $SampleID" >> $OUTDIR/$SampleID.pipeline.log

#5.Process the outputs files

cat $OUTDIR/$SampleID/final.passed.tsv  > $OUTDIR/$SampleID.final.passed.tsv
cat $OUTDIR/$SampleID/final.passed.tsv | awk '{print $1, $2, $7, $9}' | tr " " "\t" > $OUTDIR/$SampleID.MH.forAnnovar.singlemode.vcf

#6 log file

grep "input_file =" $OUTDIR/$SampleID/stdout_*.log > $OUTDIR/$SampleID.summary.log
tail  $OUTDIR/$SampleID/stdout_*.log -n16 >> $OUTDIR/$SampleID.summary.log
cat $OUTDIR/$SampleID/stdout_*.log >> $OUTDIR/$SampleID.stdout.log

echo "## TNFO: (MosaicHunter_WES_Singlemode.sh) done for $SampleID" >> $OUTDIR/$SampleID.pipeline.log

#Remove Folders
rm -r $OUTDIR/$SampleID.parameters.log
rm -r $OUTDIR/$SampleID
