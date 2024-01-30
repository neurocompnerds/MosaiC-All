#!/bin/sh
#SBATCH -J MosaicHunter_Trio.sh
#SBATCH -o /home/%u/Mosaic-All/Log/MH_Trio-slurm-%j.out

#SBATCH -p skylake,icelake
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=8:00:00
#SBATCH --mem=60GB
# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=%u@adelaide.edu.au

usage()
{
echo "
#
#README
#This script detects variants in exome data in trio mode
#There are three major steps (Prefilter-> MH variant caller-> Process output files -> Log files)
#
#DESIGNED for HPC, which needs to be executed with MasterScript only
#
#SCRIPT EXECUTE
#sbatch $0 -s \$Proband -b \$BamDIR -d \$OUTDIR -g M -f \$FatherID -m \$MotherID -c \$CONFIG_FILE
#all of this will be "specified" in MasterScript
#
#-s <string>    REQUIRED ProbandID                      (e.g A unique identifier for your sample that is also in the BAM file name)
#-b <directory> REQUIRED Directory of Bam files         (e.g /path/to/bam_files/)
#-d <directory> REQUIRED Directory of Outputs           (MH gives a folder as output, only speficy the directory intended)
#-g <string>    REQUIRED Sex of Proband                 (M or F)
#-f <string>    REQUIRED FatherID                       (A unique identifier for the father's sample that is also in the BAM file name)
#-m <string>    REQUIRED MotherID                       (A unique identifier for the mother's sample that is also in the BAM file name)
#-c <file>      REQUIRED The script configuration file  (/path/to/config_file.cfg)
#
#OUTPUT lists
#1.\$OUTDIR/\${ProbandID[$SLURM_ARRAY_TASK_ID]}                          -> Two Major Folders consisting MHexecution temp files and final.passed.tsv
#2.\$OUTDIR/\${ProbandID[$SLURM_ARRAY_TASK_ID]}.final.passed.tsv         -> Raw variants output file
#3.\$OUTDIR/\${ProbandID[$SLURM_ARRAY_TASK_ID]}.forAnnovar.triomode.vcf  -> List of variants with useful info only
#4.\$OUTDIR/\${ProbandID[$SLURM_ARRAY_TASK_ID]}.log                      -> StdOutput.lof file from MH
#5.\$OUTDIR/\${ProbandID[$SLURM_ARRAY_TASK_ID]}.summary.log              -> Logfile produced by MH
#
"
}

## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -s )                    shift
                                        ProbandID=$1
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
                -f )                    shift
                                        FatherID=$1
                                        ;;
                -m )                    shift
                                        MotherID=$1
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
for ARG in $BAMDIR $OUTDIR $CONFIG_FILE $ProbandID $Gender $MotherID $FatherID; do
    if [ -z "${ARG}" ]; then
        usage
        echo "## ERROR: Please check that you have supplied all required argunments for this script"
        exit 1
    fi
done

#define variables and directory for MosaicHunter
ProbandBamFile=$(find "$BAMDIR" -type f -name "$ProbandID.*.bam")
MotherBamFile=$(find "$BAMDIR" -type f -name "$MotherID.*.bam")
FatherBamFile=$(find "$BAMDIR" -type f -name "$FatherID.*.bam")

source $CONFIG_FILE

#module for HPCs
module purge
module use /apps/modules/all
module load Java/1.8.0_121
module load BLAT/3.5-foss-2016b

#1.prefilter

java -jar $MHDIR/build/mosaichunter.jar -C $MHDIR/conf/exome_parameters.properties \
-P reference_file=$REFGEN \
-P input_file=$ProbandBamFile \
-P heterozygous_filter.sex=$Gender \
-P output_dir=$OUTDIR/$ProbandID.parameters.log

echo "## INFO: (MosaicHunter_WES_Trio.sh) Pre-filter completed for $ProbandID" >> $OUTDIR/$ProbandID.pipeline.log

#2.define Alpha and beta value and remove white spaces before

Al=$(cat $OUTDIR/$ProbandID.parameters.log/stdout*.log | grep "alpha" | cut -d ":" -f2)
Alpha=$(echo "$Al" | sed 's/^ *//g')

Be=$(cat $OUTDIR/$ProbandID.parameters.log/stdout*.log | grep "beta" | cut -d ":" -f2)
Beta=$(echo "$Be" | sed 's/^ *//g')

Dp=$(cat $OUTDIR/$ProbandID.parameters.log/stdout*.log | grep "average depth" | cut -d ":" -f2)
Depth=$(echo "$Dp" | sed 's/^ *//g')

#3.execute mosaic variant calling

java -jar $MHDIR/build/mosaichunter.jar -C $MHDIR/conf/exome.properties \
-P reference_file=$REFGEN \
-P input_file=$ProbandBamFile \
-P depth=$Depth \
-P mosaic_filter.father_bam_file=$MotherBamFile \
-P mosaic_filter.mother_bam_file=$FatherBamFile \
-P mosaic_filter.sex=$Gender \
-P mosaic_filter.alpha_param=$Alpha \
-P mosaic_filter.beta_param=$Beta \
-P mosaic_filter.mode=trio \
-P mosaic_filter.dbsnp_file=$DBSNP \
-P repetitive_region_filter.bed_file=$REPEATS \
-P common_site_filter.bed_file=$COMMONERROR \
-P output_dir=$OUTDIR/$ProbandID

echo "## INFO: (MosaicHunter_WES_Trio.sh) Somatic variant calling completed for $ProbandID" >> $OUTDIR/$ProbandID.pipeline.log

#5.Process the outputs files

cp $OUTDIR/$ProbandID/final.passed.tsv $OUTDIR/$ProbandID.final.passed.tsv
awk '{print $1, $2, $7, $9}' $OUTDIR/$ProbandID/final.passed.tsv | tr " " "\t" > $OUTDIR/$ProbandID.MH.forAnnovar.triomode.vcf
#cat $MHDIR/Head.vcf > $OUTDIR/$ProbandID.final.triomode.vcf
#cat $OUTDIR/$ProbandID.forAnnovar.vcf | awk 'BEGIN{FS="\t"} {print $1, $2, ".", $3, $4, "50", "PASS", "NS=18;DP=10"}' | tr " " "\t" > $OUTDIR/$ProbandID.final.triomode.vcf

#6 log file
grep "input_file =" $OUTDIR/$ProbandID/stdout_*.log > $OUTDIR/$ProbandID.summary.log
tail  $OUTDIR/$ProbandID/stdout_*.log -n16 >> $OUTDIR/$ProbandID.summary.log
cat $OUTDIR/$ProbandID/stdout_*.log >> $OUTDIR/$ProbandID.stdout.log

echo "## INFO: (MosaicHunter_WES_Trio.sh) done for $ProbandID" >> $OUTDIR/$ProbandID.pipeline.log

#Remove Folders as this takes lots of memory
rm -r $OUTDIR/$ProbandID.parameters.log
rm -r $OUTDIR/$ProbandID
