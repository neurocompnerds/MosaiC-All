#!/bin/sh
#SBATCH -J MosaicHunter_Trio.sh
#SBATCH -o /hpcfs/users/%u/Mosaic-All/Log/MH_Trio-slurm-%j.out
#SBATCH -A robinson
#SBATCH -p skylake,icelake
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=8:00:00
#SBATCH --mem=60GB
# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nandini.sandran@adelaide.edu.au

usage()
{
echo "
#
#README
#This script detects variants in exome data in trio mode
#There are three major steps (Prefilter-> MH variant caller-> Process output files -> Log files)
#
#DESIGNED for HPC, which needs to be executed with MasterScript only, or else,
#for standalone-script, find in another place (need to update the folder, latest one in GA_Scripts)
#
#CHANGEME please
#Line 3 Specify slurm-log directory
#
#SCRIPT EXECUTE
#
#sbatch $SCRIPTDIR/MosaicHunter_WES_Trio.sh -s \$Proband -b \$BamDIR -d \$OUTDIR -g M -f \$FATHERID -m \$MOTHERID -c \$CONFIG_FILE
#all of this will be "specified" in MasterScript
#
#-s REQUIRED SampleID   			(e.g 004P, no suffix)
#-b REQUIRED Directory of Bam files 		(e.g $HOME/SSC/ProSib/Bam)
#-d REQUIRED Directory of Outputs   		(MH gives a folder as output, only speficy the directory intended)
#-g REQUIRED Gender of Proband 	  		(M or F)
#-f REQUIRED FatherID 	  			(004F)
#-m REQUIRED MotherID		 	  	(004M)
#-c REQUIRED Config_File     			(e.g /hpcfs/groups/phoenix-hpc-neurogenetics/Nandini/Mosaic-All/Mosaic-S/Mosaic-All.config) 
#
#OUTPUT lists
#1.\$OUTDIR/\${sample[$SLURM_ARRAY_TASK_ID]} 					-> Two Major Folders consisting MHexecution temp files and final.passed.tsv
#2.\$OUTDIR/\${sample[$SLURM_ARRAY_TASK_ID]}.final.passed.tsv			-> Raw variants output file
#3.\$OUTDIR/\${sample[$SLURM_ARRAY_TASK_ID]}.forAnnovar.triomode.vcf		-> List of variants with useful info only
#4.\$OUTDIR/\${sample[$SLURM_ARRAY_TASK_ID]}.log				  	-> StdOutput.lof file from MH
#5.\$OUTDIR/\${sample[$SLURM_ARRAY_TASK_ID]}.summary.log				-> Logfile produced by MH
#
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
                -f )                    shift
                                        FATHER=$1
                                        ;;
                -m )                    shift
                                        MOTHER=$1
                                        ;;
                -c )                    shift
                                        CONFIG_FILE=$1
                                        ;;
                * )                     usage
                                        exit 1
        esac
        shift
done

## Note: There are no tests for any required arguments

#define variables and directory for MosaicHunter
ProbandBamFile=$(/usr/bin/find "$BAMDIR" -type f -name "$SampleID.*.bam")
MotherBamFile=$(/usr/bin/find "$BAMDIR" -type f -name "$MOTHER.*.bam")
FatherBamFile=$(/usr/bin/find "$BAMDIR" -type f -name "$FATHER.*.bam")

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
-P output_dir=$OUTDIR/$SampleID

echo "Somatic variant calling completed for $SampleID" 

#5.Process the outputs files

cat $OUTDIR/$SampleID/final.passed.tsv  > $OUTDIR/$SampleID.final.passed.tsv
cat $OUTDIR/$SampleID/final.passed.tsv | awk '{print $1, $2, $7, $9}' | tr " " "\t" > $OUTDIR/$SampleID.MH.forAnnovar.triomode.vcf
#cat $MHDIR/Head.vcf > $OUTDIR/$SampleID.final.triomode.vcf
#cat $OUTDIR/$SampleID.forAnnovar.vcf | awk 'BEGIN{FS="\t"} {print $1, $2, ".", $3, $4, "50", "PASS", "NS=18;DP=10"}' | tr " " "\t" > $OUTDIR/$SampleID.final.triomode.vcf

#6 log file
grep "input_file =" $OUTDIR/$SampleID/stdout_*.log > $OUTDIR/$SampleID.summary.log
tail  $OUTDIR/$SampleID/stdout_*.log -n16 >> $OUTDIR/$SampleID.summary.log
cat $OUTDIR/$SampleID/stdout_*.log >> $OUTDIR/$SampleID.stdout.log

echo "done for $SampleID"

#Remove Folders as this takes lots of memory
rm -r $OUTDIR/$SampleID.parameters.log
rm -r $OUTDIR/$SampleID
