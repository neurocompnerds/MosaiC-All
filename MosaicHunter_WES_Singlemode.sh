#!/bin/sh
#SBATCH -J MosaicHunter_Single.sh
#SBATCH -o /hpcfs/users/%u/Mosaic-All/Log/MH_Singleton-slurm-%j.out
#SBATCH -A robinson
#SBATCH -p skylake
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=8:00:00
#SBATCH --mem=100GB
#SBATCH --gres=tmpfs:40G
# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nandini.sandran@adelaide.edu.au

usage()
{
echo "
#
#README
#This script detects variants in exome data in single mode
#Three major step (Prefilter-> MH variant caller-> Process output files -> Log files)
#
#DESIGNED for HPC, which needs to be executed with MasterScript only, or else,
#for standalone-script, find in another place (need to update the folder, latest one in GA_Scripts)
#
#CHANGEME please
#Line 3 Specify slurm-log directory
#
#SCRIPT EXECUTE
#
#sbatch $SCRIPTDIR/MosaicHunter_WES_Singlemode.sh -s $SAMPLE -b $BAMDIR -d $OUTDIR -g M -c $CONFIG_FILE 
#all of this will be "specified" in MasterScript
#
#-s REQUIRED SampleID   					(e.g 004P, or 004M, or 004F no suffix)
#-b REQUIRED Directory of Bam files 		(e.g /hpcfs/groups/phoenix-hpc-sacgf/scratch/ali/GA_bams)
#-d REQUIRED Directory of Outputs   		(/gpfs/users/a1149120/MosaicHunter_single)
#-g REQUIRED Gender of Sample 	  				(M or F)
#-c REQUIRED Config_File     				(e.g $SCRIPT.DIR/CONFIG_FILE) 
#
#OUTPUT lists
#1.$DIR/${sample[$SLURM_ARRAY_TASK_ID]} 						-> Two Folders consisting MHexecution but will removed due to storage
#2.$DIR/${sample[$SLURM_ARRAY_TASK_ID]}.final.passed.tsv			-> Raw variants output file
#3.$DIR/${sample[$SLURM_ARRAY_TASK_ID]}.forAnnovar.singlemode.vcf	-> List of variants with useful info only
#4.$DIR/${sample[$SLURM_ARRAY_TASK_ID]}.log				  		-> StdOutput.lof file from MH
#5.$DIR/${sample[$SLURM_ARRAY_TASK_ID]}.summary.log				-> Logfile produced by MH
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
                                        DIR=$1
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


#module for HPCs
module purge
module use /apps/modules/all
module load Java/1.8.0_121
module load BLAT/3.5-foss-2016b

#specify variables and directories based on config files

ProbandBamFile=$(/usr/bin/find "$BAMDIR" -type f -name "$SampleID.*.bam")
source $CONFIG_FILE



#1.prefilter

java -jar $MHDIR/build/mosaichunter.jar -C $MHDIR/conf/exome_parameters.properties \
-P reference_file=$REFGEN \
-P input_file=$ProbandBamFile \
-P heterozygous_filter.sex=$Gender \
-P output_dir=$DIR/$SampleID.parameters.log

echo "Pre-filter completed for $SampleID"

#2.define Alpha and beta value and remove white spaces before

Al=$(cat $DIR/$SampleID.parameters.log/stdout*.log | grep "alpha" | cut -d ":" -f2)
Alpha=$(echo "$Al" | sed 's/^ *//g')

Be=$(cat $DIR/$SampleID.parameters.log/stdout*.log | grep "beta" | cut -d ":" -f2)
Beta=$(echo "$Be" | sed 's/^ *//g')

Dp=$(cat $DIR/$SampleID.parameters.log/stdout*.log | grep "average depth" | cut -d ":" -f2)
Depth=$(echo "$Dp" | sed 's/^ *//g')

#3.execute mosaic variant calling

java -Djava.io.tmpdir=${TMPDIR} -jar $MHDIR/build/mosaichunter.jar -C $MHDIR/conf/exome.properties \
-P reference_file=$REFGEN \
-P input_file=$ProbandBamFile \
-P mosaic_filter.sex=$Gender \
-P depth=$Depth \
-P mosaic_filter.alpha_param=$Alpha \
-P mosaic_filter.beta_param=$Beta \
-P mosaic_filter.dbsnp_file=$DBSNP \
-P repetitive_region_filter.bed_file=$REPEATS \
-P common_site_filter.bed_file=$COMMONERROR \
-P output_dir=$DIR/$SampleID

echo "Somatic variant calling completed for $SampleID" 

#5.Process the outputs files

cat $DIR/$SampleID/final.passed.tsv  > $DIR/$SampleID.final.passed.tsv
cat $DIR/$SampleID/final.passed.tsv | awk '{print $1, $2, $7, $9}' | tr " " "\t" > $DIR/$SampleID.MH.forAnnovar.singlemode.vcf

#6 log file

grep "input_file =" $DIR/$SampleID/stdout_*.log > $DIR/$SampleID.summary.log
tail  $DIR/$SampleID/stdout_*.log -n16 >> $DIR/$SampleID.summary.log
cat $DIR/$SampleID/stdout_*.log >> $DIR/$SampleID.stdout.log

echo "done for $SampleID"

#Remove Folders
rm -r $DIR/$SampleID.parameters.log
rm -r $DIR/$SampleID
