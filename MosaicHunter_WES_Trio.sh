#!/bin/sh
#SBATCH -J MosaicHunter_Trio.sh
#SBATCH -o /hpcfs/groups/phoenix-hpc-neurogenetics/Nandini/Mosaic-All/Log/MH_Trio-slurm-%j.out
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
#sbatch $SCRIPTDIR/MosaicHunter_WES_Trio.sh -s $Proband -b $BamDIR -d $OUTDIR -g M -f $FATHERID -m $MOTHERID -r $REFGENOME 
#all of this will be "specified" in MasterScript
#
#-s REQUIRED SampleID   				(e.g 004P, no suffix)
#-b REQUIRED Directory of Bam files 		(e.g $HOME/SSC/ProSib/Bam)
#-d REQUIRED Directory of Outputs   		(MH gives a folder as output, only speficy the directory intended)
#-g REQUIRED Gender of Proband 	  		(M or F)
#-f REQUIRED FatherID 	  				(004F)
#-m REQUIRED MotherID		 	  		(004M)
#-r REQUIRED Directory/Reference Genome     	(e.g /hpcfs/groups/phoenix-hpc-sacgf/reference/hs37d5/hs37d5.fa) 
#
#OUTPUT lists
#1.$DIR/${sample[$SLURM_ARRAY_TASK_ID]} 						-> Two Major Folder consisting MHexecution temp files and final.passed.tsv
#2.$DIR/${sample[$SLURM_ARRAY_TASK_ID]}.final.passed.tsv			-> Raw variants output file
#3.$DIR/${sample[$SLURM_ARRAY_TASK_ID]}.forAnnovar.triomode.vcf		-> List of variants with useful info only
#4.$DIR/${sample[$SLURM_ARRAY_TASK_ID]}.log				  		-> StdOutput.lof file from MH
#5.$DIR/${sample[$SLURM_ARRAY_TASK_ID]}.summary.log				-> Logfile produced by MH
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
                                        DIR=$1
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
                -r )                    shift
                                        RefGen=$1
                                        ;;
                * )                     usage
                                        exit 1
        esac
        shift
done

#define variables and directory for MosaicHunter
MHDIR=/hpcfs/groups/phoenix-hpc-neurogenetics/executables/MosaicHunter-master/MosaicHunter-master
ProbandBamFile=$(/usr/bin/find "$BAMDIR" -type f -name "$SampleID.*.bam")
MotherBamFile=$(/usr/bin/find "$BAMDIR" -type f -name "$MOTHER.*.bam")
FatherBamFile=$(/usr/bin/find "$BAMDIR" -type f -name "$FATHER.*.bam")

#module for HPCs
export HOME=/hpcfs/users/$USER
module purge
module use /apps/modules/all
module load Java/1.8.0_121
module load BLAT/3.5-foss-2016b

#1.prefilter

java -jar $MHDIR/build/mosaichunter.jar -C $MHDIR/conf/exome_parameters.properties \
-P reference_file=$RefGen \
-P input_file=$ProbandBamFile \
-P heterozygous_filter.sex=$Gender \
-P output_dir=$DIR/$SampleID.parameters.log

echo "Pre-filter completed for ${sample[$SLURM_ARRAY_TASK_ID]}"

#2.define Alpha and beta value and remove white spaces before

Al=$(cat $DIR/$SampleID.parameters.log/stdout*.log | grep "alpha" | cut -d ":" -f2)
Alpha=$(echo "$Al" | sed 's/^ *//g')

Be=$(cat $DIR/$SampleID.parameters.log/stdout*.log | grep "beta" | cut -d ":" -f2)
Beta=$(echo "$Be" | sed 's/^ *//g')

Dp=$(cat $DIR/$SampleID.parameters.log/stdout*.log | grep "average depth" | cut -d ":" -f2)
Depth=$(echo "$Dp" | sed 's/^ *//g')

#3.execute mosaic variant calling

java -jar $MHDIR/build/mosaichunter.jar -C $MHDIR/conf/exome.properties \
-P reference_file=$RefGen \
-P input_file=$ProbandBamFile \
-P depth=$Depth \
-P mosaic_filter.father_bam_file=$MotherBamFile \
-P mosaic_filter.mother_bam_file=$FatherBamFile \
-P mosaic_filter.sex=$Gender \
-P mosaic_filter.alpha_param=$Alpha \
-P mosaic_filter.beta_param=$Beta \
-P mosaic_filter.mode=trio \
-P mosaic_filter.dbsnp_file=/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/GATK4/hs37d5/b37_dbsnp_138.b37.vcf \
-P repetitive_region_filter.bed_file=$MHDIR/resources/all_repeats.b37.bed \
-P common_site_filter.bed_file=$MHDIR/resources/WES_Agilent_71M.error_prone.b37.bed \
-P output_dir=$DIR/$SampleID

echo "Somatic variant calling completed for $SampleID" 

#5.Process the outputs files

cat $DIR/$SampleID/final.passed.tsv  > $DIR/$SampleID.final.passed.tsv
cat $DIR/$SampleID/final.passed.tsv | awk '{print $1, $2, $7, $9}' | tr " " "\t" > $DIR/$SampleID.forAnnovar.triomode.vcf
#cat $MHDIR/Head.vcf > $DIR/$SampleID.final.triomode.vcf
#cat $DIR/$SampleID.forAnnovar.vcf | awk 'BEGIN{FS="\t"} {print $1, $2, ".", $3, $4, "50", "PASS", "NS=18;DP=10"}' | tr " " "\t" > $DIR/$SampleID.final.triomode.vcf

#6 log file
grep "input_file =" $DIR/$SampleID/stdout_*.log > $DIR/$SampleID.summary.log
tail  $DIR/$SampleID/stdout_*.log -n16 >> $DIR/$SampleID.summary.log
mv $DIR/$SampleID/stdout_*.log $DIR

echo "done for $SampleID"

#Remove Folders as this takes lots of memory
rm -r $DIR/$SampleID.parameters.log
rm -r $DIR/$SampleID
