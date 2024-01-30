#!/bin/sh
#SBATCH -J FilterMutect.sh
#SBATCH -o /home/%u/Mosaic-All/Log/filter-slurm-%j.out

#SBATCH -p skylake,icelake
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=00:30:00
#SBATCH --mem=10GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=%u@adelaide.edu.au


usage()
{
echo "
#
#README
#This script is designed for Mosaic-All wrapper script
#To filter Mutect2 calls according to criteria set
#
#SCRIPT EXECUTION
#sbatch $0 -s \$SampleIDID -v /path/to/Mutect2_VCF_dir/ -c /path/to/config_file.cfg
#
#-s REQUIRED SampleID                       (e.g. A unique identifier for your sample)
#-v REQUIRED Directory of VCFs              (directory of outputs of Mutect2) 		
#-c REQUIRED The script configuration file  (/path/to/config_file.cfg)
# 
"
}


## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -s )                    shift
                                        SampleID=$1
                                        ;;
                -v )                    shift
                                        VCFDIR=$1  # Note: This location is the same as $OUTDIR in MASTERSCRIPT_Mosaic_All.sh
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
for ARG in $VCFDIR $CONFIG_FILE $SampleID; do
    if [ -z "${ARG}" ]; then
        usage
        echo "## ERROR: Please check that you have supplied all required argunments for this script" >> $VCFDIR/$SampleID.pipeline.log
        exit 1
    fi
done

#check and specify
source $CONFIG_FILE

if [ -z "$VCFDIR/$SampleID.$CONFIG.PON.gnomad.vcf" ]; then # If no vcfDir name specified then do not proceed
        usage
        echo "## ERROR: The file $VCFDIR/$SampleID.$CONFIG.PON.gnomad.vcf does not exist." >> $VCFDIR/$SampleID.pipeline.log
        exit 1
fi

## Load modules
module purge 
module load GATK/4.4.0.0-GCCcore-11.2.0-Java-17.0.6

## Start of the script ##
###On each sample###

#1. To "tag" each variants

gatk FilterMutectCalls \
-V $VCFDIR/$SampleID.$CONFIG.PON.gnomad.vcf \
-R $REFGEN \
--min-reads-per-strand 2 \
-O $VCFDIR/$SampleID.singlemode.filtered.vcf > $VCFDIR/$SampleID.filtered.log 2>&1


#2. To filter manually

module load BCFtools/1.17-GCC-11.2.0

bcftools view -O v -f PASS -i '(FORMAT/AD[0:1] >= 5 & FORMAT/DP[0]>=20) && (FORMAT/AF[0:0] <=0.4 || FORMAT/AF[0:0]>=0.7)' $VCFDIR/$SampleID.singlemode.filtered.vcf > $VCFDIR/$SampleID.mutect2.singlemode.PASS.aaf.vcf
