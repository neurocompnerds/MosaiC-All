#!/bin/sh
#SBATCH -J FilterMutect.sh
#SBATCH -o /hpcfs/groups/phoenix-hpc-neurogenetics/Nandini/Mosaic-All/Log/filter-slurm-%j.out
#SBATCH -A robinson
#SBATCH -p skylake,icelake
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=00:30:00
#SBATCH --mem=10GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nandini.sandran@adelaide.edu.au


usage()
{
echo "
#
#README
#This script is designed for Mosaic-All wrapper script
#For standalone script,look for the original
#To filter Mutect2 calls according to criteria set
#
#SCRIPT EXECUTION
#SCRIPTDIR=/hpcfs/groups/phoenix-hpc-neurogenetics/Mosaic-All/Mosaic-S
#sbatch $SCRIPTDIR/Mutect2.FilterMutect2.sh -s $SAMPLE -v $VCFDIR -c $CONFIG 
#all of this will be "specified" in MasterScript
#
#-s REQUIRED SampleID   					(e.g 004P, or 004M, or 004F no suffix)
#-v REQUIRED Directory of VCFs 				(directory of outputs of Mutect2) 		
#-c REQUIRED Config_File    	        	(e.g hs37d5.fa) 
"
}


## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -s )                    shift
                                        SAMPLE=$1
                                        ;;
                -v )                    shift
                                        VCFDIR=$1
                                        ;;
                -c )                    shift
                                        CONFIG_FILE=$1
                                        ;;
                * )                     usage
                                        exit 1
        esac
        shift
done

export HOME=/hpcfs/users/$USER
module purge 
module load GATK/4.4.0.0-GCCcore-11.2.0-Java-17.0.6
#module load SAMtools/1.17-GCC-11.2.0

#check and specify
source $CONFIG_FILE

if [ -z $VCFDIR/$SAMPLE.$CONFIG.PON.gnomad.vcf ]; then # If no vcfDir name specified then do not proceed
        usage
        echo "#ERROR: You need to tell me where to find the vcf files."
        exit 1
fi

## Start of the script ##
###On each sample###

#1. To "tag" each variants

gatk FilterMutectCalls \
-V $VCFDIR/$SAMPLE.$CONFIG.PON.gnomad.vcf \
-R $REFGEN \
--min-reads-per-strand 2 \
-O $VCFDIR/$SAMPLE.singlemode.filtered.vcf > $VCFDIR/$SAMPLE.filtered.log 2>&1


#2. To filter manually

module load BCFtools/1.17-GCC-11.2.0

bcftools view -O v -f PASS -i '(FORMAT/AD[0:1] >= 5 & FORMAT/DP[0]>=20) && (FORMAT/AF[0:0] <=0.4 || FORMAT/AF[0:0]>=0.7)' $VCFDIR/$SAMPLE.singlemode.filtered.vcf > $VCFDIR/$SAMPLE.mutect2.singlemode.PASS.aaf.vcf



