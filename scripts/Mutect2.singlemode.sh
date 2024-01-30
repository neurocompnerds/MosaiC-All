#!/bin/sh
#SBATCH -J Mutect2.singlemode.sh
#SBATCH -o /home/%u/Mosaic-All/Log/Mutect2-slurm-%j.out

#SBATCH -p skylake,icelake
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=8:00:00
#SBATCH --mem=100GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=%u@adelaide.edu.au

#DATE:28 JUNE 2023

usage()
{
echo "
#
#README
#This script is designed for Mosaic-All wrapper script
#
#SCRIPT EXECUTION
#sbatch $0 -b /path/to/bam_file_directory/ -s SampleID -c /path/to/config_file.cfg -o /path/to/output/  -p /path/to/panel_of_normals.vcf
#
#-s REQUIRED SampleID                       (e.g. A unique identifier for your sample that is also in the BAM file name)
#-b REQUIRED Directory of BAM files         (e.g. /path/to/bam_file_directory/)
#-c REQUIRED The script configuration file  (/path/to/config_file.cfg)
#-o REQUIRED Directory of outputs           (e.g. /path/to/output/) 
#-p REQUIRED Panel of Normals vcf file      (/path/to/panel_of_normals.vcf)
"
}

## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -b )                    shift
                                        BAMDIR=$1
                                        ;;
		-s )			shift
			               	SampleID=$1
                                        ;;
                -c )                    shift
                                        CONFIG_FILE=$1
                                        ;;
                -o )                    shift
                                        OUTDIR=$1
					;;
                -p )                    shift
                                        PON=$1
					;;
		 * )			usage
			                exit 1
        esac
        shift
done

## Check for all required arguments
for ARG in $BAMDIR $OUTDIR $CONFIG_FILE $SampleID $PON; do
    if [ -z "${ARG}" ]; then
        usage
        echo "## ERROR: Please check that you have supplied all required argunments for this script"
        exit 1
    fi
done

module purge
module use /apps/modules/all
module load SAMtools/1.8-foss-2016b
module load GATK/4.4.0.0-GCCcore-11.2.0-Java-17.0.6

source $CONFIG_FILE
ProbandBamFile=$(find "$BAMDIR" -type f -name "$SampleID.*.bam")

#execute the script
gatk Mutect2 \
-R $REFGEN \
-I $ProbandBamFile  \
--tumor-sample $SampleID \
--germline-resource $GERMLINE_RESOURCES \
--panel-of-normals $PON \
--af-of-alleles-not-in-resource -1 \
-O $OUTDIR/$SampleID.$CONFIG.PON.gnomad.vcf



