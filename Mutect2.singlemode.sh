#!/bin/sh
#SBATCH -J Mutect2.singlemode.sh
#SBATCH -o /hpcfs/groups/phoenix-hpc-neurogenetics/Nandini/Mosaic-All/Log/Mutect2-slurm-%j.out
#SBATCH -A robinson
#SBATCH -p skylake,icelake
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=8:00:00
#SBATCH --mem=100GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nandini.sandran@adelaide.edu.au

#DATE:28 JUNE 2023
#modified for wrapper script designed for HPC
#original:/hpcfs/users/a1742674/Scripts
#DATE:30 APRIL 2022
#GATK 4.2.6.1, PON(LATEST) and germlineresc
#DATE:13 January 2022
#modified script for new twins cohort
##removed somatic_b37 (commonsnps)
##genome: hs38DH.fa
#DATE:19 February 2021
## modified the script by using the inhouse PON (created from the unaffected parental samples -352)
#DATE:12 February 2021
##modified script for sbatch array
##Usage sbatch --array 0-(nSamples-1) $0  
#DATE:11 February 2021
## created test-script to run Mutect2 in single-sample mode for a single sample
## Resources 
## (i) germline variants downlaoded from gnomad website
## (ii) panel or normals - gatk (the website mentioned that gatk will be able to detect the file directly)
##choice of ref:/hpcfs/users/a1742674/2021/RefSeq/human_g1k_v37.fasta or /hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/hs37d5.fa


usage()
{
echo "
#
#README
#This script is designed for Mosaic-All wrapper script
#For standalone script,look for the original
#
#SCRIPT EXECUTION
#SCRIPTDIR=/hpcfs/groups/phoenix-hpc-neurogenetics/Mosaic-All/Mosaic-S
#sbatch $SCRIPTDIR/Mutect2.singlemode.sh -b BAMDIR -s SAMPLE -c CONFIG_FILE -o OUTDIR  
#all of this will be "specified" in MasterScript
#
#-s REQUIRED SampleID   					(e.g 004P, or 004M, or 004F no suffix)
#-b REQUIRED Directory of Bam files 		(e.g /hpcfs/groups/phoenix-hpc-sacgf/scratch/ali/GA_bams)
#-c REQUIRED Directory of Outputs   		(/gpfs/users/a1149120/MosaicHunter_single)
#-o REQUIRED Reference Genome     	        (e.g hs37d5.fa) 
"
}


## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -b )                    shift
                                        BAMDIR=$1
                                        ;;
				-s )					shift
			               			    SAMPLE=$1
                                        ;;
                -c )                    shift
                                        CONFIG_FILE=$1
                                        ;;
                -o )                    shift
                                        OUTDIR=$1
										;;
				 * )					usage
			                			exit 1
        esac
        shift
done

export HOME=/hpcfs/users/$USER
module purge

module use /apps/modules/all
#module load Java/1.8.0_121
module load SAMtools/1.8-foss-2016b
module load GATK/4.4.0.0-GCCcore-11.2.0-Java-17.0.6

source $CONFIG_FILE

#execute the script
gatk Mutect2 \
-R $REFGEN \
-I $BAMDIR/$SAMPLE.realigned.recal.sorted.bwa.hs37d5.bam  \
--tumor-sample $SAMPLE \
--germline-resource $GERMLINE_RESOURCES \
--panel-of-normals $PON \
--af-of-alleles-not-in-resource -1 \
-O $OUTDIR/$SAMPLE.$CONFIG.PON.gnomad.vcf



