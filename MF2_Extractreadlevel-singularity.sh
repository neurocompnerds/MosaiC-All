#!/bin/sh -l
#SBATCH -J MF_Extractreadlevel-singularity.sh
#SBATCH -o /hpcfs/groups/phoenix-hpc-neurogenetics/Nandini/Mosaic-All/Log/MF.extract.slurm-%j.out
#SBATCH -A robinson
#SBATCH -p skylake,icelake
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=07:00:00
#SBATCH --mem=40GB
#SBATCH --gres=tmpfs:150G
# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=a1742674@adelaide.edu.au

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
#sbatch $SCRIPTDIR/MF2_Extractreadlevel-singularity.sh -b BAMDIR -s SAMPLE -c CONFIG_FILE -o OUTDIR  
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
                -s )                    shift
                                        SAMPLE=$1
                                        ;;
                -b )                    shift
                                        BAMDIR=$1
                                        ;;
                -o )                   	shift
                                        OUTDIR=$1
										;;
				-c )                   	shift
                                        CONFIG_FILE=$1
										;;
                 * )                    usage
                                        exit 1
        esac
        shift
done


#define directories
export HOME=/gpfs/users/$USER

module purge
module use /apps/skl/modules/all
module load Singularity/3.7.1

source $CONFIG_FILE

BAMFILE=$(find "$BAMDIR -type f -name "$SAMPLE*.bam")
BAMprefix=$(basename "$BAMFILE" | sed 's/\.[^.]*$//')

/usr/bin/mkdir -p ${TMPDIR}
/usr/bin/cp -r ${BAMFILE} ${TMPDIR}
/usr/bin/cp -r ${BAMINDEX} ${TMPDIR}
/usr/bin/cp -r $OUTDIR/$SAMPLE.phasingInput.bed ${TMPDIR}


#execute the read-level extractions
singularity run -B /hpcfs -B ${TMPDIR}:/tmp $MFORECAST/mosaicforecast_0.0.1.sif ReadLevel_Features_extraction.py /tmp/$SAMPLE.phasingInput.bed /tmp/$SAMPLE.features.bed /tmp $REFGEN $MFORECAST/hg19/k24.umap.wg.bw ${SLURM_NTASKS} bam


/usr/bin/cp -r ${TMPDIR}/$SAMPLE.features.bed ${OUTDIR}/$SAMPLE.features.bed

echo "SUCCESS for $SAMPLE"
