#!/bin/sh -l
#SBATCH -J MF1_PrepareInputFile.sh
#SBATCH -o /hpcfs/users/%u/Mosaic-All/Log/MF.Input.slurm-%j.out
#SBATCH -A robinson
#SBATCH -p skylake,icelake
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=0:05:00
#SBATCH --mem=3GB

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
#sbatch $SCRIPTDIR/MF1_ProcessInput.sh -s SAMPLE -b BAMDIR -o OUTDIR -c CONFIG_FILE   
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
                -b )                   	shift
                                        BAMDIR=$1
                                        ;;
                -o )                   	shift
                                        OUTDIR=$1
                                        ;;
                -c )                   	shift
                                        CONFIG_FILE=$1
                                        ;;
                * )                     usage
                                        exit 1
        esac
        shift
done


module purge

source $CONFIG_FILE

#process Mutect2 files
BAMFILE=$(find "$BAMDIR" -type f -name "$SAMPLE*.bam")
BAMprefix=$(basename "$BAMFILE" | sed 's/\.[^.]*$//')

tmp_file=$(mktemp)  # Create a temporary file

grep -v "#" $OUTDIR/$SAMPLE.$CONFIG.PON.gnomad.vcf | tr " " "\t" > "$tmp_file"

awk -F  '\t' '{print $1, $2-1, $2, $4, $5}' "$tmp_file" > "$OUTDIR/$SAMPLE.forPhasing.bed"

rm "$tmp_file"  # Remove the temporary file

awk -v prefix="$BAMprefix" 'BEGIN{OFS="\t"} {$6 = prefix; print}' $OUTDIR/$SAMPLE.forPhasing.bed > $OUTDIR/$SAMPLE.phasingInput.bed

