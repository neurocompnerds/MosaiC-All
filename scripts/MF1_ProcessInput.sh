#!/bin/sh -l
#SBATCH -J MF1_PrepareInputFile.sh
#SBATCH -o /hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/Mosaic/MosaiC-All/TestRun/MF.Input.slurm-%j.out

#SBATCH -p skylake,icelake
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=0:05:00
#SBATCH --mem=3GB

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
#
#SCRIPT EXECUTION
#sbatch $0 -s SampleID -b BAMDIR -o OUTDIR -c CONFIG_FILE   
#all of this will be "specified" in MasterScript
#
#-s REQUIRED SampleID                       (e.g. A unique identifier for your sample that is also in the BAM file name)
#-b REQUIRED Directory of BAM files         (e.g. /path/to/bam_file_directory/)
#-c REQUIRED The script configuration file  (/path/to/config_file.cfg)
#-o REQUIRED Directory of outputs           (e.g. /path/to/output/)
#
"
}

## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -s )                    shift
                                        SampleID=$1
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

## Check for all required arguments
for ARG in $BAMDIR $OUTDIR $CONFIG_FILE $SampleID; do
    if [ -z "${ARG}" ]; then
        usage
        echo "## ERROR: Please check that you have supplied all required argunments for this script"
        exit 1
    fi
done

source $CONFIG_FILE

#process Mutect2 files
BAMFILE=$(find "$BAMDIR" -type f -name "$SampleID*.bam")
BAMprefix=$(basename "$BAMFILE" | sed 's/\.[^.]*$//')

tmp_file=$(mktemp)  # Create a temporary file

grep -v "#" $OUTDIR/$SampleID.$CONFIG.PON.gnomad.vcf | tr " " "\t" > "$tmp_file"

awk -F  '\t' '{print $1, $2-1, $2, $4, $5}' "$tmp_file" > "$OUTDIR/$SampleID.forPhasing.bed"

rm "$tmp_file"  # Remove the temporary file

awk -v prefix="$BAMprefix" 'BEGIN{OFS="\t"} {$6 = prefix; print}' $OUTDIR/$SampleID.forPhasing.bed > $OUTDIR/$SampleID.phasingInput.bed
