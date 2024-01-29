#!/bin/sh
#SBATCH -J CombineCalls
#SBATCH -o /hpcfs/users/%u/MosaiC-All/Log/CC-slurm-%j.out
#SBATCH -A robinson
#SBATCH -p skylake,icelake
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=00:05:00
#SBATCH --mem=1GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=%u@adelaide.edu.au


usage()
{
echo "
#
#README
#This script is designed to be used after Mosaic-All wrapper script,
#Post-processing MF variants
#to COMBINE all variant calls (of all samples executed by batch) in one meta-file based on each mosaic variant calling
#for M3 pipeline
#sbatch $SCRIPTDIR/CombineCalls.sh
#-s REQUIRED sampleID 
#-d REQUIRED DIRECTORY OF VARIANT LIST OF EACH SAMPLE BEFORE MERGING
"
}


## Set Variables ##
while [ "$1" != "" ]; do
    case $1 in
        -s ) shift
             SAMPLEID=$1
             ;;
	-d ) shift
             DIR=$1
             ;;
        * ) usage
            exit 1
    esac
    shift
done

# Check if input files exist
for FILE in "$DIR/$SAMPLEID.$MUT" "$DIR/$SAMPLEID.$MH" "$DIR/$SAMPLEID.$MF"; do
    if [ ! -f "$FILE" ]; then
        echo "Error: File $FILE does not exist."
        exit 1
    fi
done

##Define Suffix of each output files
MUT="mutect2.singlemode.PASS.aaf.vcf"
MH="final.passed.tsv" 
MF="Refined.dpaf.MosaicOnly.txt"

#Mutect2
awk -v ID="$SAMPLEID" '$0 !~ /^##/ {print ID "\t" "Mutect2" "\t" $0}' $DIR/$SAMPLEID.$MUT >> $DIR/Mutect2.variants.txt

#MosaicHunter
awk -v ID="$SAMPLEID" '$0 !~ /^##/ {print ID "\t" "MH" "\t" $0}' $DIR/$SAMPLEID.$MH >> $DIR/MH.variants.txt

#MosaicForecast
awk '$35=="mosaic" {OFS="\t"; print}' $DIR/$SAMPLEID.genotype.predictions.phased.singlemode.bed > $DIR/$SAMPLEID.RefinedMosaicOnly.txt
awk '$25>=20  {OFS="\t"; print}' $DIR/$SAMPLEID.RefinedMosaicOnly.txt > $DIR/$SAMPLEID.Refined.dp20.MosaicOnly.txt
awk '$24>=0.03  {OFS="\t"; print}' $DIR/$SAMPLEID.Refined.dp20.MosaicOnly.txt > $DIR/$SAMPLEID.Refined.dpaf.MosaicOnly.txt
awk -v ID="$SAMPLEID" '$0 !~ /^##/ {print ID "\t" "MF" "\t" $0}' $DIR/$SAMPLEID.$MF >> $DIR/MF.variants.txt
