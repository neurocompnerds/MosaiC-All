#!/bin/sh
#SBATCH -J CombineCalls
#SBATCH -o /hpcfs/groups//phoenix-hpc-neurogenetics/a1742674/hg38_OUTPUTS/CC-slurm-%j.out

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
#sbatch $0 -s SampleID -d /path/to/M3_pipeline/output
#-s REQUIRED sampleID 
#-d REQUIRED Directory containing the lists of variants called for each sample before merging
"
}


## Set Variables ##
while [ "$1" != "" ]; do
    case $1 in
        -s ) shift
             SampleID=$1
             ;;
        -d ) shift
             DIR=$1
             ;;
        * ) usage
            exit 1
    esac
    shift
done



## Check for all required arguments
for ARG in $SampleID $DIR; do
    if [ -z "${ARG}" ]; then
        usage
        echo "## ERROR: Please check that you have supplied all required argunments for this script" > $DIR/$SampleID.postprocessing.log
        exit 1
    fi
done

##Define Suffix of each output files
MUT="mutect2.singlemode.PASS.aaf.vcf"
MH="final.passed.tsv" 
MF="Refined.dpaf.MosaicOnly.txt"

# Check if input files exist
for FILE in "$DIR/$SampleID.$MUT" "$DIR/$SampleID.$MH" "$DIR/$SampleID.mosaicforecast.genotype.predictions.refined.bed"; do
    if [ ! -f "$FILE" ]; then
        echo "Error: File $FILE does not exist." >> $DIR/$SampleID.postprocessing.log
        exit 1
    fi
done

#Mutect2
awk -v ID="$SampleID" '$0 !~ /^##/ {print ID "\t" "Mutect2" "\t" $0}' $DIR/$SampleID.$MUT >> $DIR/Mutect2.variants.txt

#MosaicHunter
awk -v ID="$SampleID" '$0 !~ /^##/ {print ID "\t" "MH" "\t" $0}' $DIR/$SampleID.$MH >> $DIR/MH.variants.txt

#MosaicForecast
MF_output=$SampleID.mosaicforecast.genotype.predictions.refined.bed
awk '$35=="mosaic" {OFS="\t"; print}' $DIR/$MF_output > $DIR/$SampleID.RefinedMosaicOnly.txt
awk '$25>=20  {OFS="\t"; print}' $DIR/$SampleID.RefinedMosaicOnly.txt > $DIR/$SampleID.Refined.dp20.MosaicOnly.txt
awk '$24>=0.03  {OFS="\t"; print}' $DIR/$SampleID.Refined.dp20.MosaicOnly.txt > $DIR/$SampleID.Refined.dpaf.MosaicOnly.txt
awk -v ID="$SampleID" '$0 !~ /^##/ {print ID "\t" "MF" "\t" $0}' $DIR/$SampleID.$MF >> $DIR/MF.variants.txt
