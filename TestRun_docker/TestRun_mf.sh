#####################################################################   MOSAICFORECAST   ################################################################
#!/bin/bash

# Define variables
GIT="/home/neuro/Documents/Nandini/MosaiC-All"
TESTRUN_DIR="$GIT/TestRun_docker"
BAMDIR="$TESTRUN_DIR/BAM"
ProbandBamFile="$BAMDIR/1465_1024-pfc-bulk.rehead.bam"
ProbandID="1465_1024-pfc-bulk"
Gender="M"
OUTDIR="$TESTRUN_DIR/OUTPUT"
RESOURCES="$TESTRUN_DIR/Resources"
PON="${RESOURCES}/somatic-b37_Mutect2-exome-panel.vcf"
GERMLINE_RESOURCES="${RESOURCES}/somatic-b37_af-only-gnomad.raw.sites.vcf"
CONFIG="hg19"
REFGEN="/home/neuro/Public/RefSeqIndexAllPrograms"


# Create output directory
mkdir -p "$OUTDIR"

# Pull Docker image
docker image pull yanmei/mosaicforecast:0.0.1 || { echo "Failed to pull MosaicForecast Docker image."; exit 1; }

# Ensure BAM file exists and is unique
BAMFILE=$(find "$BAMDIR" -type f -name "$ProbandID*.bam")
if [ -z "$BAMFILE" ]; then
    echo "No BAM file found for $ProbandID in $BAMDIR."
    exit 1
elif [ $(echo "$BAMFILE" | wc -l) -gt 1 ]; then
    echo "Multiple BAM files found for $ProbandID in $BAMDIR. Please ensure only one exists."
    exit 1
fi

# Assign BAM prefix
BAMprefix=$(basename "$BAMFILE" | sed 's/\.[^.]*$//')

# Create temporary file
tmp_file=$(mktemp) || { echo "Failed to create temporary file."; exit 1; }

# Process Mutect2 files
grep -v "#" "$OUTDIR/$ProbandID.Mutect2.$CONFIG.PON.gnomad.vcf" | tr " " "\t" > "$tmp_file"
awk -F '\t' '{print $1, $2-1, $2, $4, $5}' "$tmp_file" > "$OUTDIR/$ProbandID.forPhasing.bed"
rm "$tmp_file"  # Clean up temporary file

# Ensure Docker image was pulled successfully
docker image pull yanmei/mosaicforecast:0.0.1 || { echo "Failed to pull MosaicForecast Docker image."; exit 1; }

# Run Docker command for ReadLevel Features extraction
docker run -it --rm \
  -v "$GIT:/mnt/mosaic-all" \
  -v "$TESTRUN_DIR:/mnt/testrun_docker" \
  -v "$REFGEN:/mnt/refgen" \
  yanmei/mosaicforecast:0.0.1 \
  /bin/bash -c " \
  MForecastDIR=/usr/MosaicForecast \
  ReadLevel_Features_extraction.py \
  /mnt/testrun_docker/OUTPUT/$ProbandID.MF.$CONFIG.phasingInput.bed \
  /mnt/testrun_docker/OUTPUT/$ProbandID.MF.$CONFIG.features.bed \
  /mnt/refgen/hs37d5.fa \
  /MForecastDIR/$CONFIG/k24.umap.wg.bw \
  2 \
  bam \
" || { echo "ExtractFeatures for MosaicForecast failed."; exit 1; }


# MF3.Genotype prediction
docker run -it --rm \
  -v "$GIT:/mnt/mosaic-all" \
  -v "$TESTRUN_DIR:/mnt/testrun_docker" \
  -v "$REFGEN:/mnt/refgen" \
  yanmei/mosaicforecast:0.0.1 \
  /bin/bash -c " \
  MForecastDIR=/usr/MosaicForecast \
  Prediction.R \
  /mnt/testrun_docker/OUTPUT/$ProbandID.MF.$CONFIG.features.bed \
  /MForecastDIR/models_trained/50xRFmodel_addRMSK_Refine.rds \
  Refine \
  /mnt/testrun_docker/OUTPUT/$ProbandID.MF.$CONFIG.genotype.predictions.refined.bed \
" || { echo "GenotypePredictions for MosaicForecast failed."; exit 1; }
