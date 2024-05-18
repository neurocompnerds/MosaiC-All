
#####################################################################   MUTECT2   #####################################################################

# Define variables
GIT="/home/neuro/Documents/Nandini/MosaiC-All"
TESTRUN_DIR="$GIT/TestRun_docker"
BAMDIR="$TESTRUN_DIR/BAM"
ProbandBamFile="$BAMDIR/1465_1024-pfc-bulk.rehead.bam"
ProbandID="1465_1024-pfc-bulk"
Gender="M"
OUTDIR="OUTPUT"
RESOURCES="$TESTRUN_DIR/Resources"
PON="Mutect2-exome-panel.vcf"
GERMLINE_RESOURCES="af-only-gnomad.raw.sites.vcf"
CONFIG="hg19"
REFGEN=/home/neuro/Public/RefSeqIndexAllPrograms

# Create output directory
mkdir -p "$TESTRUN_DIR/$OUTDIR"

#Download resources

#wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf -P $RESOURCES
#wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf.idx -P $RESOURCES
#wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf -P $RESOURCES
#wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx -P $RESOURCES

# Pull Docker image
docker image pull broadinstitute/gatk:4.4.0.0 || { echo "Failed to pull Docker image."; exit 1; }

# Mutect2.1 Variant calling
docker run -it --rm \
  -v "$GIT:/mnt/mosaic-all" \
  -v "$TESTRUN_DIR:/mnt/testrun_docker" \
  -v "$REFGEN:/mnt/refgen" \
  broadinstitute/gatk:4.4.0.0 \
  /bin/bash -c "
  gatk Mutect2 \
    -R /mnt/refgen/hs37d5.fa \
    -I /mnt/testrun_docker/BAM/1465_1024-pfc-bulk.rehead.bam  \
    --tumor-sample $ProbandID \
    --germline-resource /mnt/testrun_docker/Resources/$GERMLINE_RESOURCES \
    --panel-of-normals /mnt/testrun_docker/Resources/$PON \
    --af-of-alleles-not-in-resource -1 \
    -O /mnt/testrun_docker/$OUTDIR/$ProbandID.Mutect2.$CONFIG.PON.gnomad.vcf
" || { echo "Mutect2.1 Variant calling failed."; exit 1; }

# Mutect2.2 Tag each variant

docker run -it --rm \
  -v "$GIT:/mnt/mosaic-all" \
  -v "$TESTRUN_DIR:/mnt/testrun_docker" \
  -v "$REFGEN:/mnt/refgen" \
  broadinstitute/gatk:4.4.0.0 \
  /bin/bash -c "
  gatk FilterMutectCalls \
    -V /mnt/testrun_docker/$OUTDIR/$ProbandID.Mutect2.$CONFIG.PON.gnomad.vcf \
    -R /mnt/refgen/hs37d5.fa \
    --min-reads-per-strand 2 \
    -O /mnt/testrun_docker/$OUTDIR/$ProbandID.Mutect2.$CONFIG.singlemode.filtered.vcf > /mnt/testrun_docker/$OUTDIR/$ProbandID.Mutect2.$CONFIG.filtered.log 2>&1
" || { echo "Mutect2.2 Tag each variant failed."; exit 1; }



# Mutect2.3 Manual Filtration
bcftools view -O v -f PASS -i '(FORMAT/AD[0:1] >= 5 & FORMAT/DP[0]>=20) && (FORMAT/AF[0:0] <=0.4 || FORMAT/AF[0:0]>=0.7)' "$OUTDIR/$ProbandID.Mutect2.$CONFIG.singlemode.filtered.vcf" > "$OUTDIR/$ProbandID.Mutect2.$CONFIG.singlemode.PASS.aaf.vcf" || { echo "Mutect2.3 Manual Filtration failed."; exit 1; }



#####################################################################   MOSAICFORECAST   ################################################################

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

# Ensure BAM file exists and is unique
BAMFILE=$(find "$BAMDIR" -type f -name "$ProbandID*.bam")
if [ -z "$BAMFILE" ]; then
    echo "No BAM file found for $ProbandID in $BAMDIR."
    exit 1
elif [ $(echo "$BAMFILE" | wc -l) -gt 1 ]; then
    echo "Multiple BAM files found for $ProbandID in $BAMDIR. Please ensure only one exists."
    exit 1
fi

#Resources

# requires model_trained folder from MF which cannot be found in docker
git clone --depth 1 --branch master https://github.com/parklab/MosaicForecast.git --single-branch $RESOURCES

# Download wigtoBigWig
wget -q http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig -P "$RESOURCES"
chmod +x "$RESOURCES/wigToBigWig"

# Download fetchChromSizes and make it executable
wget -q http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes -P "$RESOURCES"
chmod +x "$RESOURCES/fetchChromSizes"

# Check if fetchChromSizes was downloaded and made executable successfully
if [ -x "$RESOURCES/fetchChromSizes" ]; then
    echo "fetchChromSizes downloaded and made executable successfully."
else
    echo "Error: Failed to download or make fetchChromSizes executable."
    exit 1
fi

# Download hg19.umap.tar.gz
wget -q https://bismap.hoffmanlab.org/raw/hg19.umap.tar.gz -P "$RESOURCES"

# Check if hg19.umap.tar.gz was downloaded successfully
if [ -e "$RESOURCES/hg19.umap.tar.gz" ]; then
    # Extract the tar.gz file
    tar -zxvf "$RESOURCES/hg19.umap.tar.gz" -C "$RESOURCES"
    
    # Check if extraction was successful
    if [ $? -eq 0 ]; then
        echo "Extraction of hg19.umap.tar.gz completed successfully."
    else
        echo "Error: Failed to extract hg19.umap.tar.gz."
        exit 1
    fi

    # Navigate to the directory containing extracted files
    cd "$RESOURCES/hg19" || exit

    # Generate chrom sizes file
    "$RESOURCES/fetchChromSizes" hg19 > "$RESOURCES/hg19"/hg19.chrom.sizes

    # Check if chrom sizes file was generated successfully
    if [ -e "$RESOURCES/hg19/hg19.chrom.sizes" ]; then
        echo "Chrom sizes file generated successfully."
        # Convert wig to bigwig
        zcat "$RESOURCES/hg19/k24.umap.wg.gz" > "$RESOURCES/hg19/k24.umap.wg"
        "$RESOURCES/wigToBigWig" "$RESOURCES/hg19/k24.umap.wg" "$RESOURCES/hg19/hg19.chrom.sizes" "$RESOURCES/hg19/k24.umap.wg.bw"
    else
        echo "Error: Failed to generate chrom sizes file."
        exit 1
    fi
else
    echo "Error: Failed to download hg19.umap.tar.gz."
    exit 1
fi


# MF1.Prepare Input files
# Assign BAM prefix
BAMprefix=$(basename "$BAMFILE" | sed 's/\.[^.]*$//')

# Create temporary file
tmp_file=$(mktemp) || { echo "Failed to create temporary file."; exit 1; }

# Process Mutect2 files
grep -v "#" "$OUTDIR/$ProbandID.Mutect2.$CONFIG.PON.gnomad.vcf" | tr " " "\t" > "$tmp_file"
awk -F '\t' '{print $1, $2-1, $2, $4, $5}' "$tmp_file" > "$OUTDIR/$ProbandID.forPhasing.bed"
awk -v prefix="$BAMprefix" 'BEGIN{OFS="\t"} {$6 = prefix; print}' $OUTDIR/$ProbandID.forPhasing.bed > $OUTDIR/$ProbandID.MF.$CONFIG.phasingInput.bed
rm "$tmp_file"  # Clean up temporary file

# Ensure Docker image was pulled successfully
docker image pull yanmei/mosaicforecast:0.0.1 || { echo "Failed to pull MosaicForecast Docker image."; exit 1; }

# MF2.Run Docker command for ReadLevel Features extraction
docker run -it --rm \
  -v "$GIT:/MF" \
  -v "$REFGEN:/mnt/refgen" \
  -v "$TESTRUN_DIR:/mnt/testrun_docker" \
  yanmei/mosaicforecast:0.0.1 \
  /bin/bash -c " \
  ReadLevel_Features_extraction.py \
  /mnt/testrun_docker/OUTPUT/$ProbandID.MF.$CONFIG.phasingInput.bed \
  /mnt/testrun_docker/OUTPUT/$ProbandID.MF.$CONFIG.features.bed \
  /mnt/testrun_docker/BAM \
  /mnt/refgen/hs37d5.fa \
  /mnt/testrun_docker/Resources/$CONFIG/k24.umap.wg.bw \
  2 \
  bam \
" || { echo "ExtractFeatures for MosaicForecast failed."; exit 1; }

# MF3.Genotype prediction
docker run -it --rm \
  -v "$GIT:/MF" \
  -v "$TESTRUN_DIR:/mnt/testrun_docker" \
  -v "$REFGEN:/mnt/refgen" \
  yanmei/mosaicforecast:0.0.1 \ 
  /bin/bash -c " \
  Prediction.R \
  /mnt/testrun_docker/OUTPUT/$ProbandID.MF.$CONFIG.features.bed \
  /mnt/testrun_docker/Resources/MF/models_trained/50xRFmodel_addRMSK_Refine.rds \
  Refine \
  /mnt/testrun_docker/OUTPUT/$ProbandID.MF.$CONFIG.genotype.predictions.refined.bed \
" || { echo "GenotypePredictions for MosaicForecast failed."; exit 1; }

