#####################################################################   MOSAICFORECAST   ################################################################
#!/bin/bash

# Define variables
#git_dir="/home/neuro/Documents/Nandini/MosaiC-All"
git_dir=/Users/nandinisandran/Desktop/MosaiC-All/TestRun_docker

# Other variables should not be changed
testrun_dir="$git_dir/TestRun_docker"
bam_dir="$testrun_dir/BAM"
proband_bam_file="$bam_dir/1465_1024-pfc-bulk.rehead.bam"
proband_id="1465_1024-pfc-bulk"
gender="M"
config="hg19"
outdir="$testrun_dir/OUTPUT" 
resources="$testrun_dir/Resources"

if [ ! -d "${outdir}" ]; then
    mkdir -p ${outdir}
fi

if [ ! -d "${resources}" ]; then
    mkdir -p ${resources}
fi

# download, extract and create index for reference genome

if [ ! -e "$resources/hs37d5.fa.gz" ]; then
    # download
    wget "https://storage.googleapis.com/genomics-public-data/references/hs37d5/hs37d5.fa.gz" -P "$resources"
    # Extract the tar.gz file
    gunzip "$resources/hs37d5.fa.gz" -c "resources"
    samtools faidx "hs37d5.fa.gz"
fi

# Ensure Docker image was pulled successfully
docker image pull yanmei/mosaicforecast:0.0.1 || { echo "Failed to pull MosaicForecast Docker image."; exit 1; }

# download resources for MosaicForecast run

# 1/2 "model_trained" folder from MF which cannot be found in docker
git clone --depth 1 --branch master https://github.com/parklab/MosaicForecast.git --single-branch $resources

# 2/2 Umap score (k=24, GRCh37/hg19)
wget -q https://bismap.hoffmanlab.org/raw/hg19.umap.tar.gz -P "$resources"

# Check if hg19.umap.tar.gz was downloaded successfully
if [ -e "$resources/hg19.umap.tar.gz" ]; then
    # Extract the tar.gz file
    tar -zxvf "$resources/hg19.umap.tar.gz" -C "$resources"
    
    # Check if extraction was successful
    if [ $? -eq 0 ]; then
        echo "Extraction of hg19.umap.tar.gz completed successfully."
    else
        echo "Error: Failed to extract hg19.umap.tar.gz."
        exit 1
    fi
fi

# Run the Docker container
docker run -it --rm \
  -v "$git:/MF" \
  -v "$resources:/mnt/resources" \
  yanmei/mosaicforecast:0.0.1 \
  /bin/bash -c ' 
    # Fetch chromosome sizes
    fetchChromSizes /mnt/resources/hg19/hg19 > /mnt/resources/hg19/hg19.chrom.sizes 
    && 
    # Decompress the Umap score file
    zcat /mnt/resources/hg19/k24.umap.wg.gz > /mnt/resources/hg19/k24.umap.wg 
    && 
    # Convert wig to BigWig format
    wigToBigWig /mnt/resources/hg19/k24.umap.wg /mnt/resources/hg19/hg19.chrom.sizes /mnt/resources/hg19/k24.umap.wg.bw 
' || { echo "Processing Umap score for MosaicForecast failed."; exit 1; }


# MF1.Prepare Input files
# Assign BAM prefix
BAMprefix=$(basename "$proband_bam_file" | sed 's/\.[^.]*$//')
echo "$BAMprefix"

# Create temporary file
tmp_file=$(mktemp) || { echo "Failed to create temporary file."; exit 1; }

# Process Mutect2 files
grep -v "#" "$outdir/$proband_id.Mutect2.$config.PON.gnomad.vcf" | tr " " "\t" > "$tmp_file"
awk -F '\t' '{print $1, $2-1, $2, $4, $5}' "$tmp_file" > "$$outdir/$proband_id.forPhasing.bed"
awk -v prefix="$BAMprefix" 'BEGIN{OFS="\t"} {$6 = prefix; print}' $outdir/$proband_id.forPhasing.bed > $outdir/$proband_id.MF.$config.phasingInput.bed
rm "$tmp_file"  # Clean up temporary file


# MF2.Run Docker command for ReadLevel Features extraction
docker run -it --rm \
  -v "$git:/MF" \
  -v "$testrun_dir:/mnt/testrun_docker" \
  yanmei/mosaicforecast:0.0.1 \
  /bin/bash -c ' \
  ReadLevel_Features_extraction.py \
  /mnt/testrun_docker/OUTPUT/$proband_id.MF.$config.phasingInput.bed \
  /mnt/testrun_docker/OUTPUT/$proband_id.MF.$config.features.bed \
  /mnt/testrun_docker/BAM \
  /mnt/testrun_docker/Resources/hs37d5.fa \
  /mnt/testrun_docker/Resources/$config/k24.umap.wg.bw \
  2 \
  bam \
  && \
  Prediction.R \
  /mnt/testrun_docker/OUTPUT/$proband_id.MF.$config.features.bed \
  /mnt/testrun_docker/Resources/MF/models_trained/50xRFmodel_addRMSK_Refine.rds \
  Refine \
  /mnt/testrun_docker/OUTPUT/$proband_id.MF.$config.genotype.predictions.refined.bed \
  ' || { echo "Variant calling for MosaicForecast failed."; exit 1; }

