#!/bin/bash

usage()
{
echo "Executing TestRun_docker.sh
#
# Usage $0 -d /path/to/MosaiC-All -s 1465_1024-pfc-bulk.rehead.bam -i 1465_1024-pfc-bulk -g M [ - h | --help ]
#
# Options
# -d git              REQUIRED: path to where the git is cloned into
# -s proband_bam      REQUIRED: full name of bamfile (which is stored in $git/BAM)
# -i proband_id       REQUIRED: proband id for output files (can be the same as previous -s option)
# -g gender           REQUIRED: gender
# -h or --help  Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
#
# See: https://github.com/neurocompnerds/MosaiC-All for history and new updates.
#
"
}

## Set Variables ##
while [ "$1" != "" ]; do
    case $1 in
            -d )                    shift
                                    git_dir=$1
                                    ;;
            -s )                    shift
                                    proband_bam=$1
                                    ;;
            -i )                    shift
                                    proband_id=$1
                                    ;;
            -g )                    shift
                                    gender=$1
                                    ;;
            -h | --help )           usage
                                    exit 0
                                    ;;
            *  )                    usage
                                    exit 1
    esac
    shift
done

# Define variables
testrun_dir="$git_dir/TestRun_docker"
bam_dir="$testrun_dir/BAM"
config="hg19"
outdir="$testrun_dir/OUTPUT"
resources="$testrun_dir/Resources"

if [ ! -d "${outdir}" ]; then
    mkdir -p ${outdir}
fi

if [ ! -d "${resources}" ]; then
    mkdir -p ${resources}
fi

# check and download samtools and bcftools modules 

if ! command -v samtools &> /dev/null
then
    echo "Samtools not found. Downloading and installing..."
    wget "https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2" -P $resources
    tar -vxjf $resources/samtools-1.9.tar.bz2
    cd $resources/samtools-1.9
    make
    echo "Samtools has been installed."
else
    echo "Samtools is already installed."
    samtools --version
fi

if ! command -v bcftools &> /dev/null
then
    echo "BCFtools not found. Downloading and installing..."
    wget "https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2" -P $resources
    tar -vxjf $resources/bcftools-1.9.tar.bz2
    cd $resources/bcftools-1.9
    make
    echo "BCFtools has been installed."
else
    echo "BCFtools is already installed."
    bcftools --version
fi

# download, extract and create index for reference genome
wget "https://storage.googleapis.com/genomics-public-data/references/hs37d5/hs37d5.fa.gz" -P "$resources"

if [ -e "$resources/hs37d5.fa.gz" ]; then
    # Extract the tar.gz file
    gunzip "$resources/hs37d5.fa.gz"
    samtools faidx "$resources/hs37d5.fa" -o "$resources/hs37d5.fa.fai"
fi


############################################################################ MOSAICHUNTER ####################################################################################################

# Resources for MosaicHunter

wget "https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz" -P "$resources"

# Check if dbsnp_138.b37.vcf.gz  was downloaded successfully
if [ -e "$resources/dbsnp_138.b37.vcf.gz" ]; then
    # Extract the tar.gz file
    gunzip "$resources/dbsnp_138.b37.vcf.gz"
fi

dbsnp="dbsnp_138.b37.vcf"
repeats="all_repeats.b37.bed" # this can be found in the docker itself
common_error="WES_Agilent_71M.error_prone.b37.bed" # this can be found in the docker itself

# Pull Docker image
docker image pull rborgesm/mosaichunter:1.0

# Run Docker container for pre-filtering
docker run -it --rm \
  -v "$git_dir:/mnt/mosaic-all" \
  -v "$testrun_dir:/mnt/testrun_docker" \
  rborgesm/mosaichunter:1.0 \
  /bin/bash -c "
  MHDIR=/usr/MosaicHunter
  /usr/local/openjdk-8/bin/java -Xmx4G -jar \$MHDIR/build/mosaichunter.jar \
    -C \$MHDIR/conf/exome_parameters.properties \
    -P reference_file=/mnt/testrun_docker/Resources/hs37d5.fa \
    -P input_file=/mnt/testrun_docker/BAM/$proband_bam \
    -P heterozygous_filter.sex=$gender \
    -P output_dir=/mnt/testrun_docker/OUTPUT/$proband_id.parameters.log
" || { echo "MosaicHunter prefiltering failed."; exit 1; }

# Extract Alpha, Beta, and Depth values from the log file
alpha=$(grep "alpha" $outdir/$proband_id.parameters.log/stdout*.log | cut -d ":" -f2 | sed 's/^ *//g')
beta=$(grep "beta" $outdir/$proband_id.parameters.log/stdout*.log | cut -d ":" -f2 | sed 's/^ *//g')
depth=$(grep "average depth" $outdir/$proband_id.parameters.log/stdout*.log | cut -d ":" -f2 | sed 's/^ *//g')

# Run Docker container for variant calling
docker run -it --rm \
  -v "$git_dir:/mnt/mosaic-all" \
  -v "$testrun_dir:/mnt/testrun_docker" \
  -e COMMONERROR="$common_error" \
  -e REPEATS="$repeats" \
  rborgesm/mosaichunter:1.0 \
  /bin/bash -c "
  MHDIR=/usr/MosaicHunter
  /usr/local/openjdk-8/bin/java -Xmx4G -jar \$MHDIR/build/mosaichunter.jar \
    -C \$MHDIR/conf/exome.properties \
    -P reference_file=/mnt/testrun_docker/Resources/hs37d5.fa \
    -P input_file=/mnt/testrun_docker/BAM/$proband_bam \
    -P mosaic_filter.sex=$gender \
    -P depth=$depth \
    -P mosaic_filter.alpha_param=$alpha \
    -P mosaic_filter.beta_param=$beta \
    -P mosaic_filter.dbsnp_file=/mnt/testrun_docker/Resources/$dbsnp \
    -P repetitive_region_filter.bed_file=\$MHDIR/resources/\$REPEATS \
    -P common_site_filter.bed_file=\$MHDIR/resources/\$COMMONERROR \
    -P output_dir=/mnt/testrun_docker/OUTPUT/$proband_id
" || { echo "MosaicHunter variant calling failed."; exit 1; }

# Process the output files
cat "$outdir/$proband_id/final.passed.tsv" > "$outdir/$proband_id.MH.$config.final.passed.tsv"
awk '{print $1, $2, $7, $9}' "$outdir/$proband_id/final.passed.tsv" | tr " " "\t" > "$outdir/$proband_id.MH.$config.forAnnovar.singlemode.vcf"

# Generate log file
grep "input_file =" $outdir/$proband_id/stdout_*.log > "$outdir/$proband_id.MH.summary.log"
tail -n 16 $outdir/$proband_id/stdout_*.log >> "$outdir/$proband_id.MH.summary.log"
cat $outdir/$proband_id/stdout_*.log >> "$outdir/$proband_id.MH.stdout.log"

#####################################################################   MUTECT2   #####################################################################

# define resources for Mutect2
pon="Mutect2-exome-panel.vcf"
germline_resources="af-only-gnomad.raw.sites.vcf"

# Download resources for Mutect2

wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf -P $resources
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf.idx -P $resources
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf -P $resources
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx -P $resources

# Pull Docker image
docker image pull broadinstitute/gatk:4.4.0.0 || { echo "Failed to pull Docker image."; exit 1; }

# create dict for hs37dH.fa
docker run -it --rm \
  -v "$testrun_dir:/mnt/testrun_docker" \
  broadinstitute/gatk:4.4.0.0 \
  /bin/bash -c "
  gatk CreateSequenceDictionary \
  -R /mnt/testrun_docker/Resources/hs37d5.fa
"

# Mutect2.1 Variant calling
docker run -it --rm \
  -v "$testrun_dir:/mnt/testrun_docker" \
  broadinstitute/gatk:4.4.0.0 \
  /bin/bash -c "
  gatk Mutect2 \
    -R /mnt/testrun_docker/Resources/hs37d5.fa \
    -I /mnt/testrun_docker/BAM/$proband_bam  \
    --tumor-sample $proband_id \
    --germline-resource /mnt/testrun_docker/Resources/$germline_resources \
    --panel-of-normals /mnt/testrun_docker/Resources/$pon \
    --af-of-alleles-not-in-resource -1 \
    -O /mnt/testrun_docker/OUTPUT/$proband_id.Mutect2.$config.PON.gnomad.vcf
" || { echo "Mutect2.1 Variant calling failed."; exit 1; }

# Mutect2.2 Tag each variant

docker run -it --rm \
  -v "$testrun_dir:/mnt/testrun_docker" \
  broadinstitute/gatk:4.4.0.0 \
  /bin/bash -c "
  gatk FilterMutectCalls \
    -V /mnt/testrun_docker/OUTPUT/$proband_id.Mutect2.$config.PON.gnomad.vcf \
    -R /mnt/testrun_docker/Resources/hs37d5.fa \
    --min-reads-per-strand 2 \
    -O /mnt/testrun_docker/OUTPUT/$proband_id.Mutect2.$config.singlemode.filtered.vcf > /mnt/testrun_docker/OUTPUT/$proband_id.Mutect2.$config.filtered.log 2>&1
" || { echo "Mutect2.2 Tag each variant failed."; exit 1; }



# Mutect2.3 Manual Filtration
bcftools view -O v -f PASS -i '(FORMAT/AD[0:1] >= 5 & FORMAT/DP[0]>=20) && (FORMAT/AF[0:0] <=0.4 || FORMAT/AF[0:0]>=0.7)' "$outdir/$proband_id.Mutect2.$config.singlemode.filtered.vcf" > "$outdir/$proband_id.Mutect2.$config.singlemode.PASS.aaf.vcf" || { echo "Mutect2.3 Manual Filtration failed."; exit 1; }



#####################################################################   MOSAICFORECAST   ################################################################
#!/bin/bash

# Ensure Docker image was pulled successfully
docker image pull yanmei/mosaicforecast:0.0.1 || { echo "Failed to pull MosaicForecast Docker image."; exit 1; }

# download resources for MosaicForecast run

# 1/2 "model_trained" folder from MF which cannot be found in docker
mkdir $resources/MFgit
git clone --depth 1 --branch master https://github.com/parklab/MosaicForecast.git --single-branch $resources/MFgit

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
  -v "$resources:/mnt/resources" \
  yanmei/mosaicforecast:0.0.1 \
  /bin/bash -c '
    cd "/mnt/resources/hg19"
    # Fetch chromosome sizes
    fetchChromSizes hg19 > hg19.chrom.sizes
    # Decompress the Umap score file
    zcat /mnt/resources/hg19/k24.umap.wg.gz > /mnt/resources/hg19/k24.umap.wg
    # Convert wig to BigWig format
    wigToBigWig /mnt/resources/hg19/k24.umap.wg /mnt/resources/hg19/hg19.chrom.sizes /mnt/resources/hg19/k24.umap.wg.bw
' || { echo "Processing Umap score for MosaicForecast failed."; exit 1; }


# MF1.Prepare Input files
# Assign BAM prefix
BAMprefix=$(basename "$bam_dir/$proband_bam" | sed 's/\.[^.]*$//')
echo "$BAMprefix"

# Create temporary file
tmp_file=$(mktemp) || { echo "Failed to create temporary file."; exit 1; }

# Process Mutect2 files
grep -v "#" "$outdir/$proband_id.Mutect2.$config.PON.gnomad.vcf" | tr " " "\t" > "$tmp_file"
awk -F '\t' '{print $1, $2-1, $2, $4, $5}' "$tmp_file" > "$outdir/$proband_id.forPhasing.bed"
awk -v prefix="$BAMprefix" 'BEGIN{OFS="\t"} {$6 = prefix; print}' $outdir/$proband_id.forPhasing.bed > $outdir/$proband_id.MF.$config.phasingInput.bed
rm "$tmp_file"  # Clean up temporary file


echo "$outdir/$proband_id.MF.$config.phasingInput.bed"
ls -l "$outdir/$proband_id.MF.$config.phasingInput.bed"


# MF2.Run Docker command for ReadLevel Features extraction
docker run -it --rm \
  -v "$git_dir:/MF" \
  -v "$testrun_dir:/mnt/testrun_docker" \
  yanmei/mosaicforecast:0.0.1 \
  /bin/bash -c " \
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
  /mnt/testrun_docker/Resources/MFgit/models_trained/50xRFmodel_addRMSK_Refine.rds \
  Refine \
  /mnt/testrun_docker/OUTPUT/$proband_id.MF.$config.genotype.predictions.refined.bed \
" || { echo "Variant calling for MosaicForecast failed."; exit 1; }



