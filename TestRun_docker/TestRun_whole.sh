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
wget "https://storage.googleapis.com/genomics-public-data/references/hs37d5/hs37d5.fa.gz" -P "$resources"

if [ -e "$resources/hs37d5.fa.gz" ]; then
    # Extract the tar.gz file
    gunzip "$resources/hs37d5.fa.gz" -c "resources"
    samtools faidx "hs37d5.fa.gz"
fi


############################################################################ MOSAICHUNTER ####################################################################################################

# Resources for MosaicHunter

wget "https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz" -P "$resources"

# Check if dbsnp_138.b37.vcf.gz  was downloaded successfully
if [ -e "$resources/dbsnp_138.b37.vcf.gz" ]; then
    # Extract the tar.gz file
    gunzip "$resources/dbsnp_138.b37.vcf.gz" -c "resources"
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
  -v "$refgen:/mnt/refgen" \
  rborgesm/mosaichunter:1.0 \
  /bin/bash -c "
  MHDIR=/usr/MosaicHunter
  java -Xmx4G -jar \$MHDIR/build/mosaichunter.jar \
    -C \$MHDIR/conf/exome_parameters.properties \
    -P reference_file=/mnt/refgen/hs37d5.fa \
    -P input_file=/mnt/testrun_docker/BAM/1465_1024-pfc-bulk.rehead.bam \
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
  -v "$refgen:/mnt/refgen" \
  -e COMMONERROR="$common_error" \
  -e REPEATS="$repeats" \
  rborgesm/mosaichunter:1.0 \
  /bin/bash -c "
  MHDIR=/usr/MosaicHunter
  java -Xmx4G -jar \$MHDIR/build/mosaichunter.jar \
    -C \$MHDIR/conf/exome.properties \
    -P reference_file=/mnt/refgen/hs37d5.fa \
    -P input_file=/mnt/testrun_docker/BAM/1465_1024-pfc-bulk.rehead.bam \
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

# Download resources fir Mutect2

wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf -P $resources
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf.idx -P $resources
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf -P $resources
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx -P $resources

# Pull Docker image
docker image pull broadinstitute/gatk:4.4.0.0 || { echo "Failed to pull Docker image."; exit 1; }

# Mutect2.1 Variant calling
docker run -it --rm \
  -v "$git:/mnt/mosaic-all" \
  -v "$testrun_dir:/mnt/testrun_docker" \
  broadinstitute/gatk:4.4.0.0 \
  /bin/bash -c "
  gatk Mutect2 \
    -R /mnt/testrun_docker/Resources/hs37d5.fa \
    -I /mnt/testrun_docker/BAM/1465_1024-pfc-bulk.rehead.bam  \
    --tumor-sample $proband_id \
    --germline-resource /mnt/testrun_docker/Resources/$germline_resources \
    --panel-of-normals /mnt/testrun_docker/Resources/$pon \
    --af-of-alleles-not-in-resource -1 \
    -O /mnt/testrun_docker/OUTPUT/$proband_id.Mutect2.$config.PON.gnomad.vcf
" || { echo "Mutect2.1 Variant calling failed."; exit 1; }

# Mutect2.2 Tag each variant

docker run -it --rm \
  -v "$git:/mnt/mosaic-all" \
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



