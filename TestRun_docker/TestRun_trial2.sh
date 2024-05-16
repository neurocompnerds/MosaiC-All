
######Download resources
#RESOURCES=${SCRIPTDIR}/TestRun_docker/Resources
#mkdir $RESOURCES
#wget https://storage.cloud.google.com/gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz -P $RESOURCES
#wget https://storage.cloud.google.com/gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz.tbi -P $RESOURCES
#wget https://storage.cloud.google.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf -P $RESOURCES
#wget https://storage.cloud.google.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf.idx -P $RESOURCES
#wget https://storage.cloud.google.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf -P $RESOURCES
#wget https://storage.cloud.google.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx -P $RESOURCES
#wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz -P $RESOURCES
#wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz.fai -P $RESOURCES

#gunzip $RESOURCES/hs37d5.fa.gz
#samtools faidx $RESOURCES/hs37d5.fa

######Installation
#docker image pull yanmei/mosaicforecast:0.0.1
#docker image pull rborgesm/mosaichunter:1.0
#docker image pull broadinstitute/gatk:latest

#####Variant Calling using respective tools

################################################################################     MOSAICHUNTER    #########################################################################
#!/bin/bash

# Define variables
GIT=/Users/nandinisandran/Desktop/MosaiC-All
TESTRUN_DIR=$GIT/TestRun_docker
BAMDIR=$TESTRUN_DIR/BAM
ProbandBamFile=$BAMDIR/1465_1024-pfc-bulk.rehead.bam
ProbandID=1465_1024-pfc-bulk
Gender=M
OUTDIR=$TESTRUN_DIR/OUTPUT
RESOURCES=$TESTRUN_DIR/Resources

# Create output directory
mkdir -p $OUTDIR

# Pull Docker image
docker image pull rborgesm/mosaichunter:1.0

# Run Docker container and execute commands within the container
docker run -it --rm \
  -v $GIT:/mnt/mosaic-all \
  -v $TESTRUN_DIR:/mnt/testrun \
  rborgesm/mosaichunter:1.0 \
  /bin/bash -c "
  MHDIR=/usr/MosaicHunter
  java -Xmx4G -jar \$MHDIR/build/mosaichunter.jar \
    -C \$MHDIR/conf/exome_parameters.properties \
    -P reference_file=/mnt/testrun/Resources/hs37d5.fa \
    -P input_file=/mnt/testrun/BAM/1465_1024-pfc-bulk.rehead.bam \
    -P heterozygous_filter.sex=$Gender \
    -P output_dir=/mnt/testrun/OUTPUT/$ProbandID.parameters.log
"

echo "Pre-filter completed for $ProbandID"

#MH2.define Alpha and beta value and remove white spaces before

Al=$(cat $OUTDIR/$ProbandID.parameters.log/stdout*.log | grep "alpha" | cut -d ":" -f2)
Alpha=$(echo "$Al" | sed 's/^ *//g')

Be=$(cat $OUTDIR/$ProbandID.parameters.log/stdout*.log | grep "beta" | cut -d ":" -f2)
Beta=$(echo "$Be" | sed 's/^ *//g')

Dp=$(cat $OUTDIR/$ProbandID.parameters.log/stdout*.log | grep "average depth" | cut -d ":" -f2)
Depth=$(echo "$Dp" | sed 's/^ *//g')

#MH3.execute mosaic variant calling

docker run -it --rm \
  -v $GIT:/mnt/mosaic-all \
  -v $TESTRUN_DIR:/mnt/testrun \
  -e COMMONERROR=$COMMONERROR \
  -e REPEATS=$REPEATS \
  rborgesm/mosaichunter:1.0 \
  /bin/bash -c "
  REPEATS=\$REPEATS
  COMMONERROR=\$COMMONERROR
  MHDIR=/usr/MosaicHunter
  java -Xmx4G -jar \$MHDIR/build/mosaichunter.jar \
   -C ./conf/exome.properties \
    -P reference_file=/mnt/testrun/Resources/hs37d5.fa \
    -P input_file=/mnt/testrun/BAM/1465_1024-pfc-bulk.rehead.bam \
    -P mosaic_filter.sex=$Gender \
    -P depth=$Depth \
    -P mosaic_filter.alpha_param=$Alpha \
    -P mosaic_filter.beta_param=$Beta \
    -P mosaic_filter.dbsnp_file=$DBSNP \
    -P repetitive_region_filter.bed_file=\$REPEATS \
    -P common_site_filter.bed_file=\$COMMONERROR \
    -P output_dir=$OUTDIR/$ProbandID
"

echo "## INFO: (MosaicHunter_WES_Singlemode.sh) Somatic variant calling completed for $ProbandID"

