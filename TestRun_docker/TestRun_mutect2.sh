

#####################################################################   MUTECT2   #####################################################################
#!/bin/bash

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
#docker image pull broadinstitute/gatk:4.4.0.0 || { echo "Failed to pull Docker image."; exit 1; }

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

exit

# Mutect2.3 Manual Filtration
bcftools view -O v -f PASS -i '(FORMAT/AD[0:1] >= 5 & FORMAT/DP[0]>=20) && (FORMAT/AF[0:0] <=0.4 || FORMAT/AF[0:0]>=0.7)' "$OUTDIR/$ProbandID.Mutect2.$CONFIG.singlemode.filtered.vcf" > "$OUTDIR/$ProbandID.Mutect2.$CONFIG.singlemode.PASS.aaf.vcf" || { echo "Mutect2.3 Manual Filtration failed."; exit 1; }



