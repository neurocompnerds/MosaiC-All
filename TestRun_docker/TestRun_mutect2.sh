

#####################################################################   MUTECT2   #####################################################################
#!/bin/bash

# Define variables
GIT="/Users/nandinisandran/Desktop/MosaiC-All"
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

# Create output directory
mkdir -p "$OUTDIR"

# Pull Docker image
docker image pull broadinstitute/gatk:4.4.0.0 || { echo "Failed to pull Docker image."; exit 1; }

# Mutect2.1 Variant calling
docker run -it --rm \
  -v "$GIT:/mnt/mosaic-all" \
  -v "$TESTRUN_DIR:/mnt/testrun_docker" \
  broadinstitute/gatk:4.4.0.0 \
  /bin/bash -c "
  gatk Mutect2 \
    -R /mnt/testrun/Resources/hs37d5.fa \
    -I /mnt/testrun/BAM/1465_1024-pfc-bulk.rehead.bam  \
    --tumor-sample $ProbandID \
    --germline-resource $GERMLINE_RESOURCES \
    --panel-of-normals $PON \
    --af-of-alleles-not-in-resource -1 \
    -O $OUTDIR/$ProbandID.Mutect2.$CONFIG.PON.gnomad.vcf
" || { echo "Mutect2.1 Variant calling failed."; exit 1; }

# Mutect2.2 Tag each variant
docker run -it --rm \
  -v "$GIT:/mnt/mosaic-all" \
  -v "$TESTRUN_DIR:/mnt/testrun_docker" \
  broadinstitute/gatk:4.4.0.0 \
  /bin/bash -c "
  gatk TagVariants \
    -V $OUTDIR/$ProbandID.Mutect2.$CONFIG.PON.gnomad.vcf \
    -R /mnt/testrun/Resources/hs37d5.fa \
    --min-reads-per-strand 2 \
    -O $OUTDIR/$ProbandID.Mutect2.$CONFIG.singlemode.filtered.vcf > $OUTDIR/$ProbandID.Mutect2.$CONFIG.filtered.log 2>&1
" || { echo "Mutect2.2 Tag each variant failed."; exit 1; }

# Mutect2.3 Manual Filtration
bcftools view -O v -f PASS -i '(FORMAT/AD[0:1] >= 5 & FORMAT/DP[0]>=20) && (FORMAT/AF[0:0] <=0.4 || FORMAT/AF[0:0]>=0.7)' "$OUTDIR/$ProbandID.Mutect2.$CONFIG.singlemode.filtered.vcf" > "$OUTDIR/$ProbandID.Mutect2.$CONFIG.singlemode.PASS.aaf.vcf" || { echo "Mutect2.3 Manual Filtration failed."; exit 1; }



