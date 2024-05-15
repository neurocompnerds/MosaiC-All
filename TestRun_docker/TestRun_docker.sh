#!/bin/bash

######specify variables and directories
GIT=/home/neuro/MosaiC-All/MosaiC-All #specify where MosaiC-All is cloned into
CONFIG_FILE=$GIT/config/Mosaic-All.TestRun.config # make sure all directories were specified as required in this config file
TESTRUN_DIR=$GIT/TestRun_docker
BAMDIR=$TESTRUN_DIR/BAM
ProbandBamFile=$BAMDIR/1465_1024-pfc-bulk.rehead.bam
ProbandID=1465_1024-pfc-bulk
Gender=M
OUTDIR=$TESTRUN_DIR/OUTPUT

source $CONFIG_FILE

######Download resources
RESOURCES=${SCRIPTDIR}/TestRun_docker/Resources
mkdir $RESOURCES
wget https://storage.cloud.google.com/gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz -P $RESOURCES
wget https://storage.cloud.google.com/gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz.tbi -P $RESOURCES
wget https://storage.cloud.google.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf -P $RESOURCES
wget https://storage.cloud.google.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf.idx -P $RESOURCES
wget https://storage.cloud.google.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf -P $RESOURCES
wget https://storage.cloud.google.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx -P $RESOURCES
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz -P $RESOURCES
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz.fai -P $RESOURCES


######Installation
docker image pull yanmei/mosaicforecast:0.0.1
docker image pull rborgesm/mosaichunter:1.0
docker image pull broadinstitute/gatk:latest


#####Variant Calling using respective tools

################################################################################     MOSAICHUNTER    #########################################################################


#MH1.prefilter

java -jar $MHDIR/build/mosaichunter.jar -C $MHDIR/conf/exome_parameters.properties \
-P reference_file=$REFGEN \
-P input_file=$ProbandBamFile \
-P heterozygous_filter.sex=$Gender \
-P output_dir=$OUTDIR/$ProbandID.parameters.log

echo "Pre-filter completed for $ProbandID"

#MH2.define Alpha and beta value and remove white spaces before

Al=$(cat $OUTDIR/$ProbandID.parameters.log/stdout*.log | grep "alpha" | cut -d ":" -f2)
Alpha=$(echo "$Al" | sed 's/^ *//g')

Be=$(cat $OUTDIR/$ProbandID.parameters.log/stdout*.log | grep "beta" | cut -d ":" -f2)
Beta=$(echo "$Be" | sed 's/^ *//g')

Dp=$(cat $OUTDIR/$ProbandID.parameters.log/stdout*.log | grep "average depth" | cut -d ":" -f2)
Depth=$(echo "$Dp" | sed 's/^ *//g')

#MH3.execute mosaic variant calling

java -Djava.io.tmpdir=${TMPDIR} -jar $MHDIR/build/mosaichunter.jar -C $MHDIR/conf/exome.properties \
-P reference_file=$REFGEN \
-P input_file=$ProbandBamFile \
-P mosaic_filter.sex=$Gender \
-P depth=$Depth \
-P mosaic_filter.alpha_param=$Alpha \
-P mosaic_filter.beta_param=$Beta \
-P mosaic_filter.dbsnp_file=$DBSNP \
-P repetitive_region_filter.bed_file=$REPEATS \
-P common_site_filter.bed_file=$COMMONERROR \
-P output_dir=$OUTDIR/$ProbandID

echo "## INFO: (MosaicHunter_WES_Singlemode.sh) Somatic variant calling completed for $ProbandID"

#MH4.Process the outputs files

cat $OUTDIR/$ProbandID/final.passed.tsv  > $OUTDIR/$ProbandID.MH.$CONFIG.final.passed.tsv
cat $OUTDIR/$ProbandID/final.passed.tsv | awk '{print $1, $2, $7, $9}' | tr " " "\t" > $OUTDIR/$ProbandID.MH.$CONFIG.forAnnovar.singlemode.vcf

#MH5.log file

grep "input_file =" $OUTDIR/$ProbandID/stdout_*.log > $OUTDIR/$ProbandID.MH.summary.log
tail  $OUTDIR/$ProbandID/stdout_*.log -n16 >> $OUTDIR/$ProbandID.MH.summary.log
cat $OUTDIR/$ProbandID/stdout_*.log >> $OUTDIR/$ProbandID.MH.stdout.log

echo "## INFO: (MosaicHunter_WES_Singlemode.sh) done for $ProbandID"

#MH6.Remove Folders
rm -r $OUTDIR/$ProbandID.parameters.log
rm -r $OUTDIR/$ProbandID



#####################################################################   MUTECT2   #####################################################################
docker image pull broadinstitute/gatk:latest
docker run -v ${your_local_directory}:/GATK --rm -it  broadinstitute/gatk:latest /bin/bash

#Mutect2.1 Variant calling
gatk Mutect2 \
-R $REFGEN \
-I $ProbandBamFile  \
--tumor-sample $ProbandID \
--germline-resource $GERMLINE_RESOURCES \
--panel-of-normals $PON \
--af-of-alleles-not-in-resource -1 \
-O $OUTDIR/$ProbandID.Mutect2.$CONFIG.PON.gnomad.vcf

#Mutect2.2 Tag each variant 

gatk FilterMutectCalls \
-V $OUTDIR/$ProbandID.Mutect2.$CONFIG.PON.gnomad.vcf \
-R $REFGEN \
--min-reads-per-strand 2 \
-O $OUTDIR/$ProbandID.Mutect2.$CONFIG.singlemode.filtered.vcf > $OUTDIR/$ProbandID.Mutect2.$CONFIG.filtered.log 2>&1


#Mutect2.3 Manual Filtration

module load BCFtools/1.17-GCC-11.2.0
bcftools view -O v -f PASS -i '(FORMAT/AD[0:1] >= 5 & FORMAT/DP[0]>=20) && (FORMAT/AF[0:0] <=0.4 || FORMAT/AF[0:0]>=0.7)' $OUTDIR/$ProbandID.Mutect2.$CONFIG.singlemode.filtered.vcf > $OUTDIR/$ProbandID.Mutect2.$CONFIG.singlemode.PASS.aaf.vcf

#####################################################################   MOSAICFORECAST   ################################################################

#MF1.process Mutect2 files
BAMFILE=$(find "$BAMDIR" -type f -name "$ProbandID*.bam")
BAMprefix=$(basename "$BAMFILE" | sed 's/\.[^.]*$//')
tmp_file=$(mktemp)  # Create a temporary file
grep -v "#" $OUTDIR/$ProbandID.Mutect2.$CONFIG.PON.gnomad.vcf | tr " " "\t" > "$tmp_file"
awk -F  '\t' '{print $1, $2-1, $2, $4, $5}' "$tmp_file" > "$OUTDIR/$ProbandID.forPhasing.bed"
rm "$tmp_file"  # Remove the temporary file
awk -v prefix="$BAMprefix" 'BEGIN{OFS="\t"} {$6 = prefix; print}' $OUTDIR/$ProbandID.forPhasing.bed > $OUTDIR/$ProbandID.MF.$CONFIG.phasingInput.bed

#MF2.Extract read features
docker run -v ${your_local_directory}:/MF --rm -it yanmei/mosaicforecast:0.0.1 /bin/bash 
ReadLevel_Features_extraction.py /tmp/$SampleID.MF.$CONFIG.phasingInput.bed /tmp/$SampleID.features.bed /tmp /RefSeq/$REFSEQ /MForecastDir/$CONFIG/k24.umap.wg.bw ${SLURM_NTASKS} bam

#-B ${REFGEN_DIR}:/RefSeq -B $MFORECAST:/MForecastDir -B ${TMPDIR}:/tmp $MFORECAST/mosaicforecast_0.0.1.sif ReadLevel_Features_extraction.py /tmp/$SampleID.MF.$CONFIG.phasingInput.bed /tmp/$SampleID.features.bed /tmp /RefSeq/$REFSEQ /MForecastDir/$CONFIG/k24.umap.wg.bw ${SLURM_NTASKS} bam
#/usr/bin/cp -r ${TMPDIR}/$SampleID.features.bed ${OUTDIR}/$SampleID.MF.$CONFIG.features.bed

#MF3.Genotype prediction
docker run $OUTDIR:/outDir -B $MFORECAST:/MForecastDir $MFORECAST/mosaicforecast_0.0.1.sif Prediction.R /outDir/$SampleID.MF.$CONFIG.features.bed /MForecastDir/models_trained/50xRFmodel_addRMSK_Refine.rds Refine /outDir/$SampleID.MF.$CON>

#-B $OUTDIR:/outDir -B $MFORECAST:/MForecastDir $MFORECAST/mosaicforecast_0.0.1.sif Prediction.R /outDir/$SampleID.MF.$CONFIG.features.bed /MForecastDir/models_trained/50xRFmodel_addRMSK_Refine.rds Refine /outDir/$SampleID.MF.$CONFIG.genotype.predictions.refined.bed

#remove this: this is sample of docker command for MF
#MosaicForecast
docker image pull yanmei/mosaicforecast:0.0.1
docker run -v ${your_local_directory}:/MF --rm -it yanmei/mosaicforecast:0.0.1 /bin/bash
gunzip hs37d5.fa.gz
Phase.py /MF/demo/ /MF/demo/phasing hs37d5.fa /MF/demo/test.input 20 k24.umap.wg.bw 4

