#!/bin/sh
#SBATCH -J Mosaicparents.sh
#SBATCH -o /home/%u/MosaiC-All/Log/pGoM.slurm-%j.out
#SBATCH -p skylake,icelake
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --time=12:30:00
#SBATCH --mem=30GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=%u@adelaide.edu.au

## Hard coded paths you may wish to change for your system
LOGDIR="/home/${USER}/MosaiC-All/Log"

usage()
{
echo "
#
#README
#This script is designed for Mosaic-All  pipeline
#It requires family VCF in trios,  in the order of
#10-father
#11-mother
#12-child (proband/sib)
#
#SCRIPT EXECUTION
#
#sbatch $0 -v /path/to/vcf_files/ -s /path/to/ID.list.txt -o /path/to/output/
#
# Options
#-v REQUIRED /path/to/vcf_files/  (directory for family VCFs)
#-s REQUIRED /path/to/ID.list.txt (familyID.txt; $FAMILYID*.vcf)
#-o REQUIRED /path/to/output/     (output dir) 
"
}

# Load modules
module purge
module use /apps/Load modules/all
module load BCFtools/1.17-GCC-11.2.0

## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -v )                    shift
                                        VCFDIR=$1
                                        ;;
                -s )                    shift
                                        FAMILYID=$1
                                        ;;
                -o )                    shift
                                        OUTDIR=$1
                                        ;;
                * )                     usage
                                        exit 1
        esac
        shift
done

## Check for all required arguments
if [ -z "$FAMILYID" ]; then
    usage
    echo "## ERROR: You need to provide a list of familyIDs that identifies the VCFs: $FAMILYID*.vcf "
	exit 1
fi

readarray -t sample < $FAMILYID

#Defined log files
if [ ! -d "$LOGDIR" ]; then
    mkdir -p $LOGDIR
	echo "## INFO: Slurm log files will be placed in this location $LOGDIR" >> $LOGDIR/$FAMILYID.pGoM.log
fi

#Check if all information are provided

if [ -z "$VCFDIR" ]; then
    usage
    echo "## ERROR: You need to provide the directory where we can find the familyVCFs" >> $LOGDIR/$FAMILYID.pGoM.log
    exit 1
fi

if [ ! -d "${OUTDIR}" ]; then
    mkdir -p ${OUTDIR}
	echo "## INFO: output directory created, you'll find all of the outputs and log files in here: ${OUTDIR}" >> $LOGDIR/$FAMILYID.pGoM.log
fi

##Inspect VCFs
VCF=$VCFDIR/${sample[$SLURM_ARRAY_TASK_ID]}*.vcf

#Check if at least 12 columns are present
#if [ "${#vcf_lines[@]}" -lt 12 ]; then
#    echo "The VCF file does not have enough columns." > $LOGDIR/FAMILYID.pGoM.log
#    exit 1
#fi

# Get the headers from the first line, and define
#headers=(${VCF[0]})
#FATHER="${headers[10]}"
#MOTHER="${headers[11]}"
#CHILD="${headers[12]}"

#Comment:
#Attempt to combine both parents, but the VCFs order will disturb.
##################################################################   FATHER    ####################################################################################################################
#if [[ -n "$FATHER"  && -n "$PROBAND"]]; then

    		#fa.all.inherited (just to check if tallies with previous count)
			bcftools view $VCF |grep "^#"  > $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.fa.all.inherited.vcf
			bcftools view $VCF |awk '$12 ~ /^0\/1/ && $10 ~ /^0\/1/ && $11 ~/^0\/0/ {OFS="\t"; print}' >> $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.fa.all.inherited.vcf

			#DP20 in father
			bcftools view -O v -f PASS -i 'FORMAT/DP[0]>=20' $VCF > $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.fa.DP20.vcf

			#aaf filter in father
			bcftools view -O v -f PASS -i '(AD[0:1]/FORMAT/DP[0])<=0.3 | (AD[0:1]/FORMAT/DP[0])>=0.7' $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.fa.DP20.vcf > $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.fa.mosaic.DP20.vcf
			VCF2fa=$OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.fa.mosaic.DP20.vcf

			#GT filter
			bcftools view $VCF2fa |grep "^#" > $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.p1.fa.mosaic.DP20.vcf
			bcftools view $VCF2fa | awk '$10 ~ /^0\/1/ && $11 ~ /^0\/0/ && $12 ~/^0\/1/ {OFS="\t"; print}' >> $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.fa.mosaic.DP20.vcf

			#make sure children is within 0.3-0.7
			bcftools view -O v -f PASS -i '(AD[2:1]/FORMAT/DP[2])>=0.3 && (AD[2:1]/FORMAT/DP[2])<=0.7' $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.fa.mosaic.DP20.vcf > $OUTDIR/$CHILD.$FATHER.fa.final.vcf

			#COUNTS

			V=$( bcftools view $VCF | grep -v "^#" | wc -l)
			A=$(grep -v "^#" $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.fa.all.inherited.vcf | wc -l)
			B=$(grep -v "^#" $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.fa.DP20.vcf| wc -l)
			C=$(grep -v "^#" $VCF2fa | wc -l)
			D=$(grep -v "^#" $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.fa.mosaic.DP20.vcf | wc -l)
			E=$(grep -v "^#" $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.fa.final.vcf  | wc -l)

						
			echo -e "${sample[$SLURM_ARRAY_TASK_ID]}\tfather\tTrios\t$V\t$A\t$B\t$C\t$D\t$E" | tr " " "\t" >> $OUTDIR/mosaic.father.trios.counts.txt

#fi

############################################################## MOTHER    ####################################################################################################################

#if [[ -n "$MOTHER" && -n "$PROBAND"]]; then

			#mo.all.inherited (just to check if tallies with previous count)
			bcftools view $VCF |grep "^#"  > $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.mo.all.inherited.vcf
			bcftools view $VCF |awk '$12 ~ /^0\/1/ && $11 ~ /^0\/1/ && $10 ~/^0\/0/ {OFS="\t"; print}' >> $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.mo.all.inherited.vcf

			#DP20 in mother
			bcftools view -O v -f PASS -i 'FORMAT/DP[1]>=20' $VCF > $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.mo.DP20.vcf

			#aaf filter in mother
			bcftools view -O v -f PASS -i '(AD[1:1]/FORMAT/DP[1])<=0.3 | (AD[1:1]/FORMAT/DP[1])>=0.7' $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.mo.DP20.vcf > $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.mo.mosaic.DP20.vcf
			VCF2mo=$OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.mo.mosaic.DP20.vcf

			#GT filter
			bcftools view $VCF2mo |	grep "^#" > $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.p1.mo.mosaic.DP20.vcf
			bcftools view $VCF2mo |	awk '$10 ~ /^0\/0/ && $11 ~ /^0\/1/ && $12 ~/^0\/1/ {OFS="\t"; print}' >> $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.mo.mosaic.DP20.vcf

			#make sure children is within 0.3-0.7
			bcftools view -O v -f PASS -i '(AD[2:1]/FORMAT/DP[2])>=0.3 && (AD[2:1]/FORMAT/DP[2])<=0.7' $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.mo.mosaic.DP20.vcf > $OUTDIR/$CHILD.$MOTHER.mo.final.vcf

			#COUNTS
			V=$( bcftools view $VCF | grep -v "^#" | wc -l)
			A=$(grep -v "^#" $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.mo.all.inherited.vcf | wc -l)
			B=$(grep -v "^#" $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.mo.DP20.vcf| wc -l)
			C=$(grep -v "^#" $VCF2mo | wc -l)
			F=$(grep -v "^#" $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.mo.mosaic.DP20.vcf | wc -l)
			G=$(grep -v "^#" $OUTDIR/${sample[$SLURM_ARRAY_TASK_ID]}.mo.final.vcf  | wc -l)

						
			echo -e "${sample[$SLURM_ARRAY_TASK_ID]}\tmother\tTrios\t$V\t$A\t$B\t$C\t$D\t$E" | tr " " "\t" >> $OUTDIR/mosaic.mother.trios.counts.txt

#fi
