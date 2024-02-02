#!/bin/sh
#SBATCH -J pGoM.sh
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
#This script is designed for MosaiC-All  pipeline
#It requires family VCF in trios
#10-father
#11-mother
#12-child (proband/sib)
#
#SCRIPT EXECUTION
#
#sbatch $0 -v /path/to/familytrio.vcf -s /path/to/ID.list.txt -o /path/to/output/
#
# Options
#-v 	<directory>		REQUIRED: input directory (not the file)
#-s 	<file>    	  	REQUIRED: A file e.g. sampleID.list (one header row and then tab-delimited columns \$BAMdir,\$ProbandID,\$Gender,\$Mother,\$Father)
#-o 	<directory> 		REQUIRED: output directory
"
}

# Load modules
module purge
module use /apps/modules/all
module load BCFtools/1.17-GCC-11.2.0

## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -v )                    shift
                                        VCFDIR=$1
                                        ;;
                -s )                    shift
                                        SAMPLELIST=$1
                                        ;;
                -o )                    shift
                                        OUTDIR=$1
                                        ;;
                * )                     usage
                                        exit 1
        esac
        shift
done


#Check if all information are provided

if [ -z "$VCFDIR" ]; then
    usage
    echo "## ERROR: You need to provide the directory where we can find all the FamilyVCFs listed in the $SAMPLELIST" > $LOGDIR/pGoM.log
    exit 1
fi

if [ -z "$SAMPLELIST" ]; then
    usage
    echo "## ERROR: You need to provide a sample list" 
    echo "#-s REQUIRED sampleID.list (one header row and then tab-delimited columns \$BAMdir,\$ProbandID,\$Gender,\$Mother,\$Father,\$FamilyVCFID)" >> $LOGDIR/pGoM.log
	exit 1
fi

if [ ! -d "${OUTDIR}" ]; then
    mkdir -p ${OUTDIR}
	echo "## INFO: output directory created, you'll find all of the outputs and log files in here: ${OUTDIR}" >> $LOGDIR/pGoM.log
fi

#Array from a list of Samples (ignoring the header of the file)
while IFS= read -r line; do
    if [[ ! "$line" =~ ^#.*$ ]]; then
        SAMPLEID+=("$line")
    fi
done < "$SAMPLELIST"

# Iteration starts here
for SAMPLEID in "${SAMPLEID[@]}"; do

    #Defining variables from each row

    	ProbandID=$(awk '{print $2}' <<< "$SAMPLEID ")
   	MotherID=$(awk '{print $4}' <<< "$SAMPLEID ")
   	FatherID=$(awk '{print $5}' <<< "$SAMPLEID ")
    	FamilyVCF=$(awk '{print $6}' <<< "$SAMPLEID ")
     	VCF=$VCFDIR/$FamilyVCF

     	echo "pGoM Pipeline for $ProbandID, $MotherID, $FatherID based on familyvcf:$VCF" >> $LOGDIR/$ProbandID.pGoM.pipeline.log

     	#Check if the VCF is available in the specified directory or else dont proceed
       	if [ -s "$VCF" ]; then
    	echo "$VCF exists and is not empty." >> $LOGDIR/$ProbandID.pGoM.pipeline.log
	else
    	echo "$VCF DOES NOT exist or is empty." >> $LOGDIR/$ProbandID.pGoM.pipeline.log
	fi

	#Find the numeric index of each sample so that  filtering can be done
	# Get numeric indices for each sample and convert to 0-based
	col_P=$(($(bcftools query -l "$VCF" | grep -n "004P" | cut -d: -f1) - 1))
	col_M=$(($(bcftools query -l "$VCF" | grep -n "004M" | cut -d: -f1) - 1))
	col_F=$(($(bcftools query -l "$VCF" | grep -n "004F" | cut -d: -f1) - 1))

	# Display the 0-based numeric indices
	echo "Numeric index for 004P: $col_P"
	echo "Numeric index for 004M: $col_M"
	echo "Numeric index for 004F: $col_F"

	#Filtering for father
	
	##Filter 1: Genotype: All inherited variants from father
	bcftools view -O v -i "GT[$col_F] = '0/1' && GT[$col_P] = '0/1' && GT[$col_M] = '0/0'" "$VCF" >>  $OUTDIR/$ProbandID.fa.all.inherited.vcf
	##Filter 2: DP20: All samples must have at least 20 supporting reads for each variant
	bcftools view -O v -f PASS -i "(FORMAT/DP[$col_F] >= 20) && (FORMAT/DP[$col_P] >= 20) && (FORMAT/DP[$col_M] >= 20)" $OUTDIR/$ProbandID.fa.all.inherited.vcf > $OUTDIR/$ProbandID.fa.DP20.vcf
	##Filter 3: AAF (father)
	bcftools view -O v -f PASS -i "(AD[$col_F:1]/FORMAT/DP[$col_F])<=0.4 | (AD[$col_F:1]/FORMAT/DP[$col_F])>=0.7" "$OUTDIR/$ProbandID.fa.DP20.vcf" > "$OUTDIR/$ProbandID.fa.mosaic.DP20.vcf"
	##Filter 4: AAF (proband)
	bcftools view -O v -f PASS -i "(AD[$col_P:1]/FORMAT/DP[$col_P])>=0.4 && (AD[$col_P:1]/FORMAT/DP[$col_P])<=0.7" "$OUTDIR/$ProbandID.fa.mosaic.DP20.vcf" > $OUTDIR/$ProbandID.fa.final.vcf

	#COUNTS
	V=$( bcftools view $VCF | grep -v "^#" | wc -l)
	A=$(grep -v "^#" $OUTDIR/$ProbandID.fa.all.inherited.vcf | wc -l)
	B=$(grep -v "^#" $OUTDIR/$ProbandID.fa.DP20.vcf| wc -l)
	C=$(grep -v "^#" $OUTDIR/$ProbandID.fa.mosaic.DP20.vcf | wc -l)
	D=$(grep -v "^#" $OUTDIR/$ProbandID.fa.final.vcf  | wc -l)
					
	echo -e "$ProbandID\tfather\tTrios\t$V\t$A\t$B\t$C\t$D" | tr " " "\t" >> $OUTDIR/mosaic.father.trios.counts.txt

	#Filtering for mother
 	
	##Filter 1: Genotype: All inherited variants from mother
	bcftools view -O v -i "GT[$col_F] = '0/0' && GT[$col_P] = '0/1' && GT[$col_M] = '0/1'" "$VCF" >>  "$OUTDIR/$ProbandID.mo.all.inherited.vcf"
	##Filter 2: DP20: All samples must have at least 20 supporting reads for each variant
	bcftools view -O v -f PASS -i "(FORMAT/DP[$col_F] >= 20) && (FORMAT/DP[$col_P] >= 20) && (FORMAT/DP[$col_M] >= 20)" "$OUTDIR/$ProbandID.mo.all.inherited.vcf" > "$OUTDIR/$ProbandID.mo.DP20.vcf"
	##Filter 3: AAF (mother)
	bcftools view -O v -f PASS -i "(AD[$col_M:1]/FORMAT/DP[$col_M])<=0.4 | (AD[$col_M:1]/FORMAT/DP[$col_M])>=0.7" "$OUTDIR/$ProbandID.mo.DP20.vcf" > "$OUTDIR/$ProbandID.mo.mosaic.DP20.vcf"
	##Filter 4: AAF (proband)
	bcftools view -O v -f PASS -i "(AD[$col_P:1]/FORMAT/DP[$col_P])>=0.4 && (AD[$col_P:1]/FORMAT/DP[$col_P])<=0.7" "$OUTDIR/$ProbandID.mo.mosaic.DP20.vcf" > "$OUTDIR/$ProbandID.fa.final.vcf"

	#COUNTS
	V=$( bcftools view $VCF | grep -v "^#" | wc -l)
	E=$(grep -v "^#" $OUTDIR/$ProbandID.mo.all.inherited.vcf | wc -l)
	F=$(grep -v "^#" $OUTDIR/$ProbandID.mo.DP20.vcf| wc -l)
	G=$(grep -v "^#" $OUTDIR/$ProbandID.mo.mosaic.DP20.vcf | wc -l)
	H=$(grep -v "^#" $OUTDIR/$ProbandID.mo.final.vcf  | wc -l)

 	echo -e "$ProbandID\tmother\tTrios\t$V\t$E\t$F\t$G\t$H" | tr " " "\t" >> $OUTDIR/mosaic.mother.trios.counts.txt

 	echo "pGoM Pipeline for $ProbandID is completed" >> $LOGDIR/$ProbandID.pGoM.pipeline.log

  	##Tidy up for postprocessing in R by removing the headers and combining both parents (of all samples) into a single file
	awk -v ID="$FatherID" '$0 !~ /^##/ {print ID "\t" "GATKHC_filtered" "\t" $0}' $OUTDIR/$ProbandID.fa.final.vcf >> $OUTDIR/pGoM_prefilter_List.txt   
	awk -v ID="$MotherID" '$0 !~ /^##/ {print ID "\t" "GATKHC_filtered" "\t" $0}' $OUTDIR/$ProbandID.mo.final.vcf >> $OUTDIR/pGoM_prefilter_List.txt
done
