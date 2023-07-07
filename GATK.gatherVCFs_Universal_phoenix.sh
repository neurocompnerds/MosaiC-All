#!/bin/bash

#SBATCH -J gatherVCFs
#SBATCH -o /hpcfs/users/%u/Mosaic-All/Log/gatherGVCFs-slurm-%j.out
#SBATCH -p skylake,icelake,skylakehm,v100cpu
#SBATCH -N 1
#SBATCH -n 3
#SBATCH --time=01:30:00
#SBATCH --mem=12GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# A script to merge gVCF files together for later genotyping
## List modules and file paths ##
scriptDir="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/map-n-call"
module purge
module use /apps/skl/modules/all
modList=("Java/17.0.6" "HTSlib/1.17-GCC-11.2.0" "SAMtools/1.17-GCC-11.2.0")

usage()
{
echo "# Script for processing and mapping Illumina 100bp pair-end sequence data and optionally plotting coverage for an interval
# Requires: BWA 0.7.x, Picard, samtools, GATKv4.x, BWA-Picard-GATK-CleanUp.sh.
# This script assumes your sequence files are gzipped
#
# Usage $0 -S sample_name [ -o /path/to/output -c /path/to/config.cfg ] | [ - h | --help ]
#
# Options
# -S	REQUIRED. Sample name if not specified then will be set the same as -p
# -c	OPTIONAL. /path/to/config.cfg. A default config will be used if this is not specified.  The config contains all of the stuff that used to be set in the top part of our scripts
# -o	OPTIONAL. Path to where you want to find your file output (if not specified current directory is used)
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# 
# Original: Derived from Illumina-Phred33-PE-FASTX-BWA-Picard-GATKv2.sh by Mark Corbett, 17/03/2014
# Modified: (Date; Name; Description)
# 24/09/2015; Mark Corbett; Fork original for genomes
# 25/09/2015; Mark Corbett; Pipe to samtools sort; Add getWGSMetrics
# 12/10/2015; Mark Corbett; Fix error collecting .bam files to merge
# 13/05/2016; Mark Corbett; Add option to specify Sample name different from OUTPREFIX.  Make seq file search explicit for *.fastq.gz
# 01/07/2016; Mark Corbett; Improve error handling
# 24/08/2016; Mark Corbett; Fork for HPC version, bring up to date with GATKv3.6
# 18/11/2016; Mark Corbett; Step down number of splits for PrintReads for higher efficiency
# 17/05/2017; Atma Ivancevic; translating for SLURM
# 16/11/2017; Mark Corbett; Fork to split scripts to see if this works on Phoenix; swap Picard for sambamba for duplicate marking
#
"
}

## Set Variables ##
while [ "$1" != "" ]; do
	case $1 in
        -c )            shift
                        Config=$1
                        ;;
        -S )            shift
                        Sample=$1
                        ;;
        -o )            shift
                        workDir=$1
                        ;;
        -h | --help )   usage
                        exit 0
                        ;;
        * )	            usage
                        exit 1
    esac
    shift
done
if [ -z "$Config" ]; then # If no Config file specified use the default
    Config=$scriptDir/configs/BWA-GATKHC.hs38DH_phoenix.cfg
    echo "## INFO: Using the default config ${Config}"
fi
source $Config
if [ -z "$Sample" ]; then # If no Sample name specified then do not proceed
	usage
	echo "## ERROR: You need to specify a Sample name that refers to your .bam file \$Sample.marked.sort.bwa.$BUILD.bam."
	exit 1
fi
if [ -z "$workDir" ]; then # If no output directory then use current directory
	workDir=/hpcfs/users/${USER}/BWA-GATK/$Sample
	echo "## INFO: Using $workDir as the output directory"
fi

tmpDir=/hpcfs/groups/phoenix-hpc-neurogenetics/tmp/${USER}/${Sample} # Use a tmp directory for all of the GATK and samtools temp files
if [ ! -d "$tmpDir" ]; then
	mkdir -p $tmpDir
fi

if [ ! -d "$gVcfFolder" ]; then
	mkdir -p $gVcfFolder
fi
echo "## INFO: All gVCF files will be output to ${gVcfFolder}"

## Load modules ##
for mod in "${modList[@]}"; do
    module load $mod
done

cd $tmpDir
cat $tmpDir/*.$Sample.pipeline.log >> $workDir/$Sample.pipeline.log
find *.$Sample.snps.g.vcf > $tmpDir/$Sample.gvcf.list.txt
sed 's,^,-I '"$tmpDir"'\/,g' $tmpDir/$Sample.gvcf.list.txt > $tmpDir/$Sample.inputGVCF.txt

java -Xmx8g -Djava.io.tmpdir=$tmpDir -jar $GATKPATH/GenomeAnalysisTK.jar GatherVcfs  \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
$(cat $tmpDir/$Sample.inputGVCF.txt) \
-O $gVcfFolder/$Sample.$BUILD.snps.g.vcf >> $workDir/$Sample.pipeline.log  2>&1

bgzip $gVcfFolder/$Sample.$BUILD.snps.g.vcf
tabix $gVcfFolder/$Sample.$BUILD.snps.g.vcf.gz

## Check for bad things and clean up
if [ ! -f "$gVcfFolder/$Sample.$BUILD.snps.g.vcf.gz.tbi" ]; then
    echo "##ERROR: Some bad things went down while this script was running please see $workDir/$Sample.pipeline.ERROR.log and prepare for disappointment."
    exit 1
fi

grep -i ERROR $workDir/$Sample.pipeline.log > $workDir/$Sample.pipeline.ERROR.log
if [ -z $(cat $workDir/$Sample.pipeline.ERROR.log) ]; then
	rm $workDir/$Sample.pipeline.ERROR.log
	rm -r $tmpDir
else 
	echo "##ERROR: Some bad things went down while this script was running please see $workDir/$Sample.pipeline.ERROR.log and prepare for disappointment."
    exit 1
fi
