#!/bin/bash

#SBATCH -J GATKHC
#SBATCH -o /hpcfs/users/%u/log/GATK4HC-slurm-%j.out
#SBATCH -p skylake,icelake,skylakehm,v100cpu
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=07:00:00
#SBATCH --mem=8GB

# Notification Configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# A script to call variants using the GATK v4.x best practices designed for the Phoenix supercomputer but will work on stand alone machines too

## List modules and file paths ##
scriptDir="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/map-n-call"
module purge
module use /apps/skl/modules/all
modList=("Java/17.0.6")

usage()
{
echo "# A script to call variants using the GATK v4.x best practices
# Requires: GATKv4.x
#
# Usage: sbatch --array 0-23 $0 -S sample_name [ -o /path/to/output -c /path/to/config.cfg ] | [ - h | --help ]
#
# Options
# -S	REQUIRED. Sample name if not specified then will be set the same as -p
# -c	OPTIONAL. /path/to/Config.cfg. A default Config will be used if this is not specified.  The Config contains all of the stuff that used to be set in the top part of our scripts
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

## Define the array for the batch job ##
bedFile=($arrIndexBedFiles)
if [ ! -d "$tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}" ]; then
    mkdir -p $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}
fi

## Load modules ##
for mod in "${modList[@]}"; do
    module load $mod
done

## Start of the script ##
 
cd $tmpDir
java -Xmx6g -Djava.io.tmpdir=$tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]} -jar $GATKPATH/GenomeAnalysisTK.jar HaplotypeCaller \
-I $workDir/$Sample.recal.sorted.bwa.$BUILD.bam \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
-L $ChrIndexPath/${bedFile[$SLURM_ARRAY_TASK_ID]} \
--dbsnp $GATKREFPATH/$BUILD/$DBSNP \
--min-base-quality-score 20 \
--native-pair-hmm-threads 10 \
-ERC GVCF \
-O $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.$Sample.snps.g.vcf >> $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.$Sample.pipeline.log 2>&1
