#!/bin/sh
#SBATCH -J CombineCalls
#SBATCH -o /home/%u/MosaiC-All/Log/CC-slurm-%j.out

#SBATCH -p skylake,icelake
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=00:05:00
#SBATCH --mem=1GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=%u@adelaide.edu.au


usage()
{
echo "
#
#README
#This script is designed for Mosaic-All wrapper script,
#to COMBINE all variant calls (of all samples executed by batch) in one meta-file
#sbatch $0 -s ID -v MH -f FILE -o /path/to/output
#-s <string>          REQUIRED: The ID of the Proband
#-v <string>          REQUIRED: The variant caller to use, typically MH
#-f /path/to/file     REQUIRED: The *.final.passed.tsv output from the variant caller   					
#-o /path/to/output/  REQUIRED: Path to the output directory
#
"
}


## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -s )                    shift
                                        ProbandID=$1
                                        ;;
                -v )                    shift
                                        VC=$1
                                        ;;
                -f )                    shift
                                        FILE=$1
                                        ;;
                -o )                    shift
                                        OUTPUT=$1
                                        ;;
                * )                     usage
                                        exit 1
        esac
        shift
done

## Check for all required arguments
for ARG in $ProbandID $VC $FILE $OUTPUT; do
    if [ -z "${ARG}" ]; then
        usage
        echo "## ERROR: Please check that you have supplied all required argunments for this script"
        exit 1
    fi
done

awk -v ID="$ProbandID" '$0 !~ /^##/ {print ID "\t" "$VC" "\t" $0}' $FILE >> $OUTPUT/$VC.calls.txt
