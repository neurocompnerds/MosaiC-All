#!/bin/sh
#SBATCH -J CombineCalls
#SBATCH -o /hpcfs/users/%u/MosaiC-All/Log/CC-slurm-%j.out
#SBATCH -A robinson
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
#sbatch $SCRIPTDIR/Mutect2.singlemode.sh -f FILE -O BIGFILE
#-s REQUIRED SAMPLEID
#-v REQUIRED VARIANTCALLER
#-f REQUIRED FILE   					
#-o REQUIRED BIGFILE
"
}


## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -s )                    shift
                                        SAMPLEID=$1
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
		             * )			              usage
			                                  exit 1
        esac
        shift
done

awk -v ID="$SAMPLEID" '$0 !~ /^##/ {print ID "\t" "$VC" "\t" $0}' $OUTPUT >> $VC.calls.txt
