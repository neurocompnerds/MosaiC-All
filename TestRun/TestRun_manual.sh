usage()
{
echo "Executing TestRun_manual.sh
#
# Usage $0 -d /path/to/MosaiC-All -s SampleID_Test -c $CONFIG_file - [ - h | --help ]
#  ..........
# Options
# -d git                REQUIRED: path to where the git is cloned into
# -s sampleID_list      REQUIRED: ID list
# -c config_file        REQUIRED: config_file that can be found in $git/config - make sure to amend accordingly
# -o output             REQUIRED
# -h or --help  Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
#
# See: https://github.com/neurocompnerds/MosaiC-All for history and new updates.
#
"
}

## Set Variables ##
while [ "$1" != "" ]; do
    case $1 in
            -d )                    shift
                                    git_dir=$1
                                    ;;
            -s )                    shift
                                    sampleID_list=$1
                                    ;;
            -c )                    shift
                                    config_file=$1
                                    ;;
            -h | --help )           usage
                                    exit 0
                                    ;;
            *  )                    usage
                                    exit 1
    esac
    shift
done

# Define variables
testrun_dir="$git_dir/TestRun"
bam_dir="$testrun_dir/BAM"
outdir="$testrun_dir/OUTPUT"
resources="$testrun_dir/Resources"

if [ ! -d "${outdir}" ]; then
    mkdir -p ${outdir}
fi

if [ ! -d "${resources}" ]; then
    mkdir -p ${resources}
fi

GIT=./MosaiC-All #specify where MosaiC-All is cloned into
CONFIG=$GIT/config/Mosaic-All.TestRun.config
$GIT/MASTERSCRIPT_MosaiC-All.sh -s $GIT/TestRun/SampleID_Test -o $GIT/TestRun/Output -c $CONFIG
