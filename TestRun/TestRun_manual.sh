usage()
{
echo "Executing TestRun_manual.sh
#
# Usage $0 -d /path/to/MosaiC-All -s $git/TestRun/SampleID_Test -c $git/config/Mosaic-All.TestRun.config -o $git/TestRun/OUTPUT
# 
# Options
# -d git                REQUIRED: path to where the git is cloned into
# -s sampleID_list      REQUIRED: ID list
# -c config_file        REQUIRED: config_file that can be found in $git/config - make sure to amend accordingly
# -o output             REQUIRED: Output directory within TestRun
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
            -o )                    shift
                                    outdir=$1
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
resources="$testrun_dir/Resources"

if [ ! -d "${outdir}" ]; then
    mkdir -p ${outdir}
fi

if [ ! -d "${resources}" ]; then
    mkdir -p ${resources}
fi

source $config_file

# check and download resources:

# 1. Reference genome

#if [ ! -e "$REFGEN" ]; then
#    # Download the file if it doesn't exist
#    wget "https://storage.googleapis.com/genomics-public-data/references/hs37d5/hs37d5.fa.gz" -P "$resources"
#    # Extract the tar.gz file
#    gunzip "$resources/hs37d5.fa.gz"
#    samtools faidx "$resources/hs37d5.fa" -o "$resources/hs37d5.fa.fai"
#fi

# 2. dbSNP

if [ ! -e "$DBSNP" ]; then
    # download
    wget "https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz" -P "$resources"
    # Extract the tar.gz file
    gunzip "$resources/dbsnp_138.b37.vcf.gz"
fi

# 3. Panel of Normals for Mutect2

if [ ! -e "$PON_A" ]; then
    # panel of normals
    wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf -P $resources
    wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf.idx -P $resources
    # germline resources
    wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf -P $resources
    wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx -P $resources
fi

# 4. Germline Resources for Mutect2

if [ ! -e "$GERMLINE_RESOURCES" ]; then
    # panel of normals
    wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf -P $resources
    wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf.idx -P $resources
    # germline resources
    wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf -P $resources
    wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx -P $resources
fi


#execute the MasterScript
$git_dir/MASTERSCRIPT_MosaiC-All.sh -s $sampleID_list -o $outdir -c $config_file
