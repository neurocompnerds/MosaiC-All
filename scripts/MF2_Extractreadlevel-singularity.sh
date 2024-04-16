#!/bin/sh -l
#SBATCH -J MF_Extractreadlevel-singularity.sh
#SBATCH -o /hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/Mosaic/MosaiC-All/TestRun/MF.extract.slurm-%j.out

#SBATCH -p skylake,a100,icelake
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=07:00:00
#SBATCH --mem=40GB
#SBATCH --gres=tmpfs:150G
# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=%u@adelaide.edu.au

usage()
{
echo "
#
#README
#This script is designed for Mosaic-All wrapper script

#
#SCRIPT EXECUTION
#sbatch $0 -b BAMDIR -s SampleID -c CONFIG_FILE -o OUTDIR  
#all of this will be "specified" in MasterScript
#
#-s REQUIRED SampleID                       (e.g. A unique identifier for your sample that is also in the BAM file name)
#-b REQUIRED Directory of BAM files         (e.g. /path/to/bam_file_directory/)
#-c REQUIRED The script configuration file  (/path/to/config_file.cfg)
#-o REQUIRED Directory of outputs           (e.g. /path/to/output/)
#
"
}


## Set Variables ##

while [ "$1" != "" ]; do
        case $1 in
                -s )                    shift
                                        SampleID=$1
                                        ;;
                -b )                    shift
                                        BAMDIR=$1
                                        ;;
                -o )                   	shift
                                        OUTDIR=$1
					;;
		-c )                   	shift
                                        CONFIG_FILE=$1
					;;
                 * )                    usage
                                        exit 1
        esac
        shift
done

## Check for all required arguments
for ARG in $BAMDIR $OUTDIR $CONFIG_FILE $SampleID; do
    if [ -z "${ARG}" ]; then
        usage
        echo "## ERROR: Please check that you have supplied all required argunments for this script"
        exit 1
    fi
done

# Load modules 
module purge
module use /apps/skl/modules/all
module load Singularity/3.7.1

source $CONFIG_FILE

BAMFILE=$(find "$BAMDIR" -type f -name "$SampleID*.bam")
BAMINDEX=$(find "$BAMDIR" -type f -name "$SampleID*.bai")
BAMprefix=$(basename "$BAMFILE" | sed 's/\.[^.]*$//')
REFGEN_DIR=$(dirname $REFGEN)
REFSEQ=$(basename $REFGEN)

/usr/bin/mkdir -p ${TMPDIR}
/usr/bin/cp -r ${BAMFILE} ${TMPDIR}
/usr/bin/cp -r ${BAMINDEX} ${TMPDIR}
/usr/bin/cp -r $OUTDIR/$SampleID.phasingInput.bed ${TMPDIR}

ls ${TMPDIR}


#execute the read-level extractions
singularity run -B ${REFGEN_DIR}:/RefSeq -B $MFORECAST:/MForecastDir -B ${TMPDIR}:/tmp $MFORECAST/mosaicforecast_0.0.1.sif ReadLevel_Features_extraction.py /tmp/$SampleID.phasingInput.bed /tmp/$SampleID.features.bed /tmp /RefSeq/$REFSEQ /MForecastDir/$CONFIG/k24.umap.wg.bw ${SLURM_NTASKS} bam

/usr/bin/cp -r ${TMPDIR}/$SampleID.features.bed ${OUTDIR}/$SampleID.features.bed

echo "## INFO: (MF2_Extractreadlevel-singularity.sh) SUCCESS for $SampleID" >> ${OUTDIR}/$SampleID.pipeline.log
