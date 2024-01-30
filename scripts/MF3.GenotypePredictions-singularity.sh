#!/bin/sh -l
#SBATCH -J MF_genotypePredictions-singularity.sh
#SBATCH -o /home/%u/Mosaic-All/Log/genotype.GA.slurm-%j.out

#SBATCH -p skylake,icelake
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=00:10:00
#SBATCH --mem=10GB
#SBATCH --gres=tmpfs:40G
#Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=a1742674@adelaide.edu.au

usage()
{
echo "
#
#README
#This script is designed for Mosaic-All wrapper script
#This process is to predict the snv to be either heterozygous, mosaic or germline etc. using the reads features extrcated in the previous step (MF-Extractreadlevel).
#
#
#SCRIPT EXECUTION
#sbatch $0 -s SampleID -c CONFIG_FILE -o OUTDIR  
#
#-s REQUIRED SampleID                       (e.g. A unique identifier for your sample that is also in the BAM file name)
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
                -o )                    shift
                                        OUTDIR=$1
                                        ;;
                -c )                    shift
                                        CONFIG_FILE=$1
                                        ;;
                * )                     usage
                                        exit 1
        esac
        shift
done

## Check for all required arguments
for ARG in $OUTDIR $CONFIG_FILE $SampleID; do
    if [ -z "${ARG}" ]; then
        usage
        echo "## ERROR: Please check that you have supplied all required argunments for this script"
        exit 1
    fi
done

# Load modules
module purge
module load Singularity/3.10.5

source $CONFIG_FILE

singularity run -B $OUTDIR:/outDir -B $MFORECAST:/MForecastDir $MFORECAST/mosaicforecast_0.0.1.sif Prediction.R /outDir/$SampleID.features.bed /MForecastDir/models_trained/50xRFmodel_addRMSK_Refine.rds Refine /outDir/$SampleID.mosaicforecast.genotype.predictions.refined.bed
