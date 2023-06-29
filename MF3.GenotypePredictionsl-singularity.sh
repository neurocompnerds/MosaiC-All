#!/bin/sh -l
#SBATCH -J MF_genotypePredictions-singularity.sh
#SBATCH -o /gpfs/users/a1149120/MosaicForecast/genotype.GA.slurm-%j.out
#SBATCH -A robinson
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
#For standalone script,look for the original
#This process is to predict the snv to be either heter, mosaic or germline etc using the reads features extrcated in the previous step (MF-Extractreadlevel)
#
#
#SCRIPT EXECUTION
#SCRIPTDIR=/hpcfs/groups/phoenix-hpc-neurogenetics/Mosaic-All/Mosaic-S
#sbatch $SCRIPTDIR/MF3.GenotypePredictionsl-singularity.sh -s SAMPLE -c CONFIG_FILE -o OUTDIR  
#all of this will be "specified" in MasterScript
#
#-s REQUIRED SampleID   					(e.g 004P, or 004M, or 004F no suffix)
#-c REQUIRED CONFIG_FILE   		(/gpfs/users/a1149120/MosaicHunter_single)
#-o REQUIRED OUTDIR     	        (e.g hs37d5.fa) 
"
}


## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -s )                    shift
                                        SAMPLE=$1
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


export HOME=/hpcfs/users/$USER

#define directories
module purge
module load Singularity/3.10.5

source $CONFIG_FILE


#singularity run -B /hpcfs -B ${TMPDIR}/tmp /hpcfs/groups/phoenix-hpc-neurogenetics/executables/MosaicForecast-master/mosaicforecast_0.0.1.sif Prediction.R /tmp/${sample[$SLURM_ARRAY_TASK_ID]}.features.bed $MFORECAST/models_trained/50xRFmodel_addRMSK_Refine.rds Refine $DIR/${sample[$SLURM_ARRAY_TASK_ID]}.genotype.predictions.refined.bed
#singularity run -B /hpcfs -B ${TMPDIR} /hpcfs/users/$USER/mosaicforecast_0.0.1.sif Prediction.R $DIR/${sample[$SLURM_ARRAY_TASK_ID]}.features.bed $MFORECAST/models_trained/50xRFmodel_addRMSK_Refine.rds Phased $DIR/${sample[$SLURM_ARRAY_TASK_ID]}.genotype.predictions.phased.bed

singularity run -B /hpcfs $MFORECAST/mosaicforecast_0.0.1.sif Prediction.R $OUTDIR/$SAMPLE.features.bed $MFORECAST/models_trained/50xRFmodel_addRMSK_Refine.rds Refine $OUTDIR/$SAMPLE.genotype.predictions.refined.bed

