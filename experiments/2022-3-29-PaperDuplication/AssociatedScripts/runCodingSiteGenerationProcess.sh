#!/bin/bash
# The interpreter used to execute the script
#SBATCH --job-name=broadTimecourseDuplicationDataAnalysis
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=clhaynes@umich.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2g
#SBATCH --time=00-08:00:00
#SBATCH --account=zamanlh0
#SBATCH --array=16-50

USERNAME=clhaynes
EXPERIMENT_ID=2022-3-29-PaperDuplication

EXPERIMENT_DIR=/scratch/zamanlh_root/zamanlh0/${USERNAME}/${EXPERIMENT_ID}

echo "Experiment: ${EXPERIMENT_ID}"

ANALYSIS_TIME=$((1000 * SLURM_ARRAY_TASK_ID))
echo $ANALYSIS_TIME

cp CodingSiteGeneratorHPCCCopy.py ${EXPERIMENT_DIR}/CodingSiteGeneratorHPCCCopy_${ANALYSIS_TIME}.py

cd ${EXPERIMENT_DIR}

module load gcc/11.2.0
python3 CodingSiteGeneratorHPCCCopy_${ANALYSIS_TIME}.py ${ANALYSIS_TIME}

rm CodingSiteGeneratorHPCCCopy_${ANALYSIS_TIME}.py