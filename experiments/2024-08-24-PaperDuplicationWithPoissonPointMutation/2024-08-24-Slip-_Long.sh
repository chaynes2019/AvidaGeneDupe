#!/bin/bash
# The interpreter used to execute the script

#SBATCH --job-name=AvidaGeneDupeFullRep_LongAncestorGenomeControl
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=clhaynes@umich.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4g
#SBATCH --time=00-04:00:00
#SBATCH --account=zamanlh0
#SBATCH --array=1-30

# -- I like to define helpful variables up top --
module load gcc/11.2.0

USERNAME=clhaynes
EXPERIMENT_ID=2024-08-24-PaperDuplicationWithPoissonPointMutation

OUTPUT_DIR=/scratch/zamanlh_root/zamanlh0/${USERNAME}/${EXPERIMENT_ID}/Slip-_Long
CONFIG_DIR=/home/${USERNAME}/Documents/AvidaGeneDupe/experiments/${EXPERIMENT_ID}/hpcc/config

SEED_OFFSET=1770

SEED=$((SEED_OFFSET + SLURM_ARRAY_TASK_ID - 1))

JOB_ID=${SLURM_ARRAY_TASK_ID}

RUN_DIR=${OUTPUT_DIR}/run_${SEED}

mkdir -p ${RUN_DIR}

cd ${RUN_DIR}

#Don't use the asterisk: actually write everything out so you know what you're working with!

cp ${CONFIG_DIR}/avida .
cp ${CONFIG_DIR}/avida.cfg .
cp ${CONFIG_DIR}/longAncestralOrganism.org .
cp ${CONFIG_DIR}/environment.cfg .
cp ${CONFIG_DIR}/eventsLongAncestralOrganism.cfg .
cp ${CONFIG_DIR}/instset-heads___sensors_NONE.cfg .
cp ${CONFIG_DIR}/analyze.cfg .

EXECUTE="avida -s ${SEED} -set COPY_MUT_PROB 0.0 -set DIVIDE_POISSON_MUT_MEAN 0.25 -set COPY_INS_PROB 0.0 -set COPY_DEL_PROB 0.0 -set DIVIDE_INS_PROB 0.00 -set DIVIDE_DEL_PROB 0.00 -set DIVIDE_SLIP_PROB 0.0 -set SLIP_FILL_MODE 0 -set EVENT_FILE eventsLongAncestralOrganism.cfg -set STERILIZE_UNSTABLE 1"
echo ${EXECUTE} > cmd.log
./${EXECUTE} > run.log
./${EXECUTE} -a > analyze.log

rm avida
rm avida.cfg
rm longAncestralOrganism.org
rm environment.cfg
rm eventsLongAncestralOrganism.cfg
rm instset-heads___sensors_NONE.cfg
rm analyze.cfg