#!/bin/sh

## MR for NMR traits on CVD outcomes using different variant sets 

#SBATCH --partition=compute
#SBATCH --job-name=nmr-cvd-mr
#SBATCH --account=sc-users
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-249
#SBATCH --output=slurm_logs/%x-%A-%2a.out

## change directory
cd /PATH/TO/PROJECT/DIRECTORY

## Use Array Index to select features
echo "[LOG] Array Job ID: $SLURM_ARRAY_JOB_ID"
echo "[LOG] Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "[LOG] Job ID: $SLURM_JOB_ID"
echo "[LOG] Node ID: $SLURM_NODEID"
echo "[LOG] Node List: $SLURM_NODELIST"
echo "[LOG] Job ID: $SLURM_ARRAY_TASK_ID"

## Do some logging
echo "[LOG] Running script with ${input} as input and ${output} as output"
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "[LOG] Starting script at: $date"

## This is the container to be used
R_CONTAINER='</PATH/TO/R/CONTAINER.sif'

# This is the script that is executed
# Get with rstudioapi::getSourceEditorContext()$path
R_SCRIPT='02_nmr_cvd_mr.R'

# Enter all directories you need, simply in a comma-separated list
BIND_DIR="</PATH/OF/NECESSARY/DIRECTORIES>"

## The container 
singularity exec \
  --bind $BIND_DIR \
  $R_CONTAINER Rscript $R_SCRIPT $SLURM_ARRAY_TASK_ID

## Some more logging
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "[LOG] Finishing script at: $date"
echo "[LOG] Done!"
exit $(echo $?)
