#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-20c96g
#SBATCH --time=48:00:00

#SBATCH --job-name=vp
#SBATCH --mem-per-cpu=90gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup7/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup7/slurm/slurm_%j.out

R CMD BATCH --no-save --no-restore Writeup7_variancePartition.R ../../../../out/kevin/Writeup7/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.txt
