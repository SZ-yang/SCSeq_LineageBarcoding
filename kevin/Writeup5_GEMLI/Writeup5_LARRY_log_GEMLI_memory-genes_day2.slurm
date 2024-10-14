#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-20c96g
#SBATCH --time=48:00:00

#SBATCH --job-name=gemli_memory
#SBATCH --mem-per-cpu=30gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/slurm/slurm_%j.out

R CMD BATCH --no-save --no-restore Writeup5_LARRY_log_GEMLI_memory-genes_day2.R ../../../../out/kevin/Writeup5/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.Rout


