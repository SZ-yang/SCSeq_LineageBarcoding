#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=48:00:00

#SBATCH --job-name=scvi
#SBATCH --mem-per-cpu=100gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup2/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup2/slurm/slurm_%j.out

source ~/miniconda3/etc/profile.d/conda.sh
conda activate scVI

conda env list
python --version
conda list

python Writeup2_LARRY_scvi.py > ../../../../out/kevin/Writeup2/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.txt