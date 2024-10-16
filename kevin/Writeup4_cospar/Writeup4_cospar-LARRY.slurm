#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=48:00:00

#SBATCH --job-name=cospar
#SBATCH --mem-per-cpu=100gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup4/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup4/slurm/slurm_%j.out

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cospar

conda env list
python --version
conda list

python Writeup4_cospar-LARRY.py > ../../../../out/kevin/Writeup4/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.txt