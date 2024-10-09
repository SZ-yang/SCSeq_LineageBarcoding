#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-20c96g
#SBATCH --time=48:00:00

#SBATCH --job-name=scvi_rosmap
#SBATCH --mem-per-cpu=80gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/subject-de/out/kevin/Writeup11/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/subject-de/out/kevin/Writeup11/slurm/slurm_%j.out

source ~/miniconda3/etc/profile.d/conda.sh
conda activate scVI

conda env list
python --version
conda list

python Writeup11_rosmap_scvi.py > ../../../../../out/kevin/Writeup11/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.txt


