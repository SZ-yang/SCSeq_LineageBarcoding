#!/bin/bash

#SBATCH --job-name=scContrastiveLearn
#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=18:00:00
#SBATCH --mem=60G
#SBATCH --array=1-9%3  # Runs 9 jobs in total, at most 3 simultaneously

# Define parameters based on SLURM_ARRAY_TASK_ID
# batch_sizes=(25 50 75)
batch_sizes=(25 25 25)
# temperatures=(0.5 0.5 0.5)
temperatures=(0.1 0.5 1.0)
index=$((SLURM_ARRAY_TASK_ID - 1))  # Adjust index for zero-based arrays
batch_size=${batch_sizes[$((index / 3))]}
temperature=${temperatures[$((index % 3))]}

# Define directory based on parameters and ensure it exists
formatted_date=$(date "+%m%d")
formatted_date=$((10#$formatted_date))
# echo $formatted_date

result_dir="/home/users/syang71/kzlinlab/projects/scContrastiveLearn/out/joshua/data/feat_${formatted_date}/feat_${formatted_date}__bs${batch_size}_temp${temperature}"
mkdir -p $result_dir

err_out = "/home/users/syang71/kzlinlab/projects/scContrastiveLearn/git/scContrastiveLearn_Joshua/slurm/err_out"
# Specify output and error files within the newly created directory
#SBATCH --output=$err_out/output_%A_%a.out
#SBATCH --error=$err_out/error_%A_%a.err

# Load environment
source ~/.bashrc
source activate linBarCdin

# Run the Python script with parameters and directory path for output

# Run your Python script
script_path="/home/users/syang71/kzlinlab/projects/scContrastiveLearn/git/scContrastiveLearn_Joshua/main/scContrastiveLearning_train_test.py"
python $script_path --batch_size $batch_size --temperature $temperature --output_dir $result_dir