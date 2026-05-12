#!/bin/bash
#SBATCH --job-name=BNKs
#SBATCH --array=1-100
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --output=BNKs_logs/%x_%A_%a.out
#SBATCH --error=BNKs_logs/%x_%A_%a.err

mkdir -p BNKs_logs

!CSV="ids.csv"

# If CSV has a header, add +1 to skip it
LINE=$((SLURM_ARRAY_TASK_ID + 1))

echo "SLURM array task: $SLURM_ARRAY_TASK_ID"
echo "Now processing: $sub_ses"

bash script.sh "$sub_ses"

echo "Finished processing: $sub_ses"