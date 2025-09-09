#!/bin/bash

# Load the same environment as the sbatch script
module purge
module load gcc/8.2.0 r/3.6.0
module load python/ondemand-jupyter-python3.10
source activate /ix1/alee/LO_LAB/Personal/Daisong/griffin_env_new
export PATH=/ix1/alee/LO_LAB/Personal/Daisong/griffin_env_new/bin:$PATH

# Change to the correct directory
cd /ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/Griffin/griffin_GC_and_mappability_correction

echo "=== UNLOCKING DIRECTORY ==="
snakemake -s griffin_GC_and_mappability_correction.snakefile --unlock --cores 1

echo " Directory unlocked!"

echo "=== RESUBMITTING BATCH 4 ==="
sbatch griffin_GC_and_mappability_correction.sbatch
