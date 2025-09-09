#!/bin/bash
# Script to run Transformer training

# General SLURM configurations
#SBATCH --job-name=TRANSFORMER_GPU_TRAIN
#SBATCH --output=transformer_gpu_train_%j.out
#SBATCH --error=transformer_gpu_train_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu 
#SBATCH --mem-per-cpu=4000        

# GPU-specific configurations
#SBATCH --nodes=1
#SBATCH --time=04:00:00             # User-specified time
#SBATCH --ntasks-per-node=1         # Based on 1 GPUs per node
#SBATCH --gres=gpu:1                # User-specified 1 GPUs per node
#SBATCH --cluster=gpu
#SBATCH --partition=a100

# Path to the persistent conda environment on IX
ENV_PREFIX="/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/envs/metdx_torch2_clean"

# Load an Anaconda-based Python module and activate the env (CRC guidelines)
echo "Loading Anaconda Python and activating env: $ENV_PREFIX ..."
module purge
module load cuda/11.7.1
module load python/ondemand-jupyter-python3.10 || module load python/anaconda3.9-2021.11
source activate "$ENV_PREFIX"

# Prevent leaking user site-packages that cause NumPy/pyarrow conflicts
export PYTHONNOUSERSITE=1

# Ensure required packages exist in the env (idempotent)
python -m pip install -q --upgrade pip >/dev/null 2>&1
python -m pip install -q "numpy==1.26.4" "pandas==2.1.4" "scikit-learn==1.3.2" "joblib==1.3.2" "matplotlib==3.8.4" "wandb==0.18.7"

# Sanity prints
echo "Python: $(which python)"
python - <<'PY'
import sys, torch
print('Python version:', sys.version)
print('Torch:', torch.__version__)
print('CUDA available:', torch.cuda.is_available())
print('Device count:', torch.cuda.device_count())
PY

# Run training
echo "Running Transformer training script on GPU..."
python "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/METDX/01_model_training_transformer.py"

echo "Script execution finished."
