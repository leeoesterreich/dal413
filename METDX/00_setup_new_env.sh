#!/bin/bash
# Script to create a new, clean conda environment for the Transformer training.
set -e # Exit immediately if a command exits with a non-zero status.

# Path to the new persistent conda environment on IX
export NEW_ENV_PREFIX="/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/envs/metdx_torch2_clean"

# Remove the old environment directory if it exists
if [ -d "$NEW_ENV_PREFIX" ]; then
    echo "Removing existing environment at $NEW_ENV_PREFIX..."
    rm -rf "$NEW_ENV_PREFIX"
fi

# Load a compatible Python module
echo "Loading Python module..."
module purge
module load python/ondemand-jupyter-python3.10 || module load python/anaconda3.9-2021.11

# Prevent leaking user site-packages that cause NumPy/pyarrow conflicts
export PYTHONNOUSERSITE=1

# Create the new environment
echo "Creating new conda environment at $NEW_ENV_PREFIX..."
conda create --prefix "$NEW_ENV_PREFIX" -y python=3.10

# Activate the new environment
echo "Activating new environment..."
source activate "$NEW_ENV_PREFIX"

# Install PyTorch with a specific CUDA version from the pytorch channel
echo "Installing PyTorch... This may take several minutes."
conda install pytorch pytorch-cuda=11.7 mkl -c pytorch -c nvidia -y

# Install other required packages
echo "Installing other packages..."
python -m pip install -q "numpy==1.26.4" "pandas==2.1.4" "scikit-learn==1.3.2" "joblib==1.3.2" "matplotlib==3.8.4" "wandb"

echo "Environment setup complete."
