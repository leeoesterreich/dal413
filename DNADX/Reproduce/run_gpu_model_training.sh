#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=a100
#SBATCH --job-name=gpu_ml_tuning
#SBATCH --output=gpu_ml_tuning_%j.out
#SBATCH --error=gpu_ml_tuning_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu # Copied from your CPU script
#SBATCH --mem-per-cpu=2000         # Copied from your CPU script, adjust if necessary for GPU nodes

# --- Environment Setup ---
echo "Job started on $(hostname) at $(date)"
echo "SLURM Job ID: ${SLURM_JOB_ID}"
echo "SLURM GPU IDs allocated: ${SLURM_JOB_GPUS:-none}"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-none}"

# Load Python module (adjust if a different version or environment is needed for GPU nodes)
echo "Loading Python module..."
module load python/3.7.0
echo "Python module loaded."

# --- Script Execution ---
PYTHON_SCRIPT_PATH="/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/scaling/model_training.py"

echo "--------------------------------------------------------------------"
echo "INFO: This job has allocated a GPU."
echo "INFO: However, the Python script: ${PYTHON_SCRIPT_PATH}"
echo "INFO: uses scikit-learn for ElasticNet, which primarily runs on the CPU."
echo "INFO: The GPU will likely NOT be used by scikit-learn for model training."
echo "INFO: To leverage GPU acceleration for model training, the Python code"
echo "INFO: would need to be modified to use libraries like RAPIDS cuML,"
echo "INFO: PyTorch, or TensorFlow."
echo "--------------------------------------------------------------------"

echo "Running Python script: ${PYTHON_SCRIPT_PATH}"
# It's good practice to navigate to the script's directory or ensure paths within the script are absolute/relative to a known point
# Assuming the script handles its paths correctly or is run from a directory where it can find its data.
# For now, directly calling the script with its full path.
python "${PYTHON_SCRIPT_PATH}"

echo "Python script execution finished at $(date)."
echo "Job finished." 