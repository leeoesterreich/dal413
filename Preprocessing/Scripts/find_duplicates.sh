#!/bin/bash
# This is a bash script used to run python script to find duplicates file in the AWS bucket.

#SBATCH --job-name=Duplicates_files
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -t 1-00:00
#SBATCH --output=duplicates_output_%j.txt  
#SBATCH --error=duplicates_error_%j.txt 

# Load necessary modules
module purge
module load gcc/12.2.0 
module load python/anaconda2.7-4.4.0_genomics

# Unset conflicting environment variables
unset PYTHONHOME
unset PYTHONPATH

# Ensure the correct Python interpreter is used
export PATH=/ihome/alee/rak373/.conda/envs/s3cmd-env/bin:$PATH

# Activate the Conda environment
source activate /ihome/alee/rak373/.conda/envs/s3cmd-env

# Correct path to the Python script
python_script="/bgfs/alee/LO_LAB/Personal/Rahul/Com_housekeeping/find_duplicates_by_md5.py"

# Grant execute permission to the Python script (if necessary)
chmod +x "$python_script"

# Run the Python script
python "$python_script"

# Deactivate the Conda environment (optional)
conda deactivate