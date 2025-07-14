#!/bin/bash

#SBATCH --job-name=run_GISTIC
#SBATCH --output=run_GISTIC_j%.out
#SBATCH --error=run_GISTIC_%j.err
#SBATCH --time=1-00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu
#SBATCH --mem-per-cpu=64

# This script is used to test GISTIC2.0
# GISTIC module identifies regions of the genome that are significantly amplified or deleted across a set of samples

# Load required modules
module purge
module load singularity/3.9.6
module load gistic/2.0.23

# Set paths to input files and output directories
base_dir="/bgfs/alee/LO_LAB/Personal/Daisong/APOLLO/Tendo"
segfile="/bgfs/alee/LO_LAB/Personal/Daisong/APOLLO/Tendo/NST.seg"
refgenefile="/bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/GISTIC2.0/gistic2-2.0.23/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat"

# Set parameters
GENEGISTIC=1 
SMALLMEM=1 
BROAD=1 
BRLEN=0.5
CONF=0.90 
ARMPEEL=1 
SAVEGENE=1 
GCM=extreme

# Run GISTIC2
gistic2 -b $base_dir \
        -seg $segfile \
        -refgene $refgenefile \
        -genegistic $GENEGISTIC \
        -smallmem $SMALLMEM \
        -broad $BROAD \
        -brlen $BRLEN \
        -conf $CONF \
        -armpeel $ARMPEEL \
        -savegene $SAVEGENE \
        -gcm $GCM 

# Print when gistic is finished
echo "GISTIC2 analysis completed" 