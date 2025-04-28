#! /bin/bash
#SBATCH --job-name=network_impact_analysis
#SBATCH --output=network_impact_analysis.out
#SBATCH --error=network_impact_analysis.err
#SBATCH --time=01:00:00
#SBATCH --partition=htc
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB


source /bgfs/alee/LO_LAB/Personal/Daisong/miniconda3/bin/activate reactome_env

python network_impact_analysis.py \
  --genes $(cat /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Superfreq/somatic_mutations_human.txt | tr '\n' ',') \
  --neo4j-uri "bolt://localhost:7687" \
  --neo4j-user "neo4j" \
  --neo4j-password "your_password" \
  --output "pathway_impact.pdf"