#!/bin/bash
#SBATCH --job-name=MICROPEP
#SBATCH --output=/storage/scratch01/groups/md/microproteins/Microproteins_pipeline/Micropeptidome/snakemake_controller_%j.out
#SBATCH --error=/storage/scratch01/groups/md/microproteins/Microproteins_pipeline/Micropeptidome/snakemake_controller_%j.err
#SBATCH --time=22:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail

source /home/sbarber/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

snakemake --use-conda --slurm -j 32 --conda-frontend conda --conda-prefix /home/sbarber/conda_envs \
  --rerun-incomplete \
  --latency-wait 60

