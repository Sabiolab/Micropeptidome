#!/bin/bash
#SBATCH --job-name=MICROPEP
#SBATCH --output=MICROPEP_%j.out
#SBATCH --error=MICROPEP_%j.err
#SBATCH --time=22:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G

set -euo pipefail

# change the conda path!
source /home/sbarber/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

snakemake --use-conda --slurm -j 32 --conda-frontend conda --conda-prefix /home/sbarber/conda_envs \
  --rerun-incomplete \
  --latency-wait 60

