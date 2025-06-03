#!/bin/bash
#SBATCH --job-name=scplus
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=300G
#SBATCH --time=72:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH -p musc3

source ~/.bashrc
conda activate scenicplus

cd /scratch/bryangranger/Haley/scenicplus/scplus_pipeline/Snakemake

snakemake --cores 20
