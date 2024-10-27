#!/bin/bash
#SBATCH --job-name=medip_seq
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeDIP/logs/medip_seq_%j.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeDIP/logs/medip_seq_%j.out"

# Load the appropriate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Set conda package cache directory
export CONDA_PKGS_DIRS="/beegfs/scratch/ric.broccoli/kubacki.michal/conda_pkgs"

# Set temporary directories for various components
export TMPDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeDIP/tmp"
export SNAKEMAKE_CONDA_PREFIX="/beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs"
export XDG_CACHE_HOME="/beegfs/scratch/ric.broccoli/kubacki.michal/cache"
export SNAKEMAKE_PERSISTENT_CACHE="/beegfs/scratch/ric.broccoli/kubacki.michal/snakemake_cache"

# Create all necessary directories
mkdir -p "$TMPDIR"
mkdir -p "$SNAKEMAKE_CONDA_PREFIX"
mkdir -p "$XDG_CACHE_HOME"
mkdir -p "$SNAKEMAKE_PERSISTENT_CACHE"

# Clean up any existing temporary files before starting
rm -rf "$TMPDIR"/*
rm -rf ~/.cache/snakemake/*

# Clear any potential locks
snakemake --unlock

# Remove existing conda environments to force recreation
rm -rf "$SNAKEMAKE_CONDA_PREFIX"/*

# Run a dry-run first to check workflow
echo "Performing dry-run..."
snakemake --dry-run -p -n

# Run Snakemake with optimal settings for MeDIP-seq pipeline
snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --use-conda \
    --conda-prefix "/beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs" \
    --cores 64 \
    --keep-going \
    --rerun-incomplete \
    --latency-wait 60 \
    --printshellcmds \
    --resources mem_mb=128000 \
    --nolock \
    2>&1 | tee "logs/snakemake_$(date +%Y%m%d_%H%M%S).log"

# Create summary of run
snakemake --summary > "logs/workflow_summary_$(date +%Y%m%d_%H%M%S).txt"

# Add cleanup at the end of your script
trap 'rm -rf "$TMPDIR"/*' EXIT
