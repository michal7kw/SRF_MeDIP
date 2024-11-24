#!/bin/bash
#SBATCH --account=kubacki.michal
#SBATCH --partition=normal
#SBATCH --ntasks={threads}
#SBATCH --mem={resources.mem_mb}M
#SBATCH --time=$(printf '%02d:%02d:00' $((${resources.runtime}/60)) $((${resources.runtime}%60)))
#SBATCH --job-name=snakemake_{rule}
#SBATCH --output=logs/cluster/{rule}_%j.out
#SBATCH --error=logs/cluster/{rule}_%j.err

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

{cmd}
