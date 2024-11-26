#!/bin/bash
#SBATCH --job-name=medip_seq
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=48:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeDIP/logs/medip_seq_%j.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeDIP/logs/medip_seq_%j.out"

# Set up environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Create directories
# mkdir -p logs/cluster
# mkdir -p output/{trimmed,aligned,dedup,peaks,fastqc,qc,func_annot,diff_meth,bigwig,meth_level}

# # Clean up and unlock
# rm -f output/trimmed/*.fastq.gz
# rm -f output/aligned/*.bam
# rm -f .snakemake/metadata/*failed*
snakemake --unlock

# Create submission script
cat > submit_job.sh <<'EOF'
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
EOF

chmod +x submit_job.sh

# Run snakemake
snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --cluster "sbatch --parsable" \
    --jobs 64 \
    --latency-wait 120 \
    --restart-times 3 \
    --keep-going \
    --rerun-incomplete \
    --default-resources \
        mem_mb=8000 \
        runtime=240 \
    2>&1 | tee "logs/snakemake_$(date +%Y%m%d_%H%M%S).log"

# Cleanup
# rm -f submit_job.sh
