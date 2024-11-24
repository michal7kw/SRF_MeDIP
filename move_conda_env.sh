bisulfite_seq            /home/kubacki.michal/.conda/envs/bisulfite_seq
jupyter_nb            *  /home/kubacki.michal/.conda/envs/jupyter_nb
r_env                    /home/kubacki.michal/.conda/envs/r_env
snakemake                /home/kubacki.michal/.conda/envs/snakemake



# Create the target directory if it doesn't exist
mkdir -p /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/

# Move each environment
for env in bisulfite_seq jupyter_nb r_env snakemake; do
    echo "Moving $env environment..."
    mv "/home/kubacki.michal/.conda/envs/$env" "/beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/"
done

# Create symbolic links to maintain functionality
for env in bisulfite_seq jupyter_nb r_env snakemake; do
    echo "Creating symlink for $env..."
    ln -s "/beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/$env" "/home/kubacki.michal/.conda/envs/$env"
done