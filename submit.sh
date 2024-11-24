#!/bin/bash
#SBATCH {cluster.account}
#SBATCH --nodes={cluster.nodes}
#SBATCH --ntasks={cluster.ntasks}
#SBATCH --mem={cluster.mem}
#SBATCH --time={cluster.time}
#SBATCH --output={cluster.output}
#SBATCH --error={cluster.error}

{exec_job}
