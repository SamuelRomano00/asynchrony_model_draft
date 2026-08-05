#!/bin/bash
#SBATCH --job-name=asynchrony-sensitivity
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/sensitivity-%j.out

# Submission template for analysis/04_sensitivity_simulations.R, the only
# script in the project that needs a cluster.
#
# detect_ncores() reads SLURM_CPUS_PER_TASK, so the R code automatically uses
# exactly the number of cores requested above: change --cpus-per-task and
# nothing else. Leaving it out makes the job fall back to the node's total
# hardware core count and oversubscribe; the code warns loudly when that
# happens.
#
#   sbatch inst/hpc/sensitivity_simulations.sbatch

set -euo pipefail

module load R          # adjust to your cluster's module names
mkdir -p logs

Rscript analysis/04_sensitivity_simulations.R
