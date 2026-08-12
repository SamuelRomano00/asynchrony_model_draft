#!/bin/bash
#SBATCH --job-name=asynchrony-sensitivity
#SBATCH --partition=scicore
#SBATCH --qos=1day
#SBATCH --cpus-per-task=128
#SBATCH --mem=64G
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

module purge
module load R/4.4.1-foss-2023b       # adjust to your cluster's module names

mkdir -p logs

# One BLAS thread per process: mclapply already provides the parallelism, and
# letting each of the 128 forks spawn its own BLAS threads oversubscribes the
# node by a factor of five and collapses throughput.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export FLEXIBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript analysis/04_sensitivity_simulations.R