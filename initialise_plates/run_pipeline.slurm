#!/bin/bash 

# === SLURM SETTINGS ===

#SBATCH --job-name=DGE-init
#SBATCH --time=96:00:00
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

# =======================


# === USER INPUT ========

# Set base paths - This is the filepath to your analysis directory:
ROOT=/mnt/odap-beegfs/a008/06_DGE_cell/initialise_tissues
BASE=/mnt/odap-beegfs/a008

# =======================


# === FILE PATHS ========
CONT=/mnt/odap-beegfs/software/singularity-images

CONT_DGE=${CONT}/dge-pipeline-container-v2.sif
SNAKEFILE=${ROOT}/Snakefile
CONFIG=${ROOT}/config.yaml

# =======================

# === RUN PIPELINE ======
singularity exec --bind ${BASE},${CONT} ${CONT_DGE} \
conda run -n DGE_pipeline snakemake -s ${SNAKEFILE} --configfile ${CONFIG} --cores 64 -p --use-conda --use-singularity --rerun-incomplete

# =======================

