#!/bin/bash

# Job name:
#SBATCH --job-name=admixture
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=24
#
# Memory per node:
#SBATCH --mem=32G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=23:50:58

# ============================================================
cd /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/

/mnt/home/dingrong/software/dist/admixture_linux-1.3.0/admixture --cv select.all.pruned.bed $K -j24 | tee log${K}.out
