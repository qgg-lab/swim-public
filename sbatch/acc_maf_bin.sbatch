#!/bin/bash

# Job name:
#SBATCH --job-name=acc_mafbin
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Memory per node:
#SBATCH --mem=6G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=00:20:00

# ============================================================

module load R

Rscript /mnt/ls15/scratch/users/dingrong/swim/pharp_acc/acc_maf_bin.R $inputname $outputname

echo $(date) done.main.proces
