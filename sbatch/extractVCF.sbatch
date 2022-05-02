#!/bin/bash

# Job name:
#SBATCH --job-name=extractVCF
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=2
#
# Memory per node:
#SBATCH --mem=32G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=8:00:00

# ============================================================

source "$env"

zcat $file1 | perl $EXTRACTVCF - $varList | gzip > $output
$BCFTOOLS index --tbi $output

echo $(date) done.main.process
