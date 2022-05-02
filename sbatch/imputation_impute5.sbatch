#!/bin/bash

# Job name:
#SBATCH --job-name=imputation_impute5
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
#SBATCH --time=4:00:00

# ============================================================

source "$env"

$IMPUTE5 \
--h "$ref_panel".vcf.gz \
--g $input.vcf.gz \
--r $chr \
--b 300 \
--ne 100 \
--o $output.vcf \
--out-gp-field --out-ap-field \
--threads $SLURM_CPUS_ON_NODE

$BCFTOOLS view $output.vcf --threads $SLURM_CPUS_ON_NODE -Oz -o $output.vcf.gz 

$BCFTOOLS index -f -t --threads $SLURM_CPUS_ON_NODE $output.vcf.gz

echo $(date) done.main.process
   
   
