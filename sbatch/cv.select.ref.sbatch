#!/bin/bash

# Job name:
#SBATCH --job-name=ref_select
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
#SBATCH --time=2:00:00

# ============================================================

/mnt/research/qgg/software/vcftools-v0.1.16/bin/vcftools \
--gzvcf $raw_gzvcf \
--remove $cv_split \
--max-missing 0.1 \
--mac 2 \
--recode \
--out $refname 

/mnt/research/qgg/software/bcftools-1.9/bcftools view $refname.recode.vcf -Oz -o $refname.vcf.gz
 
/mnt/research/qgg/software/bcftools-1.9/bcftools index -f -t --threads $SLURM_CPUS_ON_NODE $refname.vcf.gz

echo $(date) done.main.process
