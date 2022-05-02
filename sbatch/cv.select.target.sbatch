#!/bin/bash

# Job name:
#SBATCH --job-name=target_select
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

if [[ $chip == 1 ]]; then

/mnt/research/qgg/software/vcftools-v0.1.16/bin/vcftools \
--gzvcf $raw_gzvcf \
--keep $cv_split \
--snps $snp_list \
--max-missing 0.1 \
--mac 2 \
--recode \
--out $targetname 

/mnt/research/qgg/software/bcftools-1.9/bcftools view $targetname.recode.vcf -Oz -o $targetname.vcf.gz

/mnt/research/qgg/software/bcftools-1.9/bcftools index -f -t --threads $SLURM_CPUS_ON_NODE $targetname.vcf.gz

else

/mnt/research/qgg/software/vcftools-v0.1.16/bin/vcftools \
--gzvcf $raw_gzvcf \
--keep $cv_split \
--recode \
--out $targetname.raw 

/mnt/research/qgg/software/bcftools-1.9/bcftools view $targetname.raw.recode.vcf -Oz -o $targetname.raw.vcf.gz 

/mnt/research/qgg/software/bcftools-1.9/bcftools index -f -t --threads $SLURM_CPUS_ON_NODE $targetname.raw.vcf.gz 

fi


echo $(date) done.main.process
