#!/bin/bash

# Job name:
#SBATCH --job-name=finalQC
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
#SBATCH --mem=64G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=23:00:00

# ============================================================

 /mnt/research/qgg/software/vcftools-v0.1.16/bin/vcftools \
  --vcf /mnt/gs18/scratch/users/huangw53/share/imputationEligible/chr"$chr".p20.g80.het.filter.recode.vcf \
  --snps /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/chr"$chr".QC.final.sites \
  --remove /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/remove.82.ind \
  --recode \
  --recode-INFO-all \
  --out /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/chr"$chr".QC.final \
  > /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/finalQC.chr"$chr".log
  
wait

sed 's@\t.:@\t./.:@g'  /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/chr"$chr".QC.final.recode.vcf > /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/chr"$chr".QC.final.vcf	

wait

bcftools view /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/chr"$chr".QC.final.vcf -Oz -o /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/chr"$chr".QC.final.vcf.gz > /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/chr"$chr".gz.log

wait

bcftools index -f -t --threads 2 /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/chr"$chr".QC.final.vcf.gz

echo $(date) done.main.process   
