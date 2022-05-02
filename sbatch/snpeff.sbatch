#!/bin/bash

# Job name:
#SBATCH --job-name=snpeff
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
#SBATCH --mem=64G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=11:50:58

# ============================================================


cd /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/

java -jar /mnt/research/qgg/software/snpEff-5.0e/snpEff.jar \
     -config /mnt/research/qgg/software/snpEff-5.0e/snpEff.config \
     -classic -lof -o vcf -ud 5000 -s snpeff.summary.html \
     ssc11.1 chr"$chr".QC.final.vcf \
     > annotation/chr"$chr".QC.final.snpeff.vcf \
     2> annotation/snpeff."$chr".log

for chr in $(seq 1 18)
do
echo $chr
grep -v ^# annotation/chr"$chr".QC.final.snpeff.vcf|awk '{print $1,$2,$3,$4,$5,$8}' >> annotation/swim.annotation.vcf
done

echo $(date) done.main.process

