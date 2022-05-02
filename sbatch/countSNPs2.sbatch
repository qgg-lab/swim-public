#!/bin/bash

# Job name:
#SBATCH --job-name=countSNPs2
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=9
#
# Memory per node:
#SBATCH --mem=24G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=24:00:00

# ============================================================

/mnt/research/qgg/software/plink-v1.90b6.18/plink --vcf /mnt/scratch/huangw53/swim/genotype/imputationEligible/"$chr".p20.g80.het.filter.recode.vcf --keep final.id.fam --maf 0.00001 --make-bed --out "$chr".tmp > "$chr".tmp.log 2>&1

for i in `seq 1 1 19` `seq 20 5 100` `seq 200 200 1000` 1500 2000 2259
do
  /mnt/research/qgg/software/plink-v1.90b6.18/plink --bfile "$chr".tmp --keep n"$i".random.id --freq --out "$chr".tmp.count > "$chr".tmp.count.log 2>&1
echo -ne "$i\t" >> $chr.snp.count
	tail -n+2 "$chr".tmp.count.frq | awk '$5 > 0.00001 && $5 != "NA" {print $2}' | wc -l | awk '{print $1}' >> $chr.snp.count
done

rm "$chr".tmp.*

echo $(date) done.main.process
