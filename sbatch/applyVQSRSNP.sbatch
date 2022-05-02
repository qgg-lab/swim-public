#!/bin/bash

# Job name:
#SBATCH --job-name=applyVQSRSNP
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
#SBATCH --time=4:00:00

# ============================================================

source "$env"

java -Xmx"$mem" -XX:ParallelGCThreads=2 -Djava.io.tmpdir=$tmp \
  -jar $GATK ApplyVQSR \
  -V $vcf \
  -O $output \
  -R $REFFA \
  --truth-sensitivity-filter-level 99.0 \
  --tranches-file $tranches \
  --recal-file $recal \
  -mode SNP > $out/log/snp.applyVQSR.log 2>&1	

echo $(date) done.main.process
