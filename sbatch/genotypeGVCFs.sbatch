#!/bin/bash

# Job name:
#SBATCH --job-name=genotypeGVCFs
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
#SBATCH --time=48:00:00

# ============================================================

source "$env"

java -Xmx"$mem" -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$tmp \
  -jar $GATK GenotypeGVCFs \
  -R $REFFA \
  -V $input \
  -O $dir/$prefix.vcf.gz > $dir/log/$prefix.genotype.log 2>&1 

echo $(date) done.main.process
