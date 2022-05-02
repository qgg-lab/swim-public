#!/bin/bash

# Job name:
#SBATCH --job-name=imputation_beagle5.2
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

java -Xmx"$mem" -Djava.io.tmpdir=$tmp \
-jar  $BEAGLE \
ref="$ref_panel".vcf.gz \
gt=$input.vcf.gz \
out=$output \
window=100 \
ne=100 \
nthreads=$SLURM_CPUS_ON_NODE 

$BCFTOOLS index -f -t --threads $SLURM_CPUS_ON_NODE $output.vcf.gz
   
echo $(date) done.main.process
