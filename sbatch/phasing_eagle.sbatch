#!/bin/bash

# Job name:
#SBATCH --job-name=phase_eagle
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

$EAGLE \
--geneticMapFile /mnt/research/qgg/software/Eagle_v2.4/tables/genetic_map_1cMperMb.txt \
--vcf $input \
--numThreads $SLURM_CPUS_ON_NODE \
--vcfOutFormat=z \
--outPrefix $output 

$BCFTOOLS index -f -t --threads $SLURM_CPUS_ON_NODE "$output".vcf.gz


if [[ $sequencing == 1 ]]; then
   
$MINIMAC3 \
--refHaps "$output".vcf.gz \
--processReference  \
--prefix $output 
	
	wait
fi

$BCFTOOLS index -f -t --threads $SLURM_CPUS_ON_NODE "$output".m3vcf.gz

echo $(date) done.main.process

