#!/bin/bash

# Job name:
#SBATCH --job-name=phase_shapeit4
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

if [[ $sequencing!=1 ]]; then
   
	$SHAPEIT4 \
	--input $input \
	--map /mnt/ls15/scratch/users/dingrong/imputation_test/table_shapeit4/genMap_1cMperMb_"$chr".txt \
	--region $chr \
	--effective-size 100 \
	--pbwt-mdr 0 \
	--pbwt-mac 0 \
	--output "$output".vcf.gz \
	--thread $SLURM_CPUS_ON_NODE 
 
	$BCFTOOLS index -f -t --threads $SLURM_CPUS_ON_NODE "$output".vcf.gz
	
	wait

else
   
	$SHAPEIT4 \
	--input $input \
	--map /mnt/ls15/scratch/users/dingrong/imputation_test/table_shapeit4/genMap_1cMperMb_"$chr".txt \
	--region $chr \
	--effective-size 100 \
	--pbwt-mdr 0 \
	--pbwt-mac 0 \
	--output "$output".vcf.gz \
	--thread $SLURM_CPUS_ON_NODE \
	--sequencing 
	
    wait
   
	$BCFTOOLS index -f -t --threads $SLURM_CPUS_ON_NODE "$output".vcf.gz
	
	wait
 
fi

echo $(date) done.main.process
