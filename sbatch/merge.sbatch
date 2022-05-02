#!/bin/bash

# Job name:
#SBATCH --job-name=merge
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

awk '{print "'$tmp'/"$1".map.bam"}' $dir/file.info | xargs $SAMTOOLS merge $tmp/sort.bam 2> $out/log/merge.log

$SAMTOOLS index $tmp/sort.bam

echo $(date) done.main.process
