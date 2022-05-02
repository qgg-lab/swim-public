#!/bin/bash

# Job name:
#SBATCH --job-name=cov
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
#SBATCH --mem=4G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=4:00:00

source "$env"

$SAMTOOLS view -b -q 50 -F 0x400 $input 2>> $out/log/cov.log | $BEDTOOLS genomecov -ibam - > $output 2>> $out/log/cov.log

echo $(date) done.main.process
