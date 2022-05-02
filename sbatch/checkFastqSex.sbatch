#!/bin/bash

# Job name:
#SBATCH --job-name=check-fastq-sex
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
#SBATCH --mem=8G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=1:00:00

# run main program
# ============================================================

source "$env"

# add fastq-sample to deal with some fastq downloaded from SRA that apparently are sorted based on coordiantes
# could be due to fastq generated from bam file
gunzip -c $file > $tmp/tmp.fastq
$FASTQSAMPLE -n $head $tmp/tmp.fastq -o $tmp/tmp.sample

$BWA mem -t 1 $BWAIDX $tmp/tmp.sample.fastq 2> /dev/null | cut -f 3 | cat - <(echo -e "1\nX\nY") | sort | uniq -c | awk '$2 == 1 || $2 == "X" || $2 == "Y" {print $2"\t"$1}' | cut -f 2 | paste - - - | awk '{ if ($3/$2 < 0.075) { print "female\t"$0 } else if ($3/$2 > 0.15 && $3/$2 < 0.6) { print "male\t"$0 } else { print "amb\t"$0 } }' > "$outfile"

rm $tmp/tmp.fastq $tmp/tmp.sample.fastq

echo $(date +"%m-%d-%Y-%T"): done.main.process
