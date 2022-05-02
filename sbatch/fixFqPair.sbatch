#!/bin/bash

# Job name:
#SBATCH --job-name=fixFqPair
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
#SBATCH --mem=5G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=4:00:00
# ============================================================

source "$env"

while read line
do
  file=$(echo $line | cut -d " " -f 1)
  ins=$(echo $line | cut -d " " -f 2)
  run=$(echo $line | cut -d " " -f 3)
  flc=$(echo $line | cut -d " " -f 4)
  lane=$(echo $line | cut -d " " -f 5)
  rgid="$file"
	
  if [ `echo $file | grep '^pe' | wc -l` -gt 0 ]
  then
    
    # unzip and sort
    gunzip -c $dir/"$file"_1.fastq.gz | paste - - - - | sort -k1,1 > $dir/"$file"_sorted_1.fastq &
    gunzip -c $dir/"$file"_2.fastq.gz | paste - - - - | sort -k1,1 > $dir/"$file"_sorted_2.fastq &
    wait
    join -t $'\t' $dir/"$file"_sorted_1.fastq $dir/"$file"_sorted_2.fastq | cut -f 1-4 | tr '\t' '\n' | gzip > $dir/"$file"_paired_1.fastq.gz &
    join -t $'\t' $dir/"$file"_sorted_1.fastq $dir/"$file"_sorted_2.fastq | cut -f 1,5-7 | tr '\t' '\n' | gzip > $dir/"$file"_paired_2.fastq.gz &
    wait
    mv $dir/"$file"_paired_1.fastq.gz $dir/"$file"_1.fastq.gz
    mv $dir/"$file"_paired_2.fastq.gz $dir/"$file"_2.fastq.gz
    
  fi
done < $dir/file.info

echo $(date +"%m-%d-%Y-%T"):info: done.main.process
