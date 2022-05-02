#!/bin/bash

# Job name:
#SBATCH --job-name=bwa
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
		echo -n $file >> $out/seq.sum
		gunzip -c $dir/"$file"_1.fastq.gz | paste - - - - | cut -f 2 | awk '{ sum += length($1) } END { print " "NR" "sum/NR}' | tr -d '\n' >> $out/seq.sum
		gunzip -c $dir/"$file"_2.fastq.gz | paste - - - - | cut -f 2 | awk '{ sum += length($1) } END { print " "NR" "sum/NR}' >> $out/seq.sum
    $BWA mem -M -t 8 -R "@RG\tID:$rgid\tPL:illumina\tLB:$sample\tSM:$sample\tPU:$flc.$lane" $BWAIDX $dir/"$file"_1.fastq.gz $dir/"$file"_2.fastq.gz 2> $out/log/$file.mem.log | $SAMTOOLS sort -T $tmp/$file.map.bam - > $tmp/$file.map.bam 2> $out/log/$file.samsort.log
  else
		echo -n $file >> $out/seq.sum
		gunzip -c $dir/"$file"_1.fastq.gz | paste - - - - | cut -f 2 | awk '{ sum += length($1) } END { print " "NR" "sum/NR" 0 0"}' >> $out/seq.sum
    $BWA mem -M -t 8 -R "@RG\tID:$rgid\tPL:illumina\tLB:$sample\tSM:$sample\tPU:$flc.$lane" $BWAIDX $dir/"$file"_1.fastq.gz 2> $out/log/$file.mem.log | $SAMTOOLS sort -T $tmp/$file.map.bam - > $tmp/$file.map.bam 2> $out/log/$file.samsort.log
  fi
done < $dir/file.info

echo $(date) done.main.process
