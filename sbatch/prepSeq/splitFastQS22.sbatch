#!/bin/bash

# Job name:
#SBATCH --job-name=split-fastq-s22
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
#SBATCH --time=10:00:00

# ============================================================

echo "SLURM_JOB_ID: $SLURM_JOB_ID"

cd /mnt/gs18/scratch/users/huangw53/seq/s22/

# get instrument info
# ============================================================

gunzip -c /mnt/gs18/scratch/users/huangw53/seq/s22/"$sample"_1"$add".fq.gz | grep ^@ | sed 's/^@//' | awk -F ":" '{print "pe."$1"."$2"."$3"."$4" "$1" "$2" "$3" "$4}' | uniq | sort | uniq > /mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/file.info

# split fastq
# ============================================================

while read line           
do           
    # output fastq file
	ins=$(echo $line | cut -d " " -f 2)
	run=$(echo $line | cut -d " " -f 3)
	flc=$(echo $line | cut -d " " -f 4)
	lane=$(echo $line | cut -d " " -f 5)
	gunzip -c /mnt/gs18/scratch/users/huangw53/seq/s22/"$sample"_1"$add".fq.gz | grep -A 3 ^"@$ins:$run:$flc:$lane" | gzip > /mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/"pe.$ins.$run.$flc.$lane"_1.fastq.gz
	gunzip -c /mnt/gs18/scratch/users/huangw53/seq/s22/"$sample"_2"$add".fq.gz | grep -A 3 ^"@$ins:$run:$flc:$lane" | gzip > /mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/"pe.$ins.$run.$flc.$lane"_2.fastq.gz
done < /mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/file.info

chmod 440 /mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/*.fastq.gz

sacct -j $SLURM_JOB_ID --units=G -P --format="Account%15,JobID,JobName%15,CPUTime,Elapsed,MaxRSS,Start,End,NodeList,NTasks,Partition%20,NCPUS,state"
