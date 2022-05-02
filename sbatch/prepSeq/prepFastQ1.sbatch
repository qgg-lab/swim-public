#!/bin/bash

# Job name:
#SBATCH --job-name=prep-fastq
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

# deal with western_pig_MGISEQ_2000
# -rw-r----- 1 huangw53 qgg  5927833841 Jan 21 06:00 yorkshire_LTZ100602a_1.fq.gz
# -rw-r----- 1 huangw53 qgg  6041545448 Jan 21 06:00 yorkshire_LTZ100602a_2.fq.gz
# -rw-r----- 1 huangw53 qgg  7095380922 Jan 21 06:00 yorkshire_LTZ100602b_1.fq.gz
# -rw-r----- 1 huangw53 qgg  7194546127 Jan 21 06:00 yorkshire_LTZ100602b_2.fq.gz
# -rw-r----- 1 huangw53 qgg 13032390792 Jan 21 06:00 yorkshire_LTZ100608_1.fq.gz
# -rw-r----- 1 huangw53 qgg 13078073770 Jan 21 06:00 yorkshire_LTZ100608_2.fq.gz
# -rw-r----- 1 huangw53 qgg 11877490942 Jan 21 06:00 yorkshire_LTZ300_1.fq.gz
# -rw-r----- 1 huangw53 qgg 11924449150 Jan 21 06:00 yorkshire_LTZ300_2.fq.gz

echo "SLURM_JOB_ID: $SLURM_JOB_ID"

cd /mnt/research/qgg/share/

# get instrument info
# ============================================================

gunzip -c /mnt/research/qgg/share/"$file"_1"$add".fq.gz | paste - - - - | cut -f 1 | sed 's/^@//' | awk -F ":" '{print "pe."$1"."$2"."$3"."$4" "$1" "$2" "$3" "$4}' | uniq | sort | uniq > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/file.info

# split fastq
# ============================================================

while read line           
do           
    # output fastq file
	ins=$(echo $line | cut -d " " -f 2)
	run=$(echo $line | cut -d " " -f 3)
	flc=$(echo $line | cut -d " " -f 4)
	lane=$(echo $line | cut -d " " -f 5)
	gunzip -c /mnt/research/qgg/share/"$file"_1"$add".fq.gz | grep -A 3 ^"@$ins:$run:$flc:$lane" | gzip > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/"pe.$ins.$run.$flc.$lane"_1.fastq.gz
	gunzip -c /mnt/research/qgg/share/"$file"_2"$add".fq.gz | grep -A 3 ^"@$ins:$run:$flc:$lane" | gzip > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/"pe.$ins.$run.$flc.$lane"_2.fastq.gz
done < /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/file.info

chmod 440 /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/*.fastq.gz

sacct -j $SLURM_JOB_ID --units=G -P --format="Account%15,JobID,JobName%15,CPUTime,Elapsed,MaxRSS,Start,End,NodeList,NTasks,Partition%20,NCPUS,state"
