#!/bin/bash

# Job name:
#SBATCH --job-name=split-fastq-bgi
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
#SBATCH --mem=8G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=24:00:00

# ============================================================

echo "SLURM_JOB_ID: $SLURM_JOB_ID"

# prepare files, convert it into illumina
# ============================================================

echo $file1 | sed 's/|/ /' | xargs zcat | paste - - - - | perl -wne 'chomp $_; @line = split /\t/, $_; if ($line[0] =~ m/^@(.*)L(.*)C(.{3})R(.{3})(.*)\/(.?)/) { print "\@MGISEQ:-:", $1, ":", $2, ":", $5 + 1, ":", $3+0, ":", $4+0, "/", $6, "\n", join("\n", @line[1..3]), "\n"; }' > $tmp/"$sample"_tmp_1.fq &
echo $file2 | sed 's/|/ /' | xargs zcat | paste - - - - | perl -wne 'chomp $_; @line = split /\t/, $_; if ($line[0] =~ m/^@(.*)L(.*)C(.{3})R(.{3})(.*)\/(.?)/) { print "\@MGISEQ:-:", $1, ":", $2, ":", $5 + 1, ":", $3+0, ":", $4+0, "/", $6, "\n", join("\n", @line[1..3]), "\n"; }' > $tmp/"$sample"_tmp_2.fq &
wait

file1="$tmp"/"$sample"_tmp_1.fq
file2="$tmp"/"$sample"_tmp_2.fq

# check if there is unpaired lines
# ============================================================

unpline=`paste <(cat "$file1" | paste - - - - | cut -f 1 | sed 's/\/.*//' | awk '{print $1}') <(cat "$file2" | paste - - - - | cut -f 1 | sed 's/\/.*//' | awk '{print $1}') | awk '$1 != $2' | head | wc -l`

if [[ $unpline -gt 0 ]]
then
	echo $(date +"%m-%d-%Y-%T"):error: non-matching reads
	exit
fi
	
# get instrument info
# ============================================================

cat "$file1" | paste - - - - | cut -f 1 | sed 's/^@//' | awk -F ":" '{print "pe."$1"."$2"."$3"."$4" "$1" "$2" "$3" "$4}' | uniq | sort | uniq > "$dir"/file.info

# split fastq
# ============================================================

while read line           
do           
  # output fastq file
	ins=$(echo $line | cut -d " " -f 2)
	run=$(echo $line | cut -d " " -f 3)
	flc=$(echo $line | cut -d " " -f 4)
	lane=$(echo $line | cut -d " " -f 5)
	grep --no-group-separator -A 3 ^"@$ins:$run:$flc:$lane" $file1 > $dir/"pe.$ins.$run.$flc.$lane"_1.fastq &
	grep --no-group-separator -A 3 ^"@$ins:$run:$flc:$lane" $file2 > $dir/"pe.$ins.$run.$flc.$lane"_2.fastq &
  wait
  gzip $dir/"pe.$ins.$run.$flc.$lane"_1.fastq &
  gzip $dir/"pe.$ins.$run.$flc.$lane"_2.fastq &
  wait
done < "$dir"/file.info

chmod 440 $dir/*.fastq.gz
rm $tmp/"$sample"_tmp_1.fq $tmp/"$sample"_tmp_2.fq

echo $(date +"%m-%d-%Y-%T"):info: done.main.process
