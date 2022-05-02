#!/bin/bash

# Job name:
#SBATCH --job-name=haplotype-caller
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=6
#
# Memory per node:
#SBATCH --mem=60G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=32:00:00

# ============================================================

source "$env"

java -Xmx"$mem" -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$tmp \
  -jar $GATK HaplotypeCaller \
  -L $HC1 \
  -I $input \
  -R $REFFA \
  -ERC GVCF \
  -ploidy 2 \
  -O $output1 > $out/log/hc1.log 2>&1 &

java -Xmx"$mem" -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$tmp \
  -jar $GATK HaplotypeCaller \
  -I $input \
  -R $REFFA \
  -L $HC2 \
  -ERC GVCF \
  -ploidy 2 \
  -O $output2 > $out/log/hc2.log 2>&1 &

java -Xmx"$mem" -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$tmp \
  -jar $GATK HaplotypeCaller \
  -I $input \
  -R $REFFA \
  -L $HC3 \
  -ERC GVCF \
  -ploidy 2 \
  -O $output3 > $out/log/hc3.log 2>&1 &

java -Xmx"$mem" -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$tmp \
  -jar $GATK HaplotypeCaller \
  -I $input \
  -R $REFFA \
  -L $HC4 \
  -ERC GVCF \
  -ploidy 2 \
  -O $output4 > $out/log/hc4.log 2>&1 &

java -Xmx"$mem" -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$tmp \
  -jar $GATK HaplotypeCaller \
  -I $input \
  -R $REFFA \
  -L $HC5 \
  -ERC GVCF \
  -ploidy 1 \
  -O $output5 > $out/log/hc5.log 2>&1 &

wait # this is necessary because both processes need to complete for the outside call to check on logs

echo $(date) done.main.process
