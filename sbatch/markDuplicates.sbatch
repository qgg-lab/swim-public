#!/bin/bash

# Job name:
#SBATCH --job-name=mark-duplicates
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
#SBATCH --mem=64G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=4:00:00

# ============================================================

source "$env"

# the optical duplicate pixel distance is set to 2500 for hiseq and later illumina platforms
# given the dominance of this platform, it probably is OK even if unknown platforms are used

java -Xmx"$mem" -XX:ParallelGCThreads=2 -jar $PICARD MarkDuplicates \
  input=$input \
  output=$output \
  validation_stringency=LENIENT \
  OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
  metrics_file=$metric \
  remove_duplicates=false \
  assume_sorted=true \
  max_file_handles_for_read_ends_map=512 \
  create_index=true \
  tmp_dir=$tmp > $out/log/rmdup.log 2>&1

echo $(date) done.main.process
