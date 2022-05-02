#!/bin/bash

# Job name:
#SBATCH --job-name=vqsrSNP
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
#SBATCH --mem=256G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=4:00:00

# ============================================================

source "$env"

java -Xmx"$mem" -XX:ParallelGCThreads=2 -Djava.io.tmpdir=$tmp \
  -jar $GATK VariantRecalibrator \
  -V $input \
  -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 \
  -tranche 98.0 -tranche 97.0 -tranche 95.0 -tranche 90.0 \
  -R $REFFA \
  -mode SNP \
  -resource:dbsnp,known=true,training=false,truth=false,prior=5.0 $DBSNPSNV \
  -resource:dbsnp,known=false,training=true,truth=true,prior=12.0 $AFFY660K \
  -resource:dbsnp,known=false,training=true,truth=true,prior=15.0 $AFFYSP90 \
  -resource:dbsnp,known=false,training=true,truth=true,prior=15.0 $ILMN50K \
  -resource:dbsnp,known=false,training=true,truth=true,prior=15.0 $ILMN60K \
  -resource:dbsnp,known=false,training=true,truth=true,prior=15.0 $ILMN80K \
  -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
  --max-gaussians 6 \
  -O $out/"$pre".recal \
  --tranches-file $out/"$pre".tranches \
  --rscript-file $out/"$pre".tranches.plot.R > $out/log/snp.vqsr.log 2>&1	

echo $(date) done.main.process
