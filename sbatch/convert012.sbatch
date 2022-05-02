#!/bin/bash

# Job name:
#SBATCH --job-name=convert012
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
#SBATCH --mem=32G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=1:58:00

# ============================================================

if [[ $imputeSof == "minimac4" ]]; then

# Imputed set
echo "PLINK: Converting vcf to 012 format (imputation set)"

/mnt/research/qgg/software/plink-v1.90b6.18/plink  \
   --vcf $filepath/cv"$split"/target.chr"$chr".minimac4.imputation.dose.vcf.gz \
   --threads $THREADS \
   --indiv-sort f $filepath/cv"$split"/012/sampleNamesPred.txt \
   --extract $filepath/cv"$split"/012/pos_minimac4.chr"$chr".txt \
   --double-id \
   --keep-allele-order --make-bed \
   --out $filepath/cv"$split"/012/plinkPhasedPred_minimac4.chr"$chr" > /dev/null 2>&1 
	wait
	
/mnt/research/qgg/software/plink-v1.90b6.18/plink \
   --bfile $filepath/cv"$split"/012/plinkPhasedPred_minimac4.chr"$chr" \
   --threads $THREADS \
   --keep-allele-order --recode A \
   --out $filepath/cv"$split"/012/phasedPred012_minimac4.chr"$chr" > /dev/null 2>&1 
wait

elif [[ $imputeSof == "impute5" ]]; then

# Imputed set
echo "PLINK: Converting vcf to 012 format (imputation set)"

/mnt/research/qgg/software/plink-v1.90b6.18/plink  \
   --vcf $filepath/cv"$split"/target.chr"$chr".impute5.imputation.vcf.gz \
   --threads $THREADS \
   --indiv-sort f $filepath/cv"$split"/012/sampleNamesPred.txt \
   --extract $filepath/cv"$split"/012/pos_impute5.chr"$chr".txt \
   --double-id \
   --keep-allele-order --make-bed \
   --out $filepath/cv"$split"/012/plinkPhasedPred_impute5.chr"$chr" > /dev/null 2>&1 
	wait
	
/mnt/research/qgg/software/plink-v1.90b6.18/plink \
   --bfile $filepath/cv"$split"/012/plinkPhasedPred_impute5.chr"$chr" \
   --threads $THREADS \
   --keep-allele-order --recode A \
   --out $filepath/cv"$split"/012/phasedPred012_impute5.chr"$chr" > /dev/null 2>&1 
wait


elif [[ $imputeSof == "beagle" ]]; then

# Imputed set
echo "PLINK: Converting vcf to 012 format (imputation set)"

/mnt/research/qgg/software/plink-v1.90b6.18/plink  \
   --vcf $filepath/cv"$split"/target.chr"$chr".beagle.imputation.vcf.gz \
   --threads $THREADS \
   --indiv-sort f $filepath/cv"$split"/012/sampleNamesPred.txt \
   --extract $filepath/cv"$split"/012/pos_beagle.chr"$chr".txt \
   --double-id \
   --keep-allele-order --make-bed \
   --out $filepath/cv"$split"/012/plinkPhasedPred_beagle.chr"$chr" > /dev/null 2>&1 
	wait
	
/mnt/research/qgg/software/plink-v1.90b6.18/plink \
   --bfile $filepath/cv"$split"/012/plinkPhasedPred_beagle.chr"$chr" \
   --threads $THREADS \
   --keep-allele-order --recode A \
   --out $filepath/cv"$split"/012/phasedPred012_beagle.chr"$chr" > /dev/null 2>&1 
wait


fi


#validation set
echo "PLINK: Converting vcf to 012 format (validation set)"

/mnt/research/qgg/software/plink-v1.90b6.18/plink \
   --vcf $filepath/cv"$split"/target.chr"$chr".raw.vcf.gz \
   --threads $THREADS \
   --indiv-sort f $filepath/cv"$split"/012/sampleNamesPred.txt \
   --double-id \
   --extract $filepath/cv"$split"/012/pos_impute5.chr"$chr".txt \
   --keep-allele-order --make-bed \
   --out $filepath/cv"$split"/012/plinkPhasedObs.chr"$chr" > /dev/null 2>&1  
  
wait

/mnt/research/qgg/software/plink-v1.90b6.18/plink \
--bfile $filepath/cv"$split"/012/plinkPhasedObs.chr"$chr" \
   --threads $THREADS \
   --keep-allele-order --recode A \
   --out $filepath/cv"$split"/012/phasedObs012.chr"$chr" > /dev/null 2>&1 
wait
  
echo $(date) done.main.process   
