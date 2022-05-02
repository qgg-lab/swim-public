mkdir /mnt/ls15/scratch/users/dingrong/swim/imputationEligible
mkdir /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/log

###Convert to plink format

for chr in $(seq 1 18) X
do
  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=64G --time=2:00:00 \
  /mnt/research/qgg/software/vcftools-v0.1.16/bin/vcftools \
  --vcf /mnt/gs18/scratch/users/huangw53/share/imputationEligible/chr"$chr".p20.g80.het.filter.recode.vcf \
  --remove /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/remove.82.ind \
  --max-alleles 2 \
  --plink-tped \
  --out /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/chr"$chr".p20.g80.het.filter.plink \
  > /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/log/vcf2plink.chr"$chr".log 2>&1 &
done

###QC1
for chr in $(seq 1 18) X
do
  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=4:00:00 \
  /mnt/research/qgg/software/plink-v1.90b6.18/plink \
    --tfile /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/chr"$chr".p20.g80.het.filter.plink \
	--maf 0.005 \
	--geno 0.1 \
	--keep-allele-order \
	--make-bed \
	--out /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/chr"$chr".p20.g80.het.filter.QC1 \
    > /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/log/chr"$chr".filter.QC1.log 2>&1 &
done

cat /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/chr*.p20.g80.het.filter.QC1.bim |wc -l
#35703417

###QC2
awk '{if ($15=="DRC") print $1,$1 }' /mnt/ls15/scratch/users/dingrong/swim/SWIM_sampleInfo_2259.txt \
 > /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/Duroc.485.ind
awk '{if ($15=="LW") print $1,$1 }' /mnt/ls15/scratch/users/dingrong/swim/SWIM_sampleInfo_2259.txt \
> /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/Large_White.543.ind
awk '{if ($15=="LR") print $1,$1 }' /mnt/ls15/scratch/users/dingrong/swim/SWIM_sampleInfo_2259.txt \
> /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/Landrace.651.ind

for chr in $(seq 1 18) X
do
  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=2:00:00 \
  /mnt/research/qgg/software/plink-v1.90b6.18/plink \
    --bfile /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/chr"$chr".p20.g80.het.filter.QC1 \
	--keep /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/Duroc.485.ind \
	--hardy \
	--out /mnt/ls15/scratch/users/dingrong/swim/Duroc/duroc.chr"$chr".hardy \
    > /mnt/ls15/scratch/users/dingrong/swim/Duroc/log/duroc.chr"$chr".hardy.log 2>&1 &
	
	srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=2:00:00 \
  /mnt/research/qgg/software/plink-v1.90b6.18/plink \
    --bfile /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/chr"$chr".p20.g80.het.filter.QC1 \
	--keep /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/Large_White.543.ind \
	--hardy \
	--out /mnt/ls15/scratch/users/dingrong/swim/Large_White/Large_White.chr"$chr".hardy \
    > /mnt/ls15/scratch/users/dingrong/swim/Large_White/log/Large_White.chr"$chr".hardy.log 2>&1 &
	
	srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=2:00:00 \
  /mnt/research/qgg/software/plink-v1.90b6.18/plink \
    --bfile /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/chr"$chr".p20.g80.het.filter.QC1 \
	--keep /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/Landrace.651.ind \
	--hardy \
	--out /mnt/ls15/scratch/users/dingrong/swim/Landrace/Landrace.chr"$chr".hardy \
    > /mnt/ls15/scratch/users/dingrong/swim/Landrace/log/Landrace.chr"$chr".hardy.log 2>&1 &

done

###get DLY populations overlap SNP( hew <0.0000000001)
cat /mnt/ls15/scratch/users/dingrong/swim/Duroc/duroc.chr*.hardy.hwe |awk '{if ($9<0.0000000001) print $2 }' \
> /mnt/ls15/scratch/users/dingrong/swim/Duroc/duroc.hweout.txt

cat /mnt/ls15/scratch/users/dingrong/swim/Large_White/Large_White.chr*.hardy.hwe |awk '{if ($9<0.0000000001) print $2 }' \
 >/mnt/ls15/scratch/users/dingrong/swim/Large_White/Large_White.hweout.txt

cat /mnt/ls15/scratch/users/dingrong/swim/Landrace/Landrace.chr*.hardy.hwe |awk '{if ($9<0.0000000001) print $2 }' \
>/mnt/ls15/scratch/users/dingrong/swim/Landrace/Landrace.hweout.txtLandrace.hweout.txt

cat /mnt/ls15/scratch/users/dingrong/swim/Duroc/duroc.chr*.hardy.hwe /mnt/ls15/scratch/users/dingrong/swim/Large_White/Large_White.chr*.hardy.hwe /mnt/ls15/scratch/users/dingrong/swim/Landrace/Landrace.chr*.hardy.hwe \
|awk '{if ($9<0.0000000001) print $2 }' |sort | uniq -c |sort -n -r | awk '{if ($1==3) print $2}' >/mnt/ls15/scratch/users/dingrong/swim/imputationEligible/DLY.hwe.snp.list


for chr in $(seq 1 18) X
do
	srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=1:50:00 \
	/mnt/research/qgg/software/plink-v1.90b6.18/plink \
    --bfile /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/chr"$chr".p20.g80.het.filter.QC1 \
	--exclude /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/DLY.hwe.snp.list \
	--keep-allele-order \
	--make-bed \
	--out /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/chr"$chr".p20.g80.het.filter.QC.final \
    > /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/log/chr"$chr".filter.QC2.log 2>&1 &
done

cat /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/chr*[0-9].p20.g80.het.filter.QC.final.bim |wc -l
#35615773

###get final imputationEligibleQCvcf
mkdir /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf

for chr in $(seq 1 18) X
do
echo $chr
awk '{ print $2}' /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/chr"$chr".p20.g80.het.filter.QC.final.bim \
>/mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/chr"$chr".QC.final.sites
done

cd /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/

for chr in $(seq 1 18) 
do
sbatch \
--time=23:59:00 \
--export=chr="$chr" \
--output=log/chr"$chr".finalQC.out \
--error=log/chr"$chr".finalQC.err \
finalQC.sbatch
done
