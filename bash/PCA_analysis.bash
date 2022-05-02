# ================
# = PCA analysis =
# ================


###pruned
mkdir /mnt/ls15/scratch/users/dingrong/swim/PCA
mkdir /mnt/ls15/scratch/users/dingrong/swim/PCA/log

for chr in $(seq 1 18) 
do
  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=4:00:00 \
  /mnt/research/qgg/software/plink-v1.90b6.18/plink \
    --bfile /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/chr"$chr".p20.g80.het.filter.QC.final \
	--maf 0.05 \
	--geno 0.1 \
    --indep-pairwise 50 10 0.3 \
    --out /mnt/ls15/scratch/users/dingrong/swim/PCA/chr"$chr".ld \
    > /mnt/ls15/scratch/users/dingrong/swim/PCA/log/indep.pairwise.chr"$chr".log 2>&1 &
done

###filter pruned
for chr in $(seq 1 18) 
do
  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=4:00:00 \
  /mnt/research/qgg/software/plink-v1.90b6.18/plink \
    --bfile /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/chr"$chr".p20.g80.het.filter.QC.final \
    --extract /mnt/ls15/scratch/users/dingrong/swim/PCA/chr"$chr".ld.prune.in \
    --make-bed \
    --out /mnt/ls15/scratch/users/dingrong/swim/PCA/chr"$chr".p20.g80.het.filter.pruned \
    > /mnt/ls15/scratch/users/dingrong/swim/PCA/log/pruned.chr"$chr".log 2>&1 &
done


###merge bfile
rm /mnt/ls15/scratch/users/dingrong/swim/PCA/merge.batch.txt
for chr in $(seq 2 18) 
do
echo /mnt/ls15/scratch/users/dingrong/swim/PCA/chr"$chr".p20.g80.het.filter.pruned.bed /mnt/ls15/scratch/users/dingrong/swim/PCA/chr"$chr".p20.g80.het.filter.pruned.bim /mnt/ls15/scratch/users/dingrong/swim/PCA/chr"$chr".p20.g80.het.filter.pruned.fam >>/mnt/ls15/scratch/users/dingrong/swim/PCA/merge.batch.txt
done


/mnt/research/qgg/software/plink-v1.90b6.18/plink \
--bfile /mnt/ls15/scratch/users/dingrong/swim/PCA/chr1.p20.g80.het.filter.pruned \
--merge-list /mnt/ls15/scratch/users/dingrong/swim/PCA/merge.batch.txt \
--make-bed \
--out /mnt/ls15/scratch/users/dingrong/swim/PCA/all.p20.g80.het.filter.pruned

wc -l /mnt/ls15/scratch/users/dingrong/swim/PCA/all.p20.g80.het.filter.pruned.bim 
1223882

### pca analysis

/mnt/research/qgg/software/gcta_1.93.2beta/gcta64 \
--bfile /mnt/ls15/scratch/users/dingrong/swim/PCA/all.p20.g80.het.filter.pruned \
--autosome --make-grm \
--out /mnt/ls15/scratch/users/dingrong/swim/PCA/pca.p20.g80.het.filter.pruned.grm

/mnt/research/qgg/software/gcta_1.93.2beta/gcta64 \
--grm /mnt/ls15/scratch/users/dingrong/swim/PCA/pca.p20.g80.het.filter.pruned.grm \
--pca 5 \
--out /mnt/ls15/scratch/users/dingrong/swim/PCA/pca.p20.g80.het.filter.pruned.grm.pca5

