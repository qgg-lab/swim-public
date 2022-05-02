###select individual and QC

mkdir /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture
mkdir /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/log/

for chr in $(seq 1 18) 
do
  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=4:00:00 \
  /mnt/research/qgg/software/plink-v1.90b6.18/plink \
    --bfile /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/chr"$chr".p20.g80.het.filter.QC.final \
	--keep /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/admixture.select.ind \
	--maf 0.05 \
	--geno 0.1 \
	--make-bed \
    --out /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/select.chr"$chr".ind.QC \
    > /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/log/select.chr"$chr".ind.QC.log 2>&1 &
	
done

###pruned
for chr in $(seq 1 18) 
do
  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=4:00:00 \
  /mnt/research/qgg/software/plink-v1.90b6.18/plink \
    --bfile /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/select.chr"$chr".ind.QC \
    --indep-pairwise 50 10 0.3 \
    --out /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/select.chr"$chr".ld \
    > /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/log/select.indep.pairwise.chr"$chr".log 2>&1 &

done

for chr in $(seq 1 18) 
do
  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=4:00:00 \
  /mnt/research/qgg/software/plink-v1.90b6.18/plink \
    --bfile /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/select.chr"$chr".ind.QC \
    --extract /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/select.chr"$chr".ld.prune.in \
    --make-bed \
    --out /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/select.chr"$chr".pruned \
    > /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/log/pruned.chr"$chr".log 2>&1 &

 
done


###merge bfile
rm /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/merge.batch.txt
for chr in $(seq 2 18) 
do
echo /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/select.chr"$chr".pruned.bed /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/select.chr"$chr".pruned.bim /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/select.chr"$chr".pruned.fam >>/mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/merge.batch.txt
done


/mnt/research/qgg/software/plink-v1.90b6.18/plink \
--bfile /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/select.chr1.pruned \
--merge-list /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/merge.batch.txt \
--make-bed \
--out /mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/select.all.pruned

####admixture
for K in `seq 2 12`
do 
sbatch \
--time=23:58:00 \
--cpus-per-task=24 \
--export=K=$K \
--output=/mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/log/admixture.$K.out \
--error=/mnt/ls15/scratch/users/dingrong/swim/PCA/admixture/log/admixture.$K.err \
/mnt/ls15/scratch/users/dingrong/swim/admixture.sbatch
done


