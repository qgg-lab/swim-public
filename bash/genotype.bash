# ===========================
# = genotype across samples =
# ===========================

# in /mnt/home/huangw53/scratch2/swim/gvcfs/combine
# ============================================================

mkdir /mnt/home/huangw53/scratch2/swim/gvcfs/combine/log


awk '{print "/mnt/ufs18/scratch/huangw53/swim/gvcfs/"$1"/hc1.p2.g.vcf.gz"}' swim.release.2020.1.list > test.list

sbatch --export=env=/mnt/research/qgg/resource/swim/resource.env,mem=32G,tmp=/mnt/ufs18/scratch/huangw53/tmp,list=test.list,int=1,prefix=chr1,dir=/mnt/ufs18/scratch/huangw53/swim/gvcfs/combine --output=/mnt/ufs18/scratch/huangw53/swim/gvcfs/combine/log/combineGVCFs.out --error=/mnt/ufs18/scratch/huangw53/swim/gvcfs/combine/log/combineGVCFs.err /mnt/research/qgg/resource/swim/sbatch/combineGVCFs.sbatch

sbatch --export=env=/mnt/research/qgg/resource/swim/resource.env,mem=32G,tmp=/mnt/ufs18/scratch/huangw53/tmp,input=/mnt/ufs18/scratch/huangw53/swim/gvcfs/combine/chr1.g.vcf.gz,dir=/mnt/ufs18/scratch/huangw53/swim/gvcfs/combine,prefix=chr1 --output=/mnt/ufs18/scratch/huangw53/swim/gvcfs/combine/log/genotypeGVCFs.out --error=/mnt/ufs18/scratch/huangw53/swim/gvcfs/combine/log/genotypeGVCFs.err /mnt/research/qgg/resource/swim/sbatch/genotypeGVCFs.sbatch



ls /mnt/ufs18/scratch/huangw53/pigWGS/ | ~/qgg/software/parallel-20200722/bin/parallel -j 48 "bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/ufs18/scratch/huangw53/pigWGS --tmp /mnt/ufs18/scratch/huangw53/tmp --out /mnt/ufs18/scratch/huangw53/swim/gvcfs/ --sample {} > /mnt/ufs18/scratch/huangw53/swim/gvcfs/{}.out 2>&1" &

# when it's done, run a test to see if it is successful
# ============================================================

for sample in $(ls /mnt/ufs18/scratch/huangw53/pigWGS/)
do
  if [ `grep error $sample.out | wc -l` -gt 0 ] || [ `grep haplotype $sample.out | wc -l` -eq 0 ]
  then
    rm -r $sample
    echo $sample
  fi
done > 20201010.rerun.sample

cat 20201010.rerun.sample | ~/qgg/software/parallel-20200722/bin/parallel -j 48 "bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/ufs18/scratch/huangw53/pigWGS --tmp /mnt/ufs18/scratch/huangw53/tmp --out /mnt/ufs18/scratch/huangw53/swim/gvcfs/ --sample {} > /mnt/ufs18/scratch/huangw53/swim/gvcfs/{}.out 2>&1" &

# re-run ambiguous after updating the sex check sbatch
# ============================================================

for sample in $(ls *.out | sed 's/\.out//')
do
  if [ `grep amb $sample.out | wc -l` -gt 0 ] || [ `grep error $sample.out | wc -l` -gt 0 ] || [ `grep "completed haplotype" $sample.out | wc -l` -eq 0 ]
  then
    echo $sample
    rm -r $sample
    rm $sample.out
  else
    echo $sample >> gvcf.success.txt
  fi
done

ls /mnt/ufs18/scratch/huangw53/pigWGS/ | grep "yorkshire\|S21_\|S22_" | sort | comm -23 - <(sort gvcf.success.txt) > 20201214.run.sample

cat 20201214.run.sample | ~/qgg/software/parallel-20200722/bin/parallel -j 16 "bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/ufs18/scratch/huangw53/pigWGS --tmp /mnt/ufs18/scratch/huangw53/tmp --out /mnt/ufs18/scratch/huangw53/swim/gvcfs/ --sample {} > /mnt/ufs18/scratch/huangw53/swim/gvcfs/{}.out 2>&1" &

