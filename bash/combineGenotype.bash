# ==============================
# = combine and genotype GVCFs =
# ==============================

# 1. prepare list and intervals
# hc1 = chr 1-4
# hc2 = chr 5-9
# hc3 = chr 10-14
# hc4 = chr 15-18, X (ploidy = 2)
# hc5 = chr X (plody = 1)
# ============================================================

# for both sexes + autosomes, the list are those with ploidy = 2
# for X chromosome choose ploidy based on inferred sex
# ============================================================

head -n 19 /mnt/research/qgg/resource/swim/genome.fa.fai | /mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools makewindows -w 2000000 -g /dev/stdin | awk '{print $1":"$2+1"-"$3}' > haplo.call.windows

# prepare list for hc1 - hc5
awk '{print "/mnt/ufs18/scratch/huangw53/swim/gvcfs/"$1"/hc1.p2.g.vcf.gz"}' /mnt/ufs18/scratch/huangw53/swim/gvcfs/gvcf.complete.list.07252021 > hc1.auto.file.list
awk '{print "/mnt/ufs18/scratch/huangw53/swim/gvcfs/"$1"/hc2.p2.g.vcf.gz"}' /mnt/ufs18/scratch/huangw53/swim/gvcfs/gvcf.complete.list.07252021 > hc2.auto.file.list
awk '{print "/mnt/ufs18/scratch/huangw53/swim/gvcfs/"$1"/hc3.p2.g.vcf.gz"}' /mnt/ufs18/scratch/huangw53/swim/gvcfs/gvcf.complete.list.07252021 > hc3.auto.file.list
awk '{print "/mnt/ufs18/scratch/huangw53/swim/gvcfs/"$1"/hc4.p2.g.vcf.gz"}' /mnt/ufs18/scratch/huangw53/swim/gvcfs/gvcf.complete.list.07252021 > hc4.auto.file.list
awk '{ if ($2 == "female") { print "/mnt/ufs18/scratch/huangw53/swim/gvcfs/"$1"/hc4.p2.g.vcf.gz" } else { print "/mnt/ufs18/scratch/huangw53/swim/gvcfs/"$1"/hc5.p1.g.vcf.gz" } }' /mnt/ufs18/scratch/huangw53/swim/gvcfs/gvcf.complete.list.07252021 > hc5.x.file.list

# 2. run combine on each interval
# ============================================================

# prepare jobs in a text file so they can be parallelized
perl -wne 'chomp $_; @pos = split /:/, $_; if ($pos[0] eq "X") { print join(".", @pos), "\t", $_, "\t", "hc5.x.file.list", "\n"; } elsif ($pos[0] >= 1 && $pos[0] <= 4) { print join(".", @pos), "\t", $_, "\t", "hc1.auto.file.list", "\n"; } elsif ($pos[0] >= 5 && $pos[0] <= 9) { print join(".", @pos), "\t", $_, "\t", "hc2.auto.file.list", "\n"; } elsif ($pos[0] >= 10 && $pos[0] <= 14) { print join(".", @pos), "\t", $_, "\t", "hc3.auto.file.list", "\n"; } elsif ($pos[0] >= 15 && $pos[0] <= 18) { print join(".", @pos), "\t", $_, "\t", "hc4.auto.file.list", "\n"; }' haplo.call.windows > window.parallel.07262021

cat window.parallel.07262021 | tr '\t' '\n' | ~/qgg/software/parallel-20200722/bin/parallel -N 3 -j 128 "bash /mnt/research/qgg/resource/swim/combineGVCFswait.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/ufs18/scratch/huangw53/swim/combine --tmp /mnt/ufs18/scratch/huangw53/tmp --list {3} --int {2} --pre {1} > /mnt/ufs18/scratch/huangw53/swim/combine/log/{1}.out 2>&1" > /mnt/ufs18/scratch/huangw53/swim/combine/log/combineGVCFs.07262021.log 2>&1 &

# 3. find failed runs and re-run
# ============================================================

for int in $(cut -f 1 window.parallel.07262021)
do
  if [ `grep "completed combineGVCFs" log/$int.out | wc -l` -eq 0 ] || [ `grep "error" log/$int.out | wc -l` -gt 0 ] || [ ! -f $int.g.vcf.gz.tbi ] || [ `grep "error" log/$int.combineGVCF.log | wc -l` -gt 0 ]
  then
    echo $int
  fi
done | sort > window.failed.08242021

sort -k1,1 window.parallel.07262021 | join -t $'\t' - window.failed.08242021 | tr '\t' '\n' | ~/qgg/software/parallel-20200722/bin/parallel -N 3 -j 128 "bash /mnt/research/qgg/resource/swim/combineGVCFswait.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/ufs18/scratch/huangw53/swim/combine --tmp /mnt/ufs18/scratch/huangw53/tmp --list {3} --int {2} --pre {1} > /mnt/ufs18/scratch/huangw53/swim/combine/log/{1}.out 2>&1" > /mnt/ufs18/scratch/huangw53/swim/combine/log/combineGVCFs.08242021.log 2>&1 &

# one instance failed due to out of memory
# the culprit is 9:63381312-63381314

mkdir /mnt/ufs18/scratch/huangw53/tmp/9.62000001-63381311
sbatch --export=env=/mnt/research/qgg/resource/swim/resource.env,dir=/mnt/ufs18/scratch/huangw53/swim/combine,list=hc2.auto.file.list,int=9:62000001-63381311,tmp=/mnt/ufs18/scratch/huangw53/tmp,prefix=9.62000001-63381311,mem=128G --time=48:00:00 --mem=128G --output=/mnt/ufs18/scratch/huangw53/swim/combine/log/combineGVCFs.9.62000001-63381311.out --error=/mnt/ufs18/scratch/huangw53/swim/combine/log/combineGVCFs.9.62000001-63381311.err /mnt/research/qgg/resource/swim/sbatch/combineGVCFs.sbatch

mkdir /mnt/ufs18/scratch/huangw53/tmp/9.63381315-64000000
sbatch --export=env=/mnt/research/qgg/resource/swim/resource.env,dir=/mnt/ufs18/scratch/huangw53/swim/combine,list=hc2.auto.file.list,int=9:63381315-64000000,tmp=/mnt/ufs18/scratch/huangw53/tmp,prefix=9.63381315-64000000,mem=128G --time=48:00:00 --mem=128G --output=/mnt/ufs18/scratch/huangw53/swim/combine/log/combineGVCFs.9.63381315-64000000.out --error=/mnt/ufs18/scratch/huangw53/swim/combine/log/combineGVCFs.9.63381315-64000000.err /mnt/research/qgg/resource/swim/sbatch/combineGVCFs.sbatch

# 3. summarize coverage and duplication rate
# ============================================================


# 4. genotype gvcfs
# ============================================================

cat /mnt/ufs18/scratch/huangw53/swim/combine/window.parallel.07262021 <(echo -e "9.62000001-63381311\t9:62000001-63381311\thc2.auto.file.list\n9.63381315-64000000\t9:63381315-64000000\thc2.auto.file.list") | awk '$1 != "9.62000001-64000000"' | tr '\t' '\n' | ~/qgg/software/parallel-20200722/bin/parallel -N 3 -j 128 "/mnt/research/qgg/resource/swim/genotypeGVCFswait.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/gs18/scratch/users/huangw53/genotype --tmp /mnt/ufs18/scratch/huangw53/tmp --input /mnt/ufs18/scratch/huangw53/swim/combine/{1}.g.vcf.gz --pre {1} > /mnt/gs18/scratch/users/huangw53/genotype/log/{1}.out 2>&1" > /mnt/gs18/scratch/users/huangw53/genotype/log/genotypeGVCFs.08292021.log 2>&1 &

# 5. concatenace vcfs
# ============================================================

# generate file list
for chr in $(seq 1 18) X
do
  cat /mnt/ufs18/scratch/huangw53/swim/combine/window.parallel.07262021 <(echo -e "9.62000001-63381311\t9:62000001-63381311\thc2.auto.file.list\n9.63381315-64000000\t9:63381315-64000000\thc2.auto.file.list") | awk '$1 != "9.62000001-64000000"' | awk '$2 ~ /^'$chr':/ {print "/mnt/gs18/scratch/users/huangw53/genotype/"$1".vcf.gz"}' > chr"$chr".vcf.list
done

for chr in $(seq 1 18) X
do
  sbatch --export=mem=64G,env=/mnt/research/qgg/resource/swim/resource.env,input=/mnt/gs18/scratch/users/huangw53/genotype/chr"$chr".vcf.list,output=chr"$chr".vcf.gz,tmp=/mnt/ufs18/scratch/huangw53/tmp,out=/mnt/gs18/scratch/users/huangw53/genotype/,prefix=chr"$chr" --output=/mnt/gs18/scratch/users/huangw53/genotype/log/chr"$chr".mergeVCFs.out --error=/mnt/gs18/scratch/users/huangw53/genotype/log/chr"$chr".mergeVCFs.err /mnt/research/qgg/resource/swim/sbatch/mergeVCFs.sbatch
done

for chr in $(seq 1 18) X
do
  sbatch --export=mem=32G,env=/mnt/research/qgg/resource/swim/resource.env,input=chr"$chr".vcf.gz,tmp=/mnt/ufs18/scratch/huangw53/tmp,output=chr"$chr".excessHet.vcf.gz,dir=/mnt/gs18/scratch/users/huangw53/genotype/,prefix=chr"$chr" --output=/mnt/gs18/scratch/users/huangw53/genotype/log/chr"$chr".excessHet.out --error=/mnt/gs18/scratch/users/huangw53/genotype/log/chr"$chr".excessHet.err /mnt/research/qgg/resource/swim/sbatch/excessHet.sbatch
done

for chr in $(seq 1 18) X
do
  sbatch --export=mem=32G,env=/mnt/research/qgg/resource/swim/resource.env,output=chr"$chr".sites.vcf.gz,tmp=/mnt/ufs18/scratch/huangw53/tmp,input=chr"$chr".excessHet.vcf.gz,dir=/mnt/gs18/scratch/users/huangw53/genotype/,prefix=chr"$chr" --output=/mnt/gs18/scratch/users/huangw53/genotype/log/chr"$chr".makeSitesOnly.out --error=/mnt/gs18/scratch/users/huangw53/genotype/log/chr"$chr".makeSitesOnly.err /mnt/research/qgg/resource/swim/sbatch/makeSitesOnly.sbatch
done

ls chr*.sites.vcf.gz > chr.site.vcf.list
sbatch --export=mem=64G,env=/mnt/research/qgg/resource/swim/resource.env,input=/mnt/gs18/scratch/users/huangw53/genotype/chr.site.vcf.list,output=all.site.vcf.gz,tmp=/mnt/ufs18/scratch/huangw53/tmp,out=/mnt/gs18/scratch/users/huangw53/genotype/,prefix=all.site --output=/mnt/gs18/scratch/users/huangw53/genotype/log/all.site.mergeVCFs.out --error=/mnt/gs18/scratch/users/huangw53/genotype/log/all.site.mergeVCFs.err /mnt/research/qgg/resource/swim/sbatch/mergeVCFs.sbatch

# split VCF
# ============================================================

sbatch --export=mem=64G,env=/mnt/research/qgg/resource/swim/resource.env,vcf=all.site.vcf.gz,output1=all.site.snp.vcf.gz,output2=all.site.indel.vcf.gz,tmp=/mnt/ufs18/scratch/huangw53/tmp,out=/mnt/gs18/scratch/users/huangw53/genotype/ --output=/mnt/gs18/scratch/users/huangw53/genotype/log/all.site.splitVCF.out --error=/mnt/gs18/scratch/users/huangw53/genotype/log/all.site.splitVCF.err /mnt/research/qgg/resource/swim/sbatch/splitVCF.sbatch

# VQSR SNP
# ============================================================

sbatch --export=mem=256G,env=/mnt/research/qgg/resource/swim/resource.env,input=all.site.snp.vcf.gz,tmp=/mnt/ufs18/scratch/huangw53/tmp,out=/mnt/gs18/scratch/users/huangw53/genotype/,pre=all.site.snp --output=/mnt/gs18/scratch/users/huangw53/genotype/log/all.site.snp.vqsrSNP.out --error=/mnt/gs18/scratch/users/huangw53/genotype/log/all.site.snp.vqsrSNP.err /mnt/research/qgg/resource/swim/sbatch/vqsrSNP.sbatch

sbatch --export=mem=64G,env=/mnt/research/qgg/resource/swim/resource.env,vcf=all.site.snp.vcf.gz,output=all.site.snp.vqsr.vcf.gz,tmp=/mnt/ufs18/scratch/huangw53/tmp,out=/mnt/gs18/scratch/users/huangw53/genotype/,tranches=all.site.snp.tranches,recal=all.site.snp.recal --output=/mnt/gs18/scratch/users/huangw53/genotype/log/all.site.snp.applyVQSRSNP.out --error=/mnt/gs18/scratch/users/huangw53/genotype/log/all.site.snp.applyVQSRSNP.err /mnt/research/qgg/resource/swim/sbatch/applyVQSRSNP.sbatch

# Hard filter indels
# ============================================================

sbatch --export=mem=64G,env=/mnt/research/qgg/resource/swim/resource.env,vcf=all.site.indel.vcf.gz,output=all.site.indel.hardFilter.vcf.gz,tmp=/mnt/ufs18/scratch/huangw53/tmp,out=/mnt/gs18/scratch/users/huangw53/genotype/ --output=/mnt/gs18/scratch/users/huangw53/genotype/log/all.site.indel.hardFilter.out --error=/mnt/gs18/scratch/users/huangw53/genotype/log/all.site.indel.hardFilter.err /mnt/research/qgg/resource/swim/sbatch/hardFilterIndel.sbatch

# Combine SNPs and indels
# ============================================================

echo -e "all.site.snp.vqsr.vcf.gz\nall.site.indel.hardFilter.vcf.gz" > all.site.filter.vcf.list

sbatch --export=mem=64G,env=/mnt/research/qgg/resource/swim/resource.env,input=/mnt/gs18/scratch/users/huangw53/genotype/all.site.filter.vcf.list,output=all.site.add.filter.vcf.gz,tmp=/mnt/ufs18/scratch/huangw53/tmp,out=/mnt/gs18/scratch/users/huangw53/genotype/,prefix=all.site.add.filter --output=/mnt/gs18/scratch/users/huangw53/genotype/log/all.site.add.filter.mergeVCFs.out --error=/mnt/gs18/scratch/users/huangw53/genotype/log/all.site.add.filter.mergeVCFs.err /mnt/research/qgg/resource/swim/sbatch/mergeVCFs.sbatch

# Filter only PASS and remove overlapping variants
# ============================================================

srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=4:00:00 /mnt/research/qgg/software/bcftools-1.13/bcftools view -f 'PASS' -O z -o all.site.pass.vcf.gz all.site.add.filter.vcf.gz > all.site.pass.vcf.log 2>&1 &

gunzip -c all.site.pass.vcf.gz | perl filterOverlappingVCFRef.pl > all.site.pass.no.overlap.vcf 2> all.site.pass.overlap.vcf &

# remove PAR region 1M-6.5M on X chromosome
# NCBI says this: PAR: X:1029482 - 6405430；Y:200470 - 4791971
# ============================================================

srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=4:00:00 /mnt/research/qgg/software/bcftools-1.13/bcftools view -t ^X:1000001-6500000 -O v -o all.site.pass.no.overlap.no.PAR.vcf all.site.pass.no.overlap.vcf > log/all.site.pass.no.overlap.no.PAR.log 2>&1 &

# extract pass SNPs
# ============================================================

for chr in $(seq 1 18) X
do
  sbatch --export=env=/mnt/research/qgg/resource/swim/resource.env,output=chr"$chr".filter.vcf.gz,file1=chr"$chr".vcf.gz,varList=/mnt/gs18/scratch/users/huangw53/genotype/all.site.pass.no.overlap.no.PAR.vcf --time=24:00:00 --output=/mnt/gs18/scratch/users/huangw53/genotype/log/chr"$chr".extractVCF.out --error=/mnt/gs18/scratch/users/huangw53/genotype/log/chr"$chr".extractVCF.err /mnt/research/qgg/resource/swim/sbatch/extractVCF.sbatch
done
