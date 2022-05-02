# ==========
# = QC VCF =
# ==========

# calculate individual level missing rate
# ============================================================

for chr in $(seq 1 18) X
do
  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=4:00:00 /mnt/research/qgg/software/vcftools-v0.1.16/bin/vcftools --gzvcf chr"$chr".vcf.gz --missing-indv --out chr"$chr" > chr"$chr".imiss.log 2>&1 &
done

cat chr*.imiss | awk '$1 != "INDV"' | sort -k1,1 | ~/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2,4 -o sum,sum | awk '{print $1"\t"$2"\t"$3"\t"$3/$2}' > indv.all.site.imiss

# filter animals
# ============================================================

awk '$4 > 0.20 {print $1}' indv.all.site.imiss > imiss.p20.filter.indv

for chr in $(seq 1 18) X
do
  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=8:00:00 /mnt/research/qgg/software/bcftools-1.13/bcftools view -S ^imiss.p20.filter.indv -O z -o chr"$chr".p20.vcf.gz chr"$chr".filter.vcf.gz > log/chr"$chr".p20.fitler.log 2>&1 &
done

# filter vased on missing rate per site
# ============================================================

for chr in $(seq 1 18) X
do
  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=12:00:00 /mnt/research/qgg/software/vcftools-v0.1.16/bin/vcftools --gzvcf chr"$chr".p20.vcf.gz --remove-filtered-all --min-meanDP 5 --max-meanDP 500 --max-missing 0.2 --maf 0.0001 --max-alleles 2 --recode --recode-INFO-all --out chr"$chr".p20.g80 > log/chr"$chr".g80.fitler.log 2>&1 &
done

# now compute heterozygosity and filter
# ============================================================

for chr in $(seq 1 18) X
do
  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=4:00:00 /mnt/research/qgg/software/vcftools-v0.1.16/bin/vcftools --vcf chr"$chr".p20.g80.recode.vcf --het --out chr"$chr".p20.g80 > log/chr"$chr".het.log 2>&1 &
done

# summarize
cat chr*.p20.g80.het | awk '$1 != "INDV"' | sort -k1,1 | ~/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2,4 -o sum,sum | awk '{print $1"\t"$2"\t"$3"\t"1-$2/$3}' > indv.p20.g80.het

awk '$4 > 0.20 {print $1}' indv.p20.g80.het > indv.p20.g80.het.filter.ind

# also filter a duplicate animal
echo CC3_C5_illumina >> indv.p20.g80.het.filter.ind
grep korea indv.all.site.imiss | cut -f 1 >> indv.p20.g80.het.filter.ind

# final set
for chr in $(seq 1 18) X
do
  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=8:00:00 /mnt/research/qgg/software/bcftools-1.13/bcftools view -S ^indv.p20.g80.het.filter.ind -O v -o chr"$chr".p20.g20.het.vcf chr"$chr".p20.g80.recode.vcf > log/chr"$chr".het.fitler.log 2>&1 &
done

# one last filter to eliminate monomorphic sites

for chr in $(seq 1 18) X
do
  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=12:00:00 /mnt/research/qgg/software/vcftools-v0.1.16/bin/vcftools --vcf chr"$chr".p20.g20.het.vcf --remove-filtered-all --max-missing 0.2 --maf 0.0001 --max-alleles 2 --recode --recode-INFO-all --out imputationEligible/chr"$chr".p20.g80.het.filter > log/chr"$chr".het.fitler2.log 2>&1 &
done

# get IDs
awk '$4 <= 0.20 {print $1}' indv.p20.g80.het | grep -v korea | grep -v illumina | cut -f 1 | sort -k1,1 > imputationEligible/indv.id
tail -n+2 SWIM_sampleInfo_V1.txt | sort -k1,1 | join -t $'\t' - indv.id | cut -f 1,4- > indv.info

# extract individuals from korea
# ============================================================

grep korea indv.all.site.imiss | cut -f 1 > korea.ind

# final set
for chr in $(seq 1 18) X
do
  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=32G --time=8:00:00 /mnt/research/qgg/software/bcftools-1.13/bcftools view -S korea.ind -O z -o korea/korea.chr"$chr".vcf.gz chr"$chr".p20.g80.recode.vcf > log/korea.chr"$chr".log 2>&1 &
done

# prepare sequence data
# ============================================================

tail -n+2 SWIM_sampleInfo_V1.txt | sort -k1,1 | join -t $'\t' - indv.id | awk '$5 == "Landrace"' | sort -k2,2gr | awk '$2 >= 15' | tail -n+2  | cut -f 1 > test.set.id

for sample in $(cat test.set.id)
do
  if [[ -e pigWGS/$sample ]]
  then
    cp -r pigWGS/$sample test.set.seq/
  else
    cp -r pigWGS2/$sample test.set.seq/
  fi
done


