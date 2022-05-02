# =======================
# = saturation analysis =
# =======================

# 1. get list of IDs used in the final imputation dataset using plink
# ============================================================
/mnt/research/qgg/software/plink-v1.90b6.18/plink --vcf /mnt/research/qgg/share/swim/imputationQCvcf/chr18.QC.final.shapeit4.phased.vcf.gz --make-bed --out chr18.tmp > chr.tmp.log 2>&1 &

# 2. order of ids
# ============================================================

cut -d " " -f 1,2 chr18.tmp.fam > final.id.fam

gunzip -c /mnt/research/qgg/share/swim/imputationQCvcf/chr18.QC.final.shapeit4.phased.vcf.gz | head -n 10 | tail -n 1 | cut -f 10- | tr '\t' '\n' | awk '{print $1"\t"NR}' | sort -k1,1 | join -t $'\t' - <(tail -n+2 SWIM_BreedInfo.txt | cut -f 1-2 | sort -k1,1) | sort -k2,2n > final.id.info

# prepare order based on breed numbers
echo -e "Landrace\t1\nLarge_White\t2\nDuroc\t3\nEuropean_breed\t4\nChina_breed\t5\nKorean\t6\nKorean_black\t7\nWild_Boars_Asia\t8\nWild_Boars_European\t9" > breed.order

awk '{print $3"\t"$0}' final.id.info | sort -k1,1 | join -t $'\t' - <(sort -k1,1 breed.order) | sort -k3,3n | paste - final.id.fam | sort -k5,5n | awk -F "\t" '{print $6}' > final.id.ordered.fam

shuf final.id.ordered.fam > final.id.random.fam

for i in `seq 1 1 19` `seq 20 5 100` `seq 200 200 1000` 1500 2000 2259
do
  head -n $i final.id.random.fam > n"$i".random.id
done

for i in `seq 1 1 2259`
do
  head -n $i final.id.ordered.fam > n"$i".id
done

# 3. do jobs by chromosome
# ============================================================

for i in `seq 1 1 18`
do
	sbatch --export=chr="chr$i" --output=chr"$i".countSNPs.out --error=chr"$i".count.SNPs.err countSNPs.sbatch
done

# 4. summarize variants
# ============================================================

cat chr*.snp.count | sort -k1,1 | /mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2,3 -o sum,sum | sort -k1,1n > all.ids.snp.count

# 5. random
# ============================================================

rm chr*.snp.count
for i in `seq 1 1 18`
do
	sbatch --export=chr="chr$i" --output=chr"$i".countSNPs.out --error=chr"$i".count.SNPs.err countSNPs2.sbatch
done
cat chr*.snp.count | sort -k1,1 | /mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2 -o sum | sort -k1,1n > all.ids.random.snp.count



