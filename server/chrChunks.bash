# =========================================
# = make chromosome chunks for imputation =
# =========================================

for i in `seq 1 1 18`
do
	gunzip -c ~/scratch/share/imputationQCvcf/chr"$i".QC.final.shapeit4.phased.vcf.gz | tail -n+11 | awk '{print $1"\t"$2-1"\t"$2-1+length($4)}' > chr"$i".snp.bed
	awk '$1 == '$i' {print $1"\t0\t"$2}' ~/qgg/resource/sscrofa11.1/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.fai | ~/qgg/software/bedtools-2.29.2/bin/bedtools subtract -a - -b chr"$i".snp.bed | awk '$3-$2 >= 10' >> all.sub.snp.bed
	echo "$i"
done &

sort -k1,1 -k2,2n all.sub.snp.bed > all.sub.snp.sorted.bed

~/qgg/software/bedtools-2.29.2/bin/bedtools makewindows -g ~/qgg/resource/sscrofa11.1/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.fai -w 20000000 | awk '$1 >= 1 && $1 <= 18 {print $1"\t"$3-1"\t"$3}'| sort -k1,1 -k2,2n | ~/qgg/software/bedtools-2.29.2/bin/bedtools closest -a - -b all.sub.snp.sorted.bed > all.chr.window.int.bed

# finally, get windows
# ============================================================

cat <(echo -e "0\t0\t0") <(cut -f 4-6 all.chr.window.int.bed) | paste - <(cut -f 4-6 all.chr.window.int.bed) | awk 'NF > 3 { if ($1 != $4) { print $4"\t0\t"$5+5 } else { print $4"\t"$2+5"\t"$5+5 }}' | 