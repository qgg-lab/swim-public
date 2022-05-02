# ==============================
# = prepare resources for swim =
# ==============================

mkdir /mnt/research/qgg/swim/

# 1. genome sequence
# ============================================================

wget -O /mnt/research/qgg/resource/swim/genome.fa.gz ftp://ftp.ensembl.org/pub/release-100/fasta/sus_scrofa/dna//Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.gz
gunzip /mnt/research/qgg/resource/swim/genome.fa.gz
java -jar /mnt/research/qgg/software/picard-2.23.3/picard.jar CreateSequenceDictionary R=/mnt/research/qgg/resource/swim/genome.fa O=/mnt/research/qgg/resource/swim/genome.dict > /mnt/research/qgg/resource/swim/dict.log 2>&1 &

# 2. dbSNP
# ============================================================

wget -O /mnt/research/qgg/resource/swim/dbsnp.vcf.gz ftp://ftp.ensembl.org/pub/release-100/variation/vcf/sus_scrofa//sus_scrofa.vcf.gz
gunzip /mnt/research/qgg/resource/swim/dbsnp.vcf.gz
# the dbsnp vcf file contains axiom snp array snps that offends the gatk
# because they have white space in the INFO field, remove these (about 13k)
awk -F "\t" '!($8 ~ /\s/)' /mnt/research/qgg/resource/swim/dbsnp.vcf > /mnt/research/qgg/resource/swim/dbsnp.nospace.vcf
java -jar /mnt/research/qgg/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar IndexFeatureFile --input /mnt/research/qgg/resource/swim/dbsnp.nospace.vcf > /mnt/research/qgg/resource/swim/dbsnp.nospace.index.log 2>&1 &
# split dbsnp into SNP and indel
cat /mnt/research/qgg/resource/swim/dbsnp.nospace.vcf | perl -wne 'chomp $_; if ($_ =~ /^\#/) { print $_, "\n"; } else { @line = split /\t/, $_; $line[7] .= ";"; if ($line[7] =~ m/TSA=(.*?);/) { if ($1 eq "SNV") { print $_, "\n"; } } }' > /mnt/research/qgg/resource/swim/dbsnp.snv.vcf &
cat /mnt/research/qgg/resource/swim/dbsnp.nospace.vcf | perl -wne 'chomp $_; if ($_ =~ /^\#/) { print $_, "\n"; } else { @line = split /\t/, $_; $line[7] .= ";"; if ($line[7] =~ m/TSA=(.*?);/) { if ($1 eq "insertion" || $1 eq "deletion") { print $_, "\n"; } } }' > /mnt/research/qgg/resource/swim/dbsnp.indel.vcf &
java -jar /mnt/research/qgg/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar IndexFeatureFile --input /mnt/research/qgg/resource/swim/dbsnp.snv.vcf > /mnt/research/qgg/resource/swim/dbsnp.snv.index.log 2>&1 &
java -jar /mnt/research/qgg/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar IndexFeatureFile --input /mnt/research/qgg/resource/swim/dbsnp.indel.vcf > /mnt/research/qgg/resource/swim/dbsnp.indel.index.log 2>&1 &

# 3. index fa
# ============================================================

~/qgg/software/samtools-1.10/samtools faidx /mnt/research/qgg/resource/swim/genome.fa > /mnt/research/qgg/resource/swim/samtools.faidx.log 2>&1 &

# 4. index for bwa
# ============================================================

srun --mem=4G --nodes=1 --ntasks-per-node=1 --cpus-per-task=8 --time=4:00:00 ~/qgg/software/bwa-0.7.17/bwa index -p /mnt/research/qgg/resource/swim/bwa0717 /mnt/research/qgg/resource/swim/genome.fa > /mnt/research/qgg/resource/swim/bwa0717.index.log 2>&1 &

# 5. make interval lists split the genome to approximately
#    600 Mb
# ============================================================

awk '$1 >= 1 && $1 <= 4 {print $1}' /mnt/research/qgg/resource/swim/genome.fa.fai > /mnt/research/qgg/resource/swim/chr600m1.list
awk '$1 >= 5 && $1 <= 9 {print $1}' /mnt/research/qgg/resource/swim/genome.fa.fai > /mnt/research/qgg/resource/swim/chr600m2.list
awk '$1 >= 10 && $1 <= 14 {print $1}' /mnt/research/qgg/resource/swim/genome.fa.fai > /mnt/research/qgg/resource/swim/chr600m3.list
awk '!($1 <= 14) {print $1}' /mnt/research/qgg/resource/swim/genome.fa.fai > /mnt/research/qgg/resource/swim/chr600m4.list
awk '$1 == "X" || $1 == "Y" || $1 == "MT" {print $1}' /mnt/research/qgg/resource/swim/genome.fa.fai > /mnt/research/qgg/resource/swim/chr600m5.list

