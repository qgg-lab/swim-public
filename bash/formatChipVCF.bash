# Illumina 50K
# ============================================================

echo -e '##fileformat=VCFv4.1\n##fileDate=09120221\n##INFO=<ID=CHIP,Number=1,Type=String,Description="Chip type."\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' > PorcineSNP50K.vcf
tail -n+8 PorcineSNP50K_45164.gold.vcf | sed 's/ /\t/g' | cut -f 1-7 | awk '{print $0"\tCHIP=PorcineSNP50K"}' >> PorcineSNP50K.vcf

java -jar /mnt/research/qgg/software/picard-2.23.3/picard.jar SortVcf I=PorcineSNP50K.vcf O=PorcineSNP50K.sorted.vcf SD=genome.dict CREATE_INDEX=true > /mnt/research/qgg/resource/swim/PorcineSNP50K.sorted.log 2>&1 &

# Illumina 60K
# ============================================================

echo -e '##fileformat=VCFv4.1\n##fileDate=09120221\n##INFO=<ID=CHIP,Number=1,Type=String,Description="Chip type."\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' > PorcineSNP60K.vcf
tail -n+8 PorcineSNP60K_52866.gold.vcf | sed 's/ /\t/g' | cut -f 1-7 | awk '{print $0"\tCHIP=PorcineSNP60K"}' >> PorcineSNP60K.vcf

java -jar /mnt/research/qgg/software/picard-2.23.3/picard.jar SortVcf I=PorcineSNP60K.vcf O=PorcineSNP60K.sorted.vcf SD=genome.dict CREATE_INDEX=true > /mnt/research/qgg/resource/swim/PorcineSNP60K.sorted.log 2>&1 &

# Illumina 80K
# ============================================================

echo -e '##fileformat=VCFv4.1\n##fileDate=09120221\n##INFO=<ID=CHIP,Number=1,Type=String,Description="Chip type."\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' > PorcineSNP80K.vcf
tail -n+8 Affymetrix_PigSNP130K.111522.gold.vcf | sed 's/ /\t/g' | cut -f 1-7 | awk '{print $0"\tCHIP=PorcineSNP80K"}' >> PorcineSNP80K.vcf

java -jar /mnt/research/qgg/software/picard-2.23.3/picard.jar SortVcf I=PorcineSNP80K.vcf O=PorcineSNP80K.sorted.vcf SD=genome.dict CREATE_INDEX=true > /mnt/research/qgg/resource/swim/PorcineSNP80K.sorted.log 2>&1 &

# affymetrix sowpro90
# ============================================================

echo -e '##fileformat=VCFv4.1\n##fileDate=09120221\n##INFO=<ID=CHIP,Number=1,Type=String,Description="Chip type."\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' > AffySP90.vcf
tail -n+9 Affymetrix_PigSNP130K.111522.gold.vcf | sed 's/ /\t/g' | cut -f 1-7 | awk '$5 != "-" {print $0"\tCHIP=AffySNP90"}' >> AffySP90.vcf

java -jar /mnt/research/qgg/software/picard-2.23.3/picard.jar SortVcf I=AffySP90.vcf O=AffySP90.sorted.vcf SD=genome.dict CREATE_INDEX=true > /mnt/research/qgg/resource/swim/AffySP90.sorted.log 2>&1 &

# affymetrix 660k
# ============================================================

echo -e '##fileformat=VCFv4.1\n##fileDate=09120221\n##INFO=<ID=CHIP,Number=1,Type=String,Description="Chip type."\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' > Affy660K.vcf
tail -n+9 Axiom_PigSNP660K.623618.gold.vcf | sed 's/ /\t/g' | cut -f 1-7 | awk '$5 != "-" {print $0"\tCHIP=Affy660K"}' >> Affy660K.vcf

java -jar /mnt/research/qgg/software/picard-2.23.3/picard.jar SortVcf I=Affy660K.vcf O=Affy660K.sorted.vcf SD=genome.dict CREATE_INDEX=true > /mnt/research/qgg/resource/swim/Affy660K.sorted.log 2>&1 &


