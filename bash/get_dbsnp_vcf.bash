#get vcffile head
grep  ^# /mnt/research/qgg/resource/swim/dbsnp.nospace.vcf >dbsnp.reliable.snp.vcf
grep  ^# /mnt/research/qgg/resource/swim/dbsnp.nospace.vcf >dbsnp.other.snp.vcf

####filter reliable snp
grep -v ^# /mnt/research/qgg/resource/swim/dbsnp.nospace.vcf | awk -F ";" '{if( length($3)>0 && $2=="TSA=SNV" ) print $0}' >>dbsnp.reliable.snp.vcf
####filter other snp
grep -v ^# /mnt/research/qgg/resource/swim/dbsnp.nospace.vcf | awk -F ";" '{if( length($3)==0 && $2=="TSA=SNV" ) print $0}' >>dbsnp.other.snp.vcf

####keep only InDels
/mnt/research/qgg/software/vcftools_0.1.13/bin/vcftools \
--vcf /mnt/research/qgg/resource/swim/dbsnp.nospace.vcf \
--keep-only-indels \
--recode --recode-INFO-all \
--out dbsnp.onlyInDel.vcf
