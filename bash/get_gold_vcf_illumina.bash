#####################################RUN in MSU######################################
##利用芯片公司提供的SNP的Name及SNP的AlleleA_ProbeSeq筛选
##pig_50k_AlleleA_ProbeSeq.txt,pig_60k_AlleleA_ProbeSeq.txt和pig_80k_AlleleA_ProbeSeq.txt文件
##提前根据链的TOP和BOT及基因型情况情况生成info文件，50K_info.txt，60K_info.txt，80K_info.txt
##提前生成vcf的抬头信息，com_head
##具体用excel完成
########################################################################################

###################################
###1、利用probe文件生成fastq文件
###################################

cd SNP_pos_renew

##确保第2列的长度为50
awk '{if(length($2)!=50) print $0}' pig_50k_AlleleA_ProbeSeq.txt

awk '{if(length($2)!=50) print $0}' pig_60k_AlleleA_ProbeSeq.txt

awk '{if(length($2)!=50) print $0}' pig_80k_AlleleA_ProbeSeq.txt

##先生成fastq文件
awk '{if(length($2)==50) print "@"$1"\n"$2"\n""+""\n""IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"}' \
pig_50k_AlleleA_ProbeSeq.txt >pig_50k_AlleleA_ProbeSeq.fq

awk '{if(length($2)==50) print "@"$1"\n"$2"\n""+""\n""IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"}' \
pig_60k_AlleleA_ProbeSeq.txt >pig_60k_AlleleA_ProbeSeq.fq

awk '{if(length($2)==50) print "@"$1"\n"$2"\n""+""\n""IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"}' \
pig_80k_AlleleA_ProbeSeq.txt >pig_80k_AlleleA_ProbeSeq.fq


####################################
#2、利用bwa比对fastq结果
##注意bwa在mapping之前需要先index
####################################
####50K
/mnt/research/qgg/software/bwa-0.7.17/bwa aln -t 5 \
/mnt/research/qgg/resource/swim/bwa0717 pig_50k_AlleleA_ProbeSeq.fq  \
>pig_50k_AlleleA_ProbeSeq_aln.sai 2> pig_50k_AlleleA_ProbeSeq_aln_bwa.log

/mnt/research/qgg/software/bwa-0.7.17/bwa samse \
/mnt/research/qgg/resource/swim/bwa0717 pig_50k_AlleleA_ProbeSeq_aln.sai  pig_50k_AlleleA_ProbeSeq.fq\
>pig_50k_AlleleA_ProbeSeq_aln.sam 2> pig_50k_AlleleA_ProbeSeq_samse_bwa.log

grep -v ^@  pig_50k_AlleleA_ProbeSeq_aln.sam | \
grep -w "MD:Z:50" | \
awk '{if ($12 == "XT:A:U") print $0}' | \
awk 'BEGIN{print "Name Strand chr pos"} { if ($2 == 0) print $1,"+",$3,$4+50; else if($2 == 16) print $1,"-",$3,$4-1}' > pig_50k_45579_info1.txt


awk -v OFS='\t' ' {if (NR>1) print $3,$4-1,$4,$1}' pig_50k_45579_info1.txt >50K_pos.bed

/mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools getfasta -fi /mnt/research/qgg/resource/swim/genome.fa -bed 50K_pos.bed -nameOnly -fo 50K_ref

awk 'BEGIN {print "REF"} {if (NR%2==0) print $0}' 50K_ref | paste pig_50k_45579_info1.txt - >pig_50k_45579_info2.txt

####60K
/mnt/research/qgg/software/bwa-0.7.17/bwa aln -t 5 \
/mnt/research/qgg/resource/swim/bwa0717 pig_60k_AlleleA_ProbeSeq.fq  \
>pig_60k_AlleleA_ProbeSeq_aln.sai 2> pig_60k_AlleleA_ProbeSeq_aln_bwa.log

/mnt/research/qgg/software/bwa-0.7.17/bwa samse \
/mnt/research/qgg/resource/swim/bwa0717 pig_60k_AlleleA_ProbeSeq_aln.sai  pig_60k_AlleleA_ProbeSeq.fq\
>pig_60k_AlleleA_ProbeSeq_aln.sam 2> pig_60k_AlleleA_ProbeSeq_samse_bwa.log

grep -v ^@  pig_60k_AlleleA_ProbeSeq_aln.sam | \
grep -w "MD:Z:50" | \
awk '{if ($12 == "XT:A:U") print $0}' | \
awk 'BEGIN{print "Name Strand chr pos"} { if ($2 == 0) print $1,"+",$3,$4+50; else if($2 == 16) print $1,"-",$3,$4-1}'> pig_60k_53408_info1.txt

awk -v OFS='\t' ' {if (NR>1) print $3,$4-1,$4,$1}' pig_60k_53408_info1.txt >60K_pos.bed

/mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools getfasta -fi /mnt/research/qgg/resource/swim/genome.fa -bed 60K_pos.bed -nameOnly -fo 60K_ref

awk 'BEGIN {print "REF"} {if (NR%2==0) print $0}' 60K_ref | paste pig_60k_53408_info1.txt - >pig_60k_53408_info2.txt

####80K
/mnt/research/qgg/software/bwa-0.7.17/bwa aln -t 5 \
/mnt/research/qgg/resource/swim/bwa0717 pig_80k_AlleleA_ProbeSeq.fq  \
>pig_80k_AlleleA_ProbeSeq_aln.sai 2> pig_80k_AlleleA_ProbeSeq_aln_bwa.log

/mnt/research/qgg/software/bwa-0.7.17/bwa samse \
/mnt/research/qgg/resource/swim/bwa0717 pig_80k_AlleleA_ProbeSeq_aln.sai  pig_80k_AlleleA_ProbeSeq.fq\
>pig_80k_AlleleA_ProbeSeq_aln.sam 2> pig_80k_AlleleA_ProbeSeq_samse_bwa.log

grep -v ^@  pig_80k_AlleleA_ProbeSeq_aln.sam | \
grep -w "MD:Z:50" | \
awk '{if ($12 == "XT:A:U") print $0}' | \
awk 'BEGIN{print "Name Strand chr pos"} { if ($2 == 0) print $1,"+",$3,$4+50; else if($2 == 16) print $1,"-",$3,$4-1}' > pig_80k_61371_info1.txt

awk -v OFS='\t' ' {if (NR>1) print $3,$4-1,$4,$1}' pig_80k_61371_info1.txt >80K_pos.bed

/mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools getfasta -fi /mnt/research/qgg/resource/swim/genome.fa -bed 80K_pos.bed -nameOnly -fo 80K_ref

awk 'BEGIN {print "REF"} {if (NR%2==0) print $0}' 80K_ref | paste pig_80k_61371_info1.txt - >pig_80k_61371_info2.txt


####switch Ref/Alt and get final vcf file
####run in R program for 50k
illumina50K_switch_RefAlt.r

cat 50k_head PorcineSNP50K_45164.info>PorcineSNP50K_45164.gold.vcf

####run in R program for 60k
illumina60K_switch_RefAlt.r

cat 60k_head PorcineSNP60K_52866.info>PorcineSNP60K_52866.gold.vcf

####run in R program for 80k
illumina80K_switch_RefAlt.r

cat 80k_head PorcineSNP80K_60818.info>PorcineSNP80K_60818.gold.vcf


####combine 50k 60k 80k 130k snpchip info
####run in R
combine_50k_60k_80k_130k_info.r

cat com_head PorcineSNPcomb_138162.info >combinePorcineSNPCHIP.138162.gold.vcf

