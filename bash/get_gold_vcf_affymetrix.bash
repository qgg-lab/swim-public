#####################################RUN in MSU######################################
##利用芯片公司提供的SNP的Name及SNP的AlleleA_ProbeSeq筛选
##pig_660k_AlleleA_ProbeSeq,pig_660k_AlleleB_ProbeSeq.txt文件
##pig_130k_AlleleA_ProbeSeq,pig_130k_AlleleB_ProbeSeq.txt文件
##提前根据参考基因型情况情况生成info文件，660K_info.txt和130K_info.txt
##提前生成vcf的抬头信息，com_head
##具体用excel完成
########################################################################################

###################################
###1、利用probe文件生成fastq文件
###################################

##确保第2列的长度为35
awk '{if(length($2)!=35) print $0}' pig_660k_AlleleA_ProbeSeq.txt
awk '{if(length($2)!=35) print $0}' pig_660k_AlleleB_ProbeSeq.txt

awk '{if(length($2)!=35) print $0}' pig_130k_AlleleA_ProbeSeq.txt
awk '{if(length($2)!=35) print $0}' pig_130k_AlleleB_ProbeSeq.txt

##先生成fastq文件
awk '{if(length($2)==35) print "@"$1"\n"$2"\n""+""\n""IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"}' \
pig_660k_AlleleA_ProbeSeq.txt >pig_660k_AlleleA_ProbeSeq.fq

awk '{if(length($2)==35) print "@"$1"\n"$2"\n""+""\n""IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"}' \
pig_660k_AlleleB_ProbeSeq.txt >pig_660k_AlleleB_ProbeSeq.fq

awk '{if(length($2)==35) print "@"$1"\n"$2"\n""+""\n""IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"}' \
pig_130k_AlleleA_ProbeSeq.txt >pig_130k_AlleleA_ProbeSeq.fq

awk '{if(length($2)==35) print "@"$1"\n"$2"\n""+""\n""IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"}' \
pig_130k_AlleleB_ProbeSeq.txt >pig_130k_AlleleB_ProbeSeq.fq

####################################
#2、利用bwa比对fastq结果
##注意bwa在mapping之前需要先index
####################################
####660K
/mnt/research/qgg/software/bwa-0.7.17/bwa aln -t 5 \
/mnt/research/qgg/resource/swim/bwa0717 pig_660k_AlleleA_ProbeSeq.fq  \
>pig_660k_AlleleA_ProbeSeq_aln.sai 2> pig_660k_AlleleA_ProbeSeq_aln_bwa.log

/mnt/research/qgg/software/bwa-0.7.17/bwa samse \
/mnt/research/qgg/resource/swim/bwa0717 pig_660k_AlleleA_ProbeSeq_aln.sai  pig_660k_AlleleA_ProbeSeq.fq\
>pig_660k_AlleleA_ProbeSeq_aln.sam 2> pig_660k_AlleleA_ProbeSeq_samse_bwa.log

/mnt/research/qgg/software/bwa-0.7.17/bwa aln -t 5 \
/mnt/research/qgg/resource/swim/bwa0717 pig_660k_AlleleB_ProbeSeq.fq  \
>pig_660k_AlleleB_ProbeSeq_aln.sai 2> pig_660k_AlleleB_ProbeSeq_aln_bwa.log

/mnt/research/qgg/software/bwa-0.7.17/bwa samse \
/mnt/research/qgg/resource/swim/bwa0717 pig_660k_AlleleB_ProbeSeq_aln.sai  pig_660k_AlleleB_ProbeSeq.fq\
>pig_660k_AlleleB_ProbeSeq_aln.sam 2> pig_660k_AlleleB_ProbeSeq_samse_bwa.log


grep -v ^@  pig_660k_AlleleA_ProbeSeq_aln.sam | \
grep -w "MD:Z:35" | \
awk '{if ($12 == "XT:A:U") print $0}' | \
awk 'BEGIN{print "Name Strand chr pos"} { if ($2 == 0) print $1,"+",$3,$4+35; else if($2 == 16) print $1,"-",$3,$4-1}' > pig_660k_A584400_info1.txt

grep -v ^@  pig_660k_AlleleB_ProbeSeq_aln.sam | \
grep -w "MD:Z:35" | \
awk '{if ($12 == "XT:A:U") print $0}' | \
awk 'BEGIN{print "Name Strand chr pos"} { if ($2 == 0) print $1,"+",$3,$4-1; else if($2 == 16) print $1,"-",$3,$4+35}'  > pig_660k_B44766_info1.txt

####130K
/mnt/research/qgg/software/bwa-0.7.17/bwa aln -t 5 \
/mnt/research/qgg/resource/swim/bwa0717 pig_130k_AlleleA_ProbeSeq.fq  \
>pig_130k_AlleleA_ProbeSeq_aln.sai 2> pig_130k_AlleleA_ProbeSeq_aln_bwa.log

/mnt/research/qgg/software/bwa-0.7.17/bwa samse \
/mnt/research/qgg/resource/swim/bwa0717 pig_130k_AlleleA_ProbeSeq_aln.sai  pig_130k_AlleleA_ProbeSeq.fq\
>pig_130k_AlleleA_ProbeSeq_aln.sam 2> pig_130k_AlleleA_ProbeSeq_samse_bwa.log

/mnt/research/qgg/software/bwa-0.7.17/bwa aln -t 5 \
/mnt/research/qgg/resource/swim/bwa0717 pig_130k_AlleleB_ProbeSeq.fq  \
>pig_130k_AlleleB_ProbeSeq_aln.sai 2> pig_130k_AlleleB_ProbeSeq_aln_bwa.log

/mnt/research/qgg/software/bwa-0.7.17/bwa samse \
/mnt/research/qgg/resource/swim/bwa0717 pig_130k_AlleleB_ProbeSeq_aln.sai  pig_130k_AlleleB_ProbeSeq.fq\
>pig_130k_AlleleB_ProbeSeq_aln.sam 2> pig_130k_AlleleB_ProbeSeq_samse_bwa.log

grep -v ^@  pig_130k_AlleleA_ProbeSeq_aln.sam | \
grep -w "MD:Z:35" | \
awk '{if ($12 == "XT:A:U") print $0}' | \
awk 'BEGIN{print "Name Strand chr pos"} { if ($2 == 0) print $1,"+",$3,$4+35; else if($2 == 16) print $1,"-",$3,$4-1}' > pig_130k_A111753_info1.txt

grep -v ^@  pig_130k_AlleleB_ProbeSeq_aln.sam | \
grep -w "MD:Z:35" | \
awk '{if ($12 == "XT:A:U") print $0}' | \
awk 'BEGIN{print "Name Strand chr pos"} { if ($2 == 0) print $1,"+",$3,$4-1; else if($2 == 16) print $1,"-",$3,$4+35}'  > pig_130k_B5401_info1.txt



####提取660k的参考基因型，即REF
awk -v OFS='\t' ' {if (NR>1) print $3,$4-1,$4,$1}' pig_660k_A584400_info1.txt >660K_pos_A.bed

awk -v OFS='\t' ' {if (NR>1) print $3,$4-1,$4,$1}' pig_660k_B44766_info1.txt >660K_pos_B.bed

/mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools getfasta -fi /mnt/research/qgg/resource/swim/genome.fa -bed 660K_pos_A.bed -nameOnly -fo 660_refA
/mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools getfasta -fi /mnt/research/qgg/resource/swim/genome.fa -bed 660K_pos_B.bed -nameOnly -fo 660_refB

awk 'BEGIN {print "REF"} {if (NR%2==0) print $0}' 660_refA | paste pig_660k_A584400_info1.txt - >pig_660k_A584400_info2.txt
awk 'BEGIN {print "REF"} {if (NR%2==0) print $0}' 660_refB | paste pig_660k_B44766_info1.txt - >pig_660k_B44766_info2.txt

####提取130k的参考基因型，即REF
awk -v OFS='\t' ' {if (NR>1) print $3,$4-1,$4,$1}' pig_130k_A111753_info1.txt >130K_pos_A.bed

awk -v OFS='\t' ' {if (NR>1) print $3,$4-1,$4,$1}' pig_130k_B5401_info1.txt >130K_pos_B.bed

/mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools getfasta -fi /mnt/research/qgg/resource/swim/genome.fa -bed 130K_pos_A.bed -nameOnly -fo 130_refA
/mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools getfasta -fi /mnt/research/qgg/resource/swim/genome.fa -bed 130K_pos_B.bed -nameOnly -fo 130_refB

awk 'BEGIN {print "REF"} {if (NR%2==0) print $0}' 130_refA | paste pig_130k_A111753_info1.txt - >pig_130k_A111753_info2.txt
awk 'BEGIN {print "REF"} {if (NR%2==0) print $0}' 130_refB | paste pig_130k_B5401_info1.txt - >pig_130k_B5401_info2.txt


####switch Ref/Alt and get final vcf file
####run in R program for 660k
affymetrix660K_switch_RefAlt.r

cat 660k_head Axiom_PigSNP660K_623638.info > Axiom_PigSNP660K.623618.gold.vcf

####run in R program for 130k
affymetrix130K_switch_RefAlt.r

cat 130k_head Affymetrix_PigSNP130K_111522.info > Affymetrix_PigSNP130K.111522.gold.vcf

