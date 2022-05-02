####SNP chip pre-update
####Since the genotypes of SNPs may be displayed in different forms,
####we need to generate a script that automatically converts the genotypes of SNPs.

#### 50k
#### delete [A/T] and [C/G] SNP
awk '{if (($4== "A" && $5!= "T") || ($4== "T" && $5!= "A") || ($4== "G" && $5!= "C") || ($4== "C" && $5!= "G") ) print $0}' \
PorcineSNP50K_45164.info > PorcineSNP50K_45105_final.info 

#### Extract valid information for subsequent update of SNP and allele information
awk '{print $3}' PorcineSNP50K_45105_final.info  >50k_final_SNP
awk '{print $3,$2}' PorcineSNP50K_45105_final.info  >50k_final_SNPpos
awk '{print $3,$1}' PorcineSNP50K_45105_final.info  >50k_final_SNPchr
awk '{print $3,$4}' PorcineSNP50K_45105_final.info  >50k_final_SNPref

#### 60k
#### delete [A/T] and [C/G] SNP
awk '{if (($4== "A" && $5!= "T") || ($4== "T" && $5!= "A") || ($4== "G" && $5!= "C") || ($4== "C" && $5!= "G") ) print $0}' \
PorcineSNP60K_52866.info > PorcineSNP60K_52693_final.info 

#### Extract valid information for subsequent update of SNP and allele information
awk '{print $3}' PorcineSNP60K_52693_final.info  >60k_final_SNP
awk '{print $3,$2}' PorcineSNP60K_52693_final.info  >60k_final_SNPpos
awk '{print $3,$1}' PorcineSNP60K_52693_final.info  >60k_final_SNPchr
awk '{print $3,$4}' PorcineSNP60K_52693_final.info  >60k_final_SNPref

#### 80k
#### delete [A/T] and [C/G] SNP
awk '{if (($4== "A" && $5!= "T") || ($4== "T" && $5!= "A") || ($4== "G" && $5!= "C") || ($4== "C" && $5!= "G") ) print $0}' \
PorcineSNP80K_60818.info > PorcineSNP80K_60737_final.info 

#### Extract valid information for subsequent update of SNP and allele information
awk '{print $3}' PorcineSNP80K_60737_final.info  >80k_final_SNP
awk '{print $3,$2}' PorcineSNP80K_60737_final.info  >80k_final_SNPpos
awk '{print $3,$1}' PorcineSNP80K_60737_final.info  >80k_final_SNPchr
awk '{print $3,$4}' PorcineSNP80K_60737_final.info  >80k_final_SNPref

#### get  UpdateAlleles file
####test in 50k chip
####R

testname=test_100

module load R

Rscript snpchip_allele_renew.r "PorcineSNP50K_45105_final.info" "~/duroc_SNP_CHIP/test_100.bim" ${testname}


####plink
####test

####extract SNP and update position
/mnt/research/qgg/software/plink-v1.90b5.3/plink \
--bfile ${testname} \
--extract  50k_final_SNP \
--update-map 50k_final_SNPpos \
--make-bed --out ${testname}_snpPos

####update chromosome
/mnt/research/qgg/software/plink-v1.90b5.3/plink \
--bfile ${testname}_snpPos \
--update-map 50k_final_SNPchr \
--update-chr \
--make-bed --out ${testname}_snpPosChr

####update Alleles
/mnt/research/qgg/software/plink-v1.90b5.3/plink \
--bfile ${testname}_snpPosChr \
--update-alleles test_100_UpdateAlleles.txt \
--make-bed --out ${testname}_snpPosChrAllele

####update reference Alleles
/mnt/research/qgg/software/plink-v1.90b5.3/plink \
--bfile ${testname}_snpPosChrAllele \
--reference-allele 50k_final_SNPref \
--make-bed --out ${testname}_snpPosChrAlleleRef

#### get final unphase SNPchip vcf 
/mnt/research/qgg/software/plink-v1.90b5.3/plink \
--bfile ${testname}_snpPosChrAlleleRef \
--recode vcf-iid --out ${testname}_regeno
