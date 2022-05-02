mkdir /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550
mkdir /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/log
mkdir /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/tmp
cd /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550
cp /mnt/ls15/scratch/users/dingrong/imputation_test/resource.env /mnt/ls15/scratch/users/dingrong/swim/

####extract target LR100

awk '{if ($16 == "LAB" && $15 == "LR" ) print $0 }' /mnt/ls15/scratch/users/dingrong/swim/SWIM_sampleInfo_2259.txt \
|sort -k18n |head -100 | awk '{ print $1 }' >/mnt/ls15/scratch/users/dingrong/swim/imputationEligible/target.LR100.ind


for chr in $(seq 1 18) 
do
  /mnt/research/qgg/software/vcftools-v0.1.16/bin/vcftools \
  --vcf /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/chr"$chr".QC.final.vcf \
  --keep /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/target.LR100.ind \
  --recode \
  --out /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/target.LR100.chr"$chr" \
  > /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/target.LR100.chr"$chr".log 2>&1 &
done

for chr in $(seq 1 18) 
do
#  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=64G --time=1:50:00 \
  /mnt/research/qgg/software/vcftools-v0.1.16/bin/vcftools \
  --vcf /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/target.LR100.chr"$chr".recode.vcf \
  --snps /mnt/ls15/scratch/users/dingrong/imputation_test/PorcineSNP50K.list \
  --remove-indels \
  --recode \
  --out /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/target.chip.chr"$chr" \
  > /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/target.chip.chr"$chr".log 2>&1 &
done

for chr in $(seq 1 18)
do
  /mnt/research/qgg/software/plink-v1.90b6.18/plink \
  --vcf target.chip.chr"$chr".recode.vcf \
  --freq \
  --out miss_summary/target.chip.chr"$chr" \
  > log/target.chip.chr"$chr".freq.log 2>&1 &  
done


####
for chr in $(seq 1 18) 
do
echo $chr
bcftools  view target.chip.chr"$chr".recode.vcf -Oz -o target.chip.chr"$chr".vcf.gz 
bcftools index -f -t --threads 12 target.chip.chr"$chr".vcf.gz
done


#random select 550 ref individuals

grep -v -f /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/target.LR100.ind /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/Landrace.651.ind | \
shuf |awk '{print $0}'  |head -550 > /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/ref1.LR550.ind

shuf /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/ref1.LR550.ind | awk '{print $0}'  |head -250 > LR.select.250.ind

shuf /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/Duroc.485.ind | awk '{print $0}'  |head -150 \
> DUC.select.150.ind

shuf /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/Large_White.543.ind | awk '{print $0}'  |head -150 \
> LW.select.150.ind

cat LR.select.250.ind DUC.select.150.ind LW.select.150.ind >ref2.DLY550.ind

awk '{if (NR >1 && $12!="European_breed") print $1 }' /mnt/ls15/scratch/users/dingrong/swim/SWIM_sampleInfo_2259.txt \
 |shuf |awk '{print $0}'  |head -300 > noDLY.select.300.ind

cat LR.select.250.ind noDLY.select.300.ind >ref3.noDLY550.ind

# select 
for chr in $(seq 1 18) 
do
  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=64G --time=1:59:00 \
  /mnt/research/qgg/software/vcftools-v0.1.16/bin/vcftools \
  --vcf /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/chr"$chr".QC.final.vcf \
  --keep /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/ref1.LR550.ind \
  --recode \
  --out /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/ref1.LR550.chr"$chr" \
  > log/ref1.LR550.chr"$chr".filter.log 2>&1 &

  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=64G --time=1:59:00 \
  /mnt/research/qgg/software/vcftools-v0.1.16/bin/vcftools \
  --vcf /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/chr"$chr".QC.final.vcf \
  --keep /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/ref2.DLY550.ind \
  --recode \
  --out /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/ref2.DLY550.chr"$chr" \
  > log/ref2.DLY550.chr"$chr".filter.log 2>&1 &

  srun --cpus-per-task=1 --ntasks-per-node=1 --mem=64G --time=5:00:00 \
  /mnt/research/qgg/software/vcftools-v0.1.16/bin/vcftools \
  --vcf /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/chr"$chr".QC.final.vcf \
  --keep /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/ref3.noDLY550.ind \
  --recode \
  --out /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/ref3.noDLY550.chr"$chr" \
  > log/ref3.noDLY550.chr"$chr".filter.log 2>&1 &  
   
done

cd /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/

for chr in $(seq 1 18) 
do
echo $chr
bcftools view ref1.LR550.chr"$chr".recode.vcf -Oz -o ref1.LR550.chr"$chr".vcf.gz 
bcftools index -f -t --threads 24 ref1.LR550.chr"$chr".vcf.gz
done


for chr in $(seq 1 18) 
do
echo $chr
bcftools view ref2.DLY550.chr"$chr".recode.vcf -Oz -o ref2.DLY550.chr"$chr".vcf.gz 
bcftools index -f -t --threads 24 ref2.DLY550.chr"$chr".vcf.gz
done


for chr in $(seq 1 18) 
do
echo $chr	
bcftools  view ref3.noDLY550.chr"$chr".recode.vcf -Oz -o ref3.noDLY550.chr"$chr".vcf.gz 
bcftools index -f -t --threads 24 ref3.noDLY550.chr"$chr".vcf.gz
done

###phasing
filepath=/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550


for ref in ref1.LR550 ref2.DLY550 ref3.noDLY550
do
for chr in $(seq 1 18) 
do
###eagle_ref
sbatch --cpus-per-task=8 \
--mem=128G \
--time=23:00:00 \
--export=env=/mnt/ls15/scratch/users/dingrong/swim/resource.env,filepath=$filepath,input="$ref".chr"$chr".vcf.gz,output="$ref".chr"$chr".eagle.phased \
--output=log/"$ref".chr"$chr".eagle.phasing.out \
--error=log/"$ref".chr"$chr".eagle.phasing.err \
/mnt/ls15/scratch/users/dingrong/swim/phasing_eagle_ref.sbatch

###beagle_ref
sbatch --cpus-per-task=12 \
--time=6:00:00 \
--mem=64G \
--export=env=/mnt/ls15/scratch/users/dingrong/swim/resource.env,filepath=$filepath,input="$ref".chr"$chr".vcf.gz,output="$ref".chr"$chr".beagle.phased,tmp=$filepath/tmp,mem=64G \
--output=log/"$ref".chr"$chr".beagle.phasing.out \
--error=log/"$ref".chr"$chr".beagle.phasing.err \
/mnt/ls15/scratch/users/dingrong/swim/phasing_beagle.sbatch

###shapeit4_ref
sbatch -C intel16 \
--cpus-per-task=12 \
--time=8:50:00 \
--mem=64G \
--export=env=/mnt/ls15/scratch/users/dingrong/swim/resource.env,filepath=$filepath,input="$ref".chr"$chr".vcf.gz,output="$ref".chr"$chr".shapeit4.phased,chr="$chr" \
--output=log/"$ref".chr"$chr".shapeit4.phasing.out \
--error=log/"$ref".chr"$chr".shapeit4.phasing.err \
/mnt/ls15/scratch/users/dingrong/swim/phasing_shapeit4_ref.sbatch

done
done


for ref in ref1.LR550 ref2.DLY550 ref3.noDLY550
do
for chr in $(seq 1 18) 
do
#eagle_target
sbatch --cpus-per-task=24 \
--mem=64G \
--time=1:59:00 \
--export=env=/mnt/ls15/scratch/users/dingrong/swim/resource.env,filepath=$filepath,input=target.chip.chr"$chr".vcf.gz,refinput="$ref".chr"$chr".eagle.phased.vcf.gz,output=target.chip.chr"$chr"."$ref".eagle.phased \
--output=log/target.chr"$chr"."$ref".eagle.phasing.out \
--error=log/target.chr"$chr"."$ref".eagle.phasing.err \
/mnt/ls15/scratch/users/dingrong/swim/phasing_eagle_target.sbatch


#shapeit4_target
sbatch -C intel16 \
--cpus-per-task=24 \
--export=env=/mnt/ls15/scratch/users/dingrong/swim/resource.env,filepath=$filepath,input=target.chip.chr"$chr".vcf.gz,refinput="$ref".chr"$chr".shapeit4.phased.vcf.gz,output=target.chr"$chr"."$ref".shapeit4.phased,chr="$chr" \
--time=1:59:00 \
--mem=32G \
--output=log/target.chr"$chr"."$ref".shapeit4.phasing.out \
--error=log/target.chr"$chr"."$ref".shapeit4.phasing.err \
/mnt/ls15/scratch/users/dingrong/swim/phasing_shapeit4_targetf.sbatch

done
done

filepath=/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550

for chr in $(seq 1 18) 
do
#beagle_target
sbatch \
--cpus-per-task=6 \
--time=1:59:00 \
--export=env=/mnt/ls15/scratch/users/dingrong/swim/resource.env,filepath=$filepath,input=target.chip.chr"$chr".vcf.gz,output=target.chip.chr"$chr".beagle.phased,tmp=$filepath/tmp,mem=32G \
--output=log/target.chr"$chr".beagle.phasing.out \
--error=log/target.chr"$chr".beagle.phasing.err \
/mnt/ls15/scratch/users/dingrong/swim/phasing_beagle.sbatch

done


###imputation

mkdir imputation_result

filepath=/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550

for ref in ref1.LR550 ref2.DLY550 ref3.noDLY550
do
for chr in $(seq 1 18) 
do

###impute5
sbatch \
--time=1:50:00 \
--export=env=/mnt/ls15/scratch/users/dingrong/swim/resource.env,filepath=$filepath,ref_panel="$ref".chr"$chr".shapeit4.phased,input=target.chr"$chr"."$ref".shapeit4.phased,output=$filepath/imputation_result/target.chr"$chr"."$ref".impute5.imputation,chr="$chr" \
--output=log/target.chr"$chr"."$ref".impute5.imputation.out \
--error=log/target.chr"$chr"."$ref".impute5.imputation.err \
/mnt/ls15/scratch/users/dingrong/swim/imputation_impute5.sbatch


###beagle5.2
sbatch \
--time=1:50:00 \
--export=env=/mnt/ls15/scratch/users/dingrong/swim/resource.env,filepath=$filepath,ref_panel="$ref".chr"$chr".beagle.phased,input=target.chip.chr"$chr".beagle.phased,output=$filepath/imputation_result/target.chr"$chr"."$ref".beagle.imputation,tmp=$filepath/tmp,mem=64G \
--output=log/target.chr"$chr"."$ref".beagle.imputation.out \
--error=log/target.chr"$chr"."$ref".beagle.imputation.err \
/mnt/ls15/scratch/users/dingrong/swim/imputation_beagle.sbatch


###minimac4
sbatch \
--time=1:50:00 \
--export=env=/mnt/ls15/scratch/users/dingrong/swim/resource.env,filepath=$filepath,ref_panel="$ref".chr"$chr".eagle.phased,input=target.chip.chr"$chr"."$ref".eagle.phased,output=$filepath/imputation_result/target.chr"$chr"."$ref".minimac4.imputation \
--output=log/target.chr"$chr"."$ref".minimac4.imputation.out \
--error=log/target.chr"$chr"."$ref".minimac4.imputation.err \
/mnt/ls15/scratch/users/dingrong/swim/imputation_minimac4.sbatch

done
done

############################################################################################
##################three ref population maf>0################################################
############################################################################################
mkdir imputation_result/012

# SORT SAMPLE NAMES
sort /mnt/ls15/scratch/users/dingrong/swim/imputationEligible/target.LR100.ind | awk '{print $0,$0}'> imputation_result/012/sampleNamesPred.txt

### SORT SAMPLE NAMES and filter imputed SNP ID
for chr in $(seq 1 18) 
do
echo $chr
#filter imputed SNP ID
bcftools query -i 'IMP=1' -f '%ID\n' imputation_result/target.chr"$chr".ref1.LR550.impute5.imputation.vcf.gz > imputation_result/012/pos_impute5.chr"$chr".txt
bcftools query -i 'IMPUTED=1' -f '%ID\n' imputation_result/target.chr"$chr".ref1.LR550.minimac4.imputation.dose.vcf.gz > imputation_result/012/pos_minimac4.chr"$chr".txt
awk -F ':'  '{print $1"_"$2"_"$3}' imputation_result/012/pos_minimac4.chr"$chr".txt > imputation_result/012/pos_minimac4.chr"$chr"_2.txt
bcftools query -i 'IMP=1' -f '%ID\n' imputation_result/target.chr"$chr".ref1.LR550.beagle.imputation.vcf.gz > imputation_result/012/pos_beagle.chr"$chr".txt
done


### get target SNP maf
mkdir miss_summary

for chr in $(seq 1 18)
do
  /mnt/research/qgg/software/plink-v1.90b6.18/plink \
  --vcf target.LR100.chr"$chr".recode.vcf \
  --freq \
  --out miss_summary/target.LR100.chr"$chr" \
  > log/target.LR100.chr"$chr".log 2>&1 &  
done

###Extract the position where maf > 0 in the target panel
for chr in $(seq 1 18)
do
echo $chr
awk '{if($5>0 && NR>1) print $2}' miss_summary/target.LR100.chr"$chr".frq >miss_summary/target_pos.LR100.chr"$chr".mafM0.txt
done

### get ref SNP maf

for chr in $(seq 1 18)
do 
  /mnt/research/qgg/software/plink-v1.90b6.18/plink \
  --vcf ref1.LR550.chr"$chr".recode.vcf \
  --freq \
  --out miss_summary/ref1.LR550.chr"$chr" \
  > log/ref1.LR550.chr"$chr".log 2>&1 &
  
  /mnt/research/qgg/software/plink-v1.90b6.18/plink \
  --vcf ref2.DLY550.chr"$chr".recode.vcf \
  --freq \
  --out miss_summary/ref2.DLY550.chr"$chr" \
  > log/ref2.DLY550.chr"$chr".log 2>&1 &
  
  /mnt/research/qgg/software/plink-v1.90b6.18/plink \
  --vcf ref3.noDLY550.chr"$chr".recode.vcf \
  --freq \
  --out miss_summary/ref3.noDLY550.chr"$chr" \
  > log/ref3.noDLY550.chr"$chr".log 2>&1 &
  
done


### Extract the positions where the maf > 0 in the 3 ref panels
for chr in $(seq 1 18)
do
echo $chr
awk '{if(NR>1 && $5>0) print $2}' miss_summary/ref1.LR550.chr"$chr".frq > miss_summary/ref1.LR550.mafM0.chr"$chr".txt
awk '{if(NR>1 && $5>0) print $2}' miss_summary/ref2.DLY550.chr"$chr".frq > miss_summary/ref2.DLY550.mafM0.chr"$chr".txt
awk '{if(NR>1 && $5>0) print $2}' miss_summary/ref3.noDLY550.chr"$chr".frq > miss_summary/ref3.noDLY550.mafM0.chr"$chr".txt

cat miss_summary/ref1.LR550.mafM0.chr"$chr".txt miss_summary/ref2.DLY550.mafM0.chr"$chr".txt miss_summary/ref3.noDLY550.mafM0.chr"$chr".txt \
 |sort |uniq -c |sort -n -r | awk '{if ($1==3) print $2}' \
 | awk -F '_'  '{print $2,$0}' |sort -k1n |awk '{print $2}' >miss_summary/final.mafM0.chr"$chr".txt
 
done

### Further extract the positions where maf is also greater than 0 in the target population
for chr in $(seq 1 18)
do
echo $chr

 cat miss_summary/final.mafM0.chr"$chr".txt  miss_summary/target_pos.LR100.chr"$chr".mafM0.txt imputation_result/012/pos_impute5.chr"$chr".txt \
 |sort |uniq -c |sort -n -r | awk '{if ($1==3) print $2}' \
 | awk -F '_'  '{print $2,$0}' |sort -k1n |awk '{print $2}' >imputation_result/012/pos_impute5.mafM0.chr"$chr".txt
 
 cat miss_summary/final.mafM0.chr"$chr".txt miss_summary/target_pos.LR100.chr"$chr".mafM0.txt imputation_result/012/pos_beagle.chr"$chr".txt \
 |sort |uniq -c |sort -n -r | awk '{if ($1==3) print $2}' \
 | awk -F '_'  '{print $2,$0}' |sort -k1n |awk '{print $2}' >imputation_result/012/pos_beagle.mafM0.chr"$chr".txt

cat miss_summary/final.mafM0.chr"$chr".txt  miss_summar/target_pos.LR100.chr"$chr".mafM0.txt imputation_result/012/pos_minimac4.chr"$chr"_2.txt \
 |sort |uniq -c |sort -n -r | awk '{if ($1==3) print $2}' | awk -F '_'  '{print $2,$0}' |sort -k1n |awk '{print $2}' \
 |  awk -F '_'  '{print $1":"$2":"$3}' | grep -F -f - imputation_result/012/pos_minimac4.chr"$chr".txt \
 >imputation_result/012/pos_minimac4.mafM0.chr"$chr".txt

cat miss_summary/final.mafM0.chr"$chr".txt  miss_summary/target_pos.LR100.chr"$chr".mafM0.txt imputation_result/012/pos_minimac4.chr"$chr"_2.txt \
 |sort |uniq -c |sort -n -r | awk '{if ($1==3) print $2}' | awk -F '_'  '{print $2,$0}' |sort -k1n |awk '{print $2}' \
 >imputation_result/012/pos_minimac4.mafM0.chr"$chr"_2.txt
done

### The maf corresponding to the filtered position
for chr in $(seq 1 18) 
do
echo $chr
###filter maf

grep -F -f imputation_result/012/pos_beagle.mafM0.chr"$chr".txt miss_summary/target.LR100.chr"$chr".frq \
|awk '{print $2,$5}' |awk -F '_'  '{print $2,$0}' |sort -k1n |awk '{print $3}'> imputation_result/012/maf_beagle.mafM0.chr"$chr".target.txt

grep -F -f imputation_result/012/pos_impute5.mafM0.chr"$chr".txt miss_summary/target.LR100.chr"$chr".frq \
|awk '{print $2,$5}' |awk -F '_'  '{print $2,$0}' |sort -k1n |awk '{print $3}'> imputation_result/012/maf_impute5.mafM0.chr"$chr".target.txt

grep -F -f imputation_result/012/pos_minimac4.mafM0.chr"$chr"_2.txt miss_summary/target.LR100.chr"$chr".frq \
|awk '{print $2,$5}' |awk -F '_'  '{print $2,$0}' |sort -k1n |awk '{print $3}'> imputation_result/012/maf_minimac4.mafM0.chr"$chr".target.txt

done

### CONVERT 012 for ref
for ref in ref1.LR550 ref2.DLY550 ref3.noDLY550
do
for chr in $(seq 1 18) 
do
sbatch \
--time=1:50:00 \
--export=imputeSof="beagle",ref=$ref,chr="$chr",THREADS=2 \
--output=log/convert012_impute.mafM0.chr"$chr"."$ref".beagle.out \
--error=log/convert012_impute.mafM0.chr"$chr"."$ref".beagle.err \
convert012_mafM0.sbatch

sbatch \
--time=1:50:00 \
--export=imputeSof="impute5",ref=$ref,chr="$chr",THREADS=2 \
--output=log/convert012_impute.mafM0.chr"$chr"."$ref".impute5.out \
--error=log/convert012_impute.mafM0.chr"$chr"."$ref".impute5.err \
convert012_mafM0.sbatch

sbatch \
--time=1:50:00 \
--export=imputeSof="minimac4",ref=$ref,chr="$chr",THREADS=2 \
--output=log/convert012_impute.mafM0.chr"$chr"."$ref".minimac4.out \
--error=log/convert012_impute.mafM0.chr"$chr"."$ref".minimac4.err \
convert012_mafM0.sbatch

done
done

### CONVERT 012 for target
for chr in $(seq 1 18) 
do
sbatch \
--time=0:30:00 \
--export=chr="$chr",THREADS=2 \
--output=log/convert012.target."$chr".out \
--error=log/convert012.target."$chr".err \
convert012.mafM0.sbatch
done

### get acc
for ref in ref1.LR550 ref2.DLY550 ref3.noDLY550
do
for chr in $(seq 1 18) 
do
sbatch \
--time=0:50:00 \
--mem=64G \
--export=imputeSof="beagle",ref="$ref",chr="$chr",numRowObs=100 \
--output=log/acc_imputation_target.chr"$chr"."$ref".beagle.out \
--error=log/acc_imputation_target.chr"$chr"."$ref".beagle.err \
acc_imputation.mafM0_target.sbatch

sbatch \
--time=0:50:00 \
--mem=64G \
--export=imputeSof="impute5",ref="$ref",chr="$chr",numRowObs=100 \
--output=log/acc_imputation_target.chr"$chr"."$ref".impute5.out \
--error=log/acc_imputation_target.chr"$chr"."$ref".impute5.err \
acc_imputation.mafM0_target.sbatch

sbatch \
--time=0:50:00 \
--mem=64G \
--export=imputeSof="minimac4",ref="$ref",chr="$chr",numRowObs=100 \
--output=log/acc_imputation_target.chr"$chr"."$ref".minimac4.out \
--error=log/acc_imputation_target.chr"$chr"."$ref".minimac4.err \
acc_imputation.mafM0_target.sbatch

done
done

###get SSC result

for ref in ref1.LR550 ref2.DLY550 ref3.noDLY550
do
sbatch \
--export=imputeSof="beagle",ref="$ref" \
--output=log/final_acc_imputation.mafM0.SSC."$ref".beagle.out \
--error=log/final_acc_imputation.mafM0.SSC."$ref".beagle.err \
final_acc_imputation.mafM0.SSC.sbatch

sbatch \
--export=imputeSof="impute5",ref="$ref" \
--output=log/final_acc_imputation.mafM0.SSC."$ref".impute5.out \
--error=log/final_acc_imputation.mafM0.SSC."$ref".impute5.err \
final_acc_imputation.mafM0.SSC.sbatch

sbatch \
--export=imputeSof="minimac4",ref="$ref" \
--output=log/final_acc_imputation.mafM0.SSC."$ref".minimac4.out \
--error=log/final_acc_imputation.mafM0.SSC."$ref".minimac4.err \
final_acc_imputation.mafM0.SSC.sbatch

done


######################
###get maf result
cat  imputation_result/summary/accuracy_target_maf.impute5.mafM0.chr*.ref1.LR550.txt |grep -v ^maf >imputation_result/summary/final_accuracy_maf.impute5.mafM0.ref1.LR550.txt
cat  imputation_result/summary/accuracy_target_maf.minimac4.mafM0.chr*.ref1.LR550.txt |grep -v ^maf >imputation_result/summary/final_accuracy_maf.minimac4.mafM0.ref1.LR550.txt
cat  imputation_result/summary/accuracy_target_maf.beagle.mafM0.chr*.ref1.LR550.txt |grep -v ^maf >imputation_result/summary/final_accuracy_maf.beagle.mafM0.ref1.LR550.txt

cat  imputation_result/summary/accuracy_target_maf.impute5.mafM0.chr*.ref2.DLY550.txt |grep -v ^maf >imputation_result/summary/final_accuracy_maf.impute5.mafM0.ref2.DLY550.txt
cat  imputation_result/summary/accuracy_target_maf.minimac4.mafM0.chr*.ref2.DLY550.txt |grep -v ^maf >imputation_result/summary/final_accuracy_maf.minimac4.mafM0.ref2.DLY550.txt
cat  imputation_result/summary/accuracy_target_maf.beagle.mafM0.chr*.ref2.DLY550.txt |grep -v ^maf >imputation_result/summary/final_accuracy_maf.beagle.mafM0.ref2.DLY550.txt

cat  imputation_result/summary/accuracy_target_maf.impute5.mafM0.chr*.ref3.noDLY550.txt |grep -v ^maf >imputation_result/summary/final_accuracy_maf.impute5.mafM0.ref3.noDLY550.txt
cat  imputation_result/summary/accuracy_target_maf.minimac4.mafM0.chr*.ref3.noDLY550.txt |grep -v ^maf >imputation_result/summary/final_accuracy_maf.minimac4.mafM0.ref3.noDLY550.txt
cat  imputation_result/summary/accuracy_target_maf.beagle.mafM0.chr*.ref3.noDLY550.txt |grep -v ^maf >imputation_result/summary/final_accuracy_maf.beagle.mafM0.ref3.noDLY550.txt


###get mafbin result
###ref1
sbatch \
--export=inputname="/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/summary/final_accuracy_maf.impute5.mafM0.ref1.LR550",outputname="/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/summary/final_accuracy_maf_bin.impute5.mafM0.ref1.LR550" \
--output=log/acc_maf_bin.swim.out \
--error=log/acc_maf_bin.swim.err \
/mnt/ls15/scratch/users/dingrong/swim/pharp_acc/acc_maf_bin.sbatch

sbatch \
--export=inputname="/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/summary/final_accuracy_maf.minimac4.mafM0.ref1.LR550",outputname="/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/summary/final_accuracy_maf_bin.minimac4.mafM0.ref1.LR550" \
--output=log/acc_maf_bin.swim.out \
--error=log/acc_maf_bin.swim.err \
/mnt/ls15/scratch/users/dingrong/swim/pharp_acc/acc_maf_bin.sbatch

sbatch \
--export=inputname="/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/summary/final_accuracy_maf.beagle.mafM0.ref1.LR550",outputname="/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/summary/final_accuracy_maf_bin.beagle.mafM0.ref1.LR550" \
--output=log/acc_maf_bin.swim.out \
--error=log/acc_maf_bin.swim.err \
/mnt/ls15/scratch/users/dingrong/swim/pharp_acc/acc_maf_bin.sbatch

###ref2
sbatch \
--export=inputname="/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/summary/final_accuracy_maf.impute5.mafM0.ref2.DLY550",outputname="/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/summary/final_accuracy_maf_bin.impute5.mafM0.ref2.DLY550" \
--output=log/acc_maf_bin.swim.out \
--error=log/acc_maf_bin.swim.err \
/mnt/ls15/scratch/users/dingrong/swim/pharp_acc/acc_maf_bin.sbatch

sbatch \
--export=inputname="/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/summary/final_accuracy_maf.minimac4.mafM0.ref2.DLY550",outputname="/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/summary/final_accuracy_maf_bin.minimac4.mafM0.ref2.DLY550" \
--output=log/acc_maf_bin.swim.out \
--error=log/acc_maf_bin.swim.err \
/mnt/ls15/scratch/users/dingrong/swim/pharp_acc/acc_maf_bin.sbatch

sbatch \
--export=inputname="/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/summary/final_accuracy_maf.beagle.mafM0.ref2.DLY550",outputname="/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/summary/final_accuracy_maf_bin.beagle.mafM0.ref2.DLY550" \
--output=log/acc_maf_bin.swim.out \
--error=log/acc_maf_bin.swim.err \
/mnt/ls15/scratch/users/dingrong/swim/pharp_acc/acc_maf_bin.sbatch

###ref3
sbatch \
--export=inputname="/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/summary/final_accuracy_maf.impute5.mafM0.ref3.noDLY550",outputname="/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/summary/final_accuracy_maf_bin.impute5.mafM0.ref3.noDLY550" \
--output=log/acc_maf_bin.swim.out \
--error=log/acc_maf_bin.swim.err \
/mnt/ls15/scratch/users/dingrong/swim/pharp_acc/acc_maf_bin.sbatch

sbatch \
--export=inputname="/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/summary/final_accuracy_maf.minimac4.mafM0.ref3.noDLY550",outputname="/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/summary/final_accuracy_maf_bin.minimac4.mafM0.ref3.noDLY550" \
--output=log/acc_maf_bin.swim.out \
--error=log/acc_maf_bin.swim.err \
/mnt/ls15/scratch/users/dingrong/swim/pharp_acc/acc_maf_bin.sbatch

sbatch \
--export=inputname="/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/summary/final_accuracy_maf.beagle.mafM0.ref3.noDLY550",outputname="/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/summary/final_accuracy_maf_bin.beagle.mafM0.ref3.noDLY550" \
--output=log/acc_maf_bin.swim.out \
--error=log/acc_maf_bin.swim.err \
/mnt/ls15/scratch/users/dingrong/swim/pharp_acc/acc_maf_bin.sbatch


