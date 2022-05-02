### get remove100LR panel

for chr in $(seq 1 18) 
do
sbatch \
--export=chr=$chr \
--output=/mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/log/remove100LR.chr"$chr".out \
--error=/mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/log/remove100LR.chr"$chr".err \
/mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/remove100LR.sbatch
done

for chr in $(seq 1 18) 
do
echo $chr
bcftools index -f -t --threads 24 /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/chr"$chr".QC.final.remove100LR.recode.vcf.gz 
done

###all sample phase
filepath=/mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf

for chr in $(seq 1 18) 
do
###shapeit4_ref
sbatch -C intel16 \
--cpus-per-task=12 \
--mem=128G \
--export=env=/mnt/ls15/scratch/users/dingrong/swim/resource.env,filepath=$filepath,input=chr"$chr".QC.final.vcf.gz,output=chr"$chr".QC.final.shapeit4.phased,chr="$chr" \
--output=log/chr"$chr".QC.final.shapeit4.phasing.out \
--error=log/chr"$chr".QC.final.shapeit4.phasing.err \
/mnt/ls15/scratch/users/dingrong/swim/phasing_shapeit4_ref.sbatch
done

for chr in $(seq 1 18) 
do
#shapeit4_target
sbatch -C intel16 \
--cpus-per-task=24 \
--export=env=/mnt/ls15/scratch/users/dingrong/swim/resource.env,filepath=$filepath,input=/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/target.chip.chr"$chr".vcf.gz,refinput=chr"$chr".QC.final.shapeit4.phased.vcf.gz,output=target.chr"$chr".shapeit4.phased,chr="$chr" \
--time=1:59:00 \
--mem=64G \
--output=log/target.chr"$chr".shapeit4.phasing.out \
--error=log/target.chr"$chr".shapeit4.phasing.err \
/mnt/ls15/scratch/users/dingrong/swim/phasing_shapeit4_targetf.sbatch
done

### remove100LR phase

for chr in $(seq 1 18) 
do
###shapeit4_ref
sbatch -C intel16 \
--cpus-per-task=12 \
--time=12:59:58 \
--mem=128G \
--export=env=/mnt/ls15/scratch/users/dingrong/swim/resource.env,filepath=$filepath,input=chr"$chr".QC.final.remove100LR.recode.vcf.gz,output=chr"$chr".QC.final.remove100LR.shapeit4.phased,chr="$chr" \
--output=$filepath/log/chr"$chr".QC.final.remove100LR.shapeit4.phasing.out \
--error=$filepath/log/chr"$chr".QC.final.remove100LR.shapeit4.phasing.err \
/mnt/ls15/scratch/users/dingrong/swim/phasing_shapeit4_ref.sbatch
done

for chr in $(seq 1 18) 
do
#shapeit4_target
sbatch -C intel16 \
--cpus-per-task=24 \
--export=env=/mnt/ls15/scratch/users/dingrong/swim/resource.env,filepath=$filepath,input=/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/target.chip.chr"$chr".vcf.gz,refinput=chr"$chr".QC.final.remove100LR.shapeit4.phased.vcf.gz,output=target.chr"$chr".remove100LR.shapeit4.phased,chr="$chr" \
--time=1:59:00 \
--mem=64G \
--output=log/target.chr"$chr".remove100LR.shapeit4.phasing.out \
--error=log/target.chr"$chr".remove100LR.shapeit4.phasing.err \
/mnt/ls15/scratch/users/dingrong/swim/phasing_shapeit4_targetf.sbatch
done

### all sample imputation
mkdir imputation_result

for chr in $(seq 1 18) 
do
###impute5
sbatch \
--time=1:50:00 \
--export=env=/mnt/ls15/scratch/users/dingrong/swim/resource.env,filepath=$filepath,ref_panel=chr"$chr".QC.final.shapeit4.phased,input=target.chr"$chr".shapeit4.phased,output=imputation_result/target.chr"$chr".impute5.imputation,chr="$chr" \
--output=log/target.chr"$chr".impute5.imputation.out \
--error=log/target.chr"$chr".impute5.imputation.err \
/mnt/ls15/scratch/users/dingrong/swim/imputation_impute5.sbatch

done

### remove100LR imputation

for chr in $(seq 1 18) 
do

###impute5
sbatch \
--time=1:50:00 \
--export=env=/mnt/ls15/scratch/users/dingrong/swim/resource.env,filepath=$filepath,ref_panel=chr"$chr".QC.final.remove100LR.shapeit4.phased,input=target.chr"$chr".remove100LR.shapeit4.phased,output=imputation_result/target.chr"$chr".remove100LR.impute5.imputation,chr="$chr" \
--output=log/target.chr"$chr".remove100LR.impute5.imputation.out \
--error=log/target.chr"$chr".remove100LR.impute5.imputation.err \
/mnt/ls15/scratch/users/dingrong/swim/imputation_impute5.sbatch

done


########################imputation in PHARP server#####################
#######################################################################

###ready for imputation
###conver target head and SSCID
cd /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/

for chr in $(seq 1 18) 
do
echo $chr
bcftools reheader -h head2target.txt target.chip.chr"$chr".vcf -o target.chip.head.chr"$chr".vcf
bcftools annotate --rename-chrs add_chrname.txt target.chip.head.chr"$chr".vcf -o target.chip.headAchr.chr"$chr".vcf
bcftools  view target.chip.headAchr.chr"$chr".vcf -Oz -o target.chip.headAchr.chr"$chr".vcf.gz 
done

###imputation in pharp server

cd /mnt/ls15/scratch/users/dingrong/swim/pharp_acc

awk '{print $2,$1}' /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/add_chrname.txt reduce_chrname.txt

for chr in $(seq 1 18) 
do
echo $chr
bcftools annotate --rename-chrs reduce_chrname.txt \
pharp_result_raw/target.chip.headAchr.chr"$chr".vcf.gz.impute.dose.vcf.gz \
-o target.chip.headAchr.chr"$chr".impute.dose.vcf.gz
done

#####################################################################################################
################################all SNV acc for SWIM and SWIM_remove100LR############################
#####################################################################################################
# SORT SAMPLE NAMES

cd /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/imputation_result/

cp  /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/012/sampleNamesPred.txt  ./012

###filter imputed SNP ID
for chr in $(seq 1 18) 
do
echo $chr

bcftools query -i 'IMP=1' -f '%ID\n' target.chr"$chr".impute5.imputation.vcf.gz > 012/pos.swim.chr"$chr".txt

bcftools query -i 'IMP=1' -f '%ID\n' target.chr"$chr".remove100LR.impute5.imputation.vcf.gz > 012/pos.swim_remove100LR.chr"$chr".txt
done

###maf
for chr in $(seq 1 18) 
do

echo $chr
###all_2259
grep -F -f 012/pos.swim.chr"$chr".txt /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/miss_summary/target.LR100.chr"$chr".frq \
|awk '{print $2,$5}' |awk -F '_'  '{print $2,$0}' |sort -k1n |awk '{print $3}'> 012/maf.swim.chr"$chr".target.txt
###remove100LR_2159
grep -F -f 012/pos.swim_remove100LR.chr"$chr".txt /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/miss_summary/target.LR100.chr"$chr".frq \
|awk '{print $2,$5}' |awk -F '_'  '{print $2,$0}' |sort -k1n |awk '{print $3}'> 012/maf.swim_remove100LR.chr"$chr".target.txt

done

### CONVERT 012
for chr in $(seq 1 18) 
do
sbatch \
--time=1:50:00 \
--export=chr="$chr" \
--output=log/convert012.chr"$chr".out \
--error=log/convert012.chr"$chr".err \
convert012.sbatch
done

### imputation acc
for chr in $(seq 1 18) 
do
###all_2259
sbatch \
--time=0:59:00 \
--mem=64G \
--export=server="swim",chr="$chr" \
--output=log/acc_imputation.swim.chr"$chr".out \
--error=log/acc_imputation.swim.chr"$chr".err \
acc_imputation.sbatch

###remove100LR_2159
sbatch \
--time=0:59:00 \
--mem=64G \
--export=server="swim_remove100LR",chr="$chr" \
--output=log/acc_imputation.swim_remove100LR.chr"$chr".out \
--error=log/acc_imputation.swim_remove100LR.chr"$chr".err \
acc_imputation.sbatch
done


###get final acc
sbatch \
--export=server="swim_remove100LR" \
--output=log/final_acc_imputation.SSC.swim_remove100LR.out \
--error=log/final_acc_imputation.SSC.swim_remove100LR.err \
final_acc_imputation.mafM0.SSC.sbatch

sbatch \
--export=server="swim" \
--output=log/final_acc_imputation.SSC.swim.out \
--error=log/final_acc_imputation.SSC.swim.err \
final_acc_imputation.mafM0.SSC.sbatch


######################
###get maf result
cat  summary/accuracy_maf.swim.chr*.txt |grep -v ^maf >summary/final_accuracy_maf.swim.txt
cat  summary/accuracy_maf.swim_remove100LR.chr*.txt |grep -v ^maf >summary/final_accuracy_maf.swim_remove100LR.txt


###get mafbin result
###all_2259
sbatch \
--export=inputname="/mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/imputation_result/summary/final_accuracy_maf.swim",outputname="/mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/imputation_result/summary/final_accuracy_maf_bin.swim.all_SNV" \
--output=log/acc_maf_bin.swim.out \
--error=log/acc_maf_bin.swim.err \
/mnt/ls15/scratch/users/dingrong/swim/pharp_acc/acc_maf_bin.sbatch

###remove100LR_2159
sbatch \
--export=inputname="/mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/imputation_result/summary/final_accuracy_maf.swim_remove100LR",outputname="/mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/imputation_result/summary/final_accuracy_maf_bin.swim_remove100LR.all_SNV" \
--output=log/acc_maf_bin.swim_remove100LR.out \
--error=log/acc_maf_bin.swim_remove100LR.err \
/mnt/ls15/scratch/users/dingrong/swim/pharp_acc/acc_maf_bin.sbatch

############################################################################################################
################################overlap SNP acc for SWIM, SWIM_remove100LR,PHARP############################
############################################################################################################
# SORT SAMPLE NAMES
cd /mnt/ls15/scratch/users/dingrong/swim/pharp_acc/

cp  /mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/imputation_result/012/sampleNamesPred.txt  ./012

###filter imputed SNP ID
for chr in $(seq 1 18) 
do
echo $chr

bcftools query -i 'IMP=1' -f '%CHROM\t%POS\t%REF\t%ALT\n' /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/imputation_result/target.chr"$chr".impute5.imputation.vcf.gz \
|awk '{print "chr"$1":"$2":"$3":"$4}'  > 012/pos_raw.swim.chr"$chr".txt

bcftools query -i 'IMP=1' -f '%CHROM\t%POS\t%REF\t%ALT\n' /mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/imputation_result/target.chr"$chr".remove100LR.impute5.imputation.vcf.gz \
|awk '{print "chr"$1":"$2":"$3":"$4}'  > 012/pos_raw.swim_remove100LR.chr"$chr".txt

bcftools query -i 'IMPUTED=1' -f '%ID\n' target.chip.headAchr.chr"$chr".impute.dose.vcf.gz \
>012/pos_raw.pharp.chr"$chr".txt

###PHARP server
cat 012/pos_raw.swim.chr"$chr".txt 012/pos_raw.swim_remove100LR.chr"$chr".txt \
012/pos_raw.pharp.chr"$chr".txt |sort |uniq -c |sort -n -r \
| awk '{if ($1==3) print $2}'| awk -F ':'  '{print $2,$0}' |sort -k1n |awk '{print $2}' \
>012/pos_final.pharp.chr"$chr".txt

###all_2259
awk -F ':'  '{print $1"_"$2"_"$3}'  012/pos_final.pharp.chr"$chr".txt | sed 's/chr//g'  >012/pos_final.swim.chr"$chr".txt

###remove100LR_2159
cp 012/pos_final.swim.chr"$chr".txt 012/pos_final.swim_remove100LR.chr"$chr".txt
done

###maf
for chr in $(seq 1 18) 
do
echo $chr
###PHARP server
grep -F -f 012/pos_final.swim.chr"$chr".txt \
/mnt/ls15/scratch/users/dingrong/swim/Landrace/imputation_LR550/miss_summary/target.LR100.chr"$chr".frq \
|awk '{print $2,$5}' |awk -F '_'  '{print $2,$0}' |sort -k1n |awk '{print $3}'> 012/maf_final.swim.chr"$chr".target.txt
###all_2259
cp 012/maf_final.swim.chr"$chr".target.txt 012/maf_final.pharp.chr"$chr".target.txt
###remove100LR_2159
cp 012/maf_final.swim.chr"$chr".target.txt 012/maf_final.swim_remove100LR.chr"$chr".target.txt
done


### CONVERT 012
for chr in $(seq 1 18) 
do
sbatch \
--time=1:50:00 \
--export=chr="$chr" \
--output=log/convert012.chr"$chr".out \
--error=log/convert012.chr"$chr".err \
/mnt/ls15/scratch/users/dingrong/swim/pharp_acc/convert012.sbatch
done

### imputation acc
for chr in $(seq 1 18) 
do
###all_2259
sbatch \
--time=0:59:00 \
--mem=64G \
--export=server="swim",chr="$chr" \
--output=log/acc_imputation.swim.chr"$chr".out \
--error=log/acc_imputation.swim.chr"$chr".err \
acc_imputation.sbatch

###remove100LR_2159
sbatch \
--time=0:59:00 \
--mem=64G \
--export=server="swim_remove100LR",chr="$chr" \
--output=log/acc_imputation.swim_remove100LR.chr"$chr".out \
--error=log/acc_imputation.swim_remove100LR.chr"$chr".err \
acc_imputation.sbatch

###PHARP server
sbatch \
--time=0:50:00 \
--mem=64G \
--export=server="pharp",chr="$chr" \
--output=log/acc_imputation.pharp.chr"$chr".out \
--error=log/acc_imputation.pharp.chr"$chr".err \
acc_imputation.sbatch
done

###get final acc
sbatch \
--export=server="pharp" \
--output=log/final_acc_imputation.SSC.pharp.out \
--error=log/final_acc_imputation.SSC.pharp.err \
final_acc_imputation.mafM0.SSC.sbatch

sbatch \
--export=server="swim" \
--output=log/final_acc_imputation.SSC.swim.out \
--error=log/final_acc_imputation.SSC.swim.err \
final_acc_imputation.mafM0.SSC.sbatch

sbatch \
--export=server="swim_remove100LR" \
--output=log/final_acc_imputation.SSC.swim_remove100LR.out \
--error=log/final_acc_imputation.SSC.swim_remove100LR.err \
final_acc_imputation.mafM0.SSC.sbatch


######################
###get maf result
cat  summary/accuracy_maf.pharp.chr*.txt |grep -v ^maf >summary/final_accuracy_maf.pharp.txt
cat  summary/accuracy_maf.swim.chr*.txt |grep -v ^maf >summary/final_accuracy_maf.swim.txt
cat  summary/accuracy_maf.swim_remove100LR.chr*.txt |grep -v ^maf >summary/final_accuracy_maf.swim_remove100LR.txt


###get mafbin result
###all_2259
sbatch \
--export=inputname="summary/final_accuracy_maf.swim",outputname="summary/final_accuracy_maf_bin.swim.overlap_SNP" \
--output=log/acc_maf_bin.swim.out \
--error=log/acc_maf_bin.swim.err \
acc_maf_bin.sbatch

###remove100LR_2159
sbatch \
--export=inputname="summary/final_accuracy_maf.swim_remove100LR",outputname="summary/final_accuracy_maf_bin.swim_remove100LR.overlap_SNP" \
--output=log/acc_maf_bin.swim_remove100LR.out \
--error=log/acc_maf_bin.swim_remove100LR.err \
acc_maf_bin.sbatch

###PHARP server
sbatch \
--export=inputname="summary/final_accuracy_maf.pharp",outputname="summary/final_accuracy_maf_bin.pharp.overlap_SNP" \
--output=log/acc_maf_bin.pharp.out \
--error=log/acc_maf_bin.pharp.err \
acc_maf_bin.sbatch
