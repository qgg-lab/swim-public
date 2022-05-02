#! /bin/bash
# ============================================================

# get command line options
# ============================================================

args=`getopt -o "r:,p:,m:,t:,o:" -l "resource:,ped:,map:,tmp:,out:" -- "$@"`
>&2 echo "command arguments given: $args"
eval set -- "$args"

# parse arguments
# ============================================================

while true;
do
  case $1 in

    -r|--resource)
      resource=$2
      shift 2;;

    -p|--ped)
			ped=$2
			shift 2;; 

  	-m|--map)
			map=$2
			shift 2;;

    -t|--tmp)
			tmp=$2
			shift 2;;

    -o|--out)
		  out=$2
			shift 2;;
    
    --)
      shift
      break;;

  esac
done

# start time
# ============================================================

startTime=$(date)
echo -n "" > $out/runStatus.out

# prepare directory
# ============================================================

if [[ -e $tmp ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):error: $tmp exists."
  exit 1
else
  mkdir $tmp
fi

if [[ -e $out ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):error: $out exists."
  exit 1
else
  mkdir $out
fi

# check for file availability
# ============================================================

if [[ -e $resource ]]
then
  echo "$(date +"%m-%d-%Y-%T"):info: found resource dir."
else
  echo "$(date +"%m-%d-%Y-%T"):error: cannot find resource dir." >> $out/runStatus.out
  exit 1
fi

source $resource/swim.resource.env

# exit when command fails
# ============================================================

set -o pipefail

# =============================
# = 1. convert ped/map to vcf =
# =============================

# 1a) get plink binary format
$PLINK --silent --allow-no-sex --allow-extra-chr --ped $ped --map $map --make-bed --out $tmp/plink > $tmp/plink.log 2>&1

if [ $? -gt 0 ]
then
  echo "plink error 1:" `grep Error tmp/plink.log | sed 's/Error: //'` >> $out/runStatus.out
  exit 1
fi

nind1=`wc -l $tmp/plink.fam | awk '{print $1}'`
nsnp1=`wc -l $tmp/plink.bim | awk '{print $1}'`

echo "$(date +"%m-%d-%Y-%T"):progress: $nind1 individuals and $nsnp1 variants are imported."

# 1b) eliminate snps and individuals with high missing rate
$PLINK --silent --allow-no-sex --allow-extra-chr --bfile $tmp/plink --mind 0.2 --geno 0.2 --maf 1e-20 --missing --make-bed --out $tmp/plink.filter1 > $tmp/plink.filter1.log 2>&1

if [ $? -gt 0 ]
then
  echo "plink error 2:" `grep Error tmp/plink.filer1.log | sed 's/Error: //'` >> $out/runStatus.out
  exit 1
fi

nind2=`wc -l $tmp/plink.filter1.fam | awk '{print $1}'`
nsnp2=`wc -l $tmp/plink.filter1.bim | awk '{print $1}'`

genoSNP=`grep "variants removed due to missing genotype data" $tmp/plink.filter1.log | awk '{print $1}'`
mafSNP=`grep "variants removed due to minor allele threshold" $tmp/plink.filter1.log | awk '{print $1}'`

echo "$(date +"%m-%d-%Y-%T"):progress: $nind2 individuals remain after removing missing rate > 0.2."
echo "$(date +"%m-%d-%Y-%T"):progress: $genoSNP variants removed due to missing rate > 0.2, $mafSNP variants removed due to MAF == 0."

# 1c) now remove duplicates (retain highest call rate one) or variants that do not conform with the alleles

awk '{print $2"\t"$5"\t"$6}' $tmp/plink.filter1.bim | sort -k1,1 | join -t $'\t' - <(tail -n+2 $tmp/plink.filter1.lmiss | awk '{print $2"\t"$5}' | sort -k1,1) | join -t $'\t' - $CHIP | awk '{print $5"_"$6"_"$7"_"$8"\t"$0}' | awk '$6 >=1 && $6 <= 18' | perl $CHECKALLELES | cut -f 2- | sort -k1,1 -k4,4g | $BEDTOOLS groupby -g 1 -c 2,3,4,5,6,7,8,9 -o first,first,first,first,first,first,first,first 2> $tmp/alleleCheck.bedtools.log | sort -k1,1 > $tmp/plink.filter1.alleles.retain

# bedtools does not exit with status > 0

if [ `grep ERROR $tmp/alleleCheck.bedtools.log | wc -l | awk '{print $1}'` -gt 0 ]
then
  echo "check allele error: allele check bedtools error." >> $out/runStatus.out
  exit 1
fi

# 1d) flip alleles

cut -f 1 $tmp/plink.filter1.alleles.retain > $tmp/plink.filter1.retain.snp
awk '$9 == "flip" {print $1}' $tmp/plink.filter1.alleles.retain > $tmp/plink.filter1.flip.snp

$PLINK --silent --allow-no-sex --allow-extra-chr --bfile $tmp/plink.filter1 --extract $tmp/plink.filter1.retain.snp --flip $tmp/plink.filter1.flip.snp --make-bed --out $tmp/plink.filter2

if [ $? -gt 0 ]
then
  echo "plink error 3:" `grep Error tmp/plink.filter2.log | sed 's/Error: //'` >> $out/runStatus.out
  exit 1
fi

# fix bim for reference position
cat -n $tmp/plink.filter2.bim | awk '{print $3"\t"$0}' | sort -k1,1 | join -t $'\t' - $tmp/plink.filter1.alleles.retain | sort -k2,2n | awk '{print $12"\t"$1"\t0\t"$13"\t"$7"\t"$8}' > $tmp/tmp.bim
mv $tmp/tmp.bim $tmp/plink.filter2.bim 
# fix split chr problem
$PLINK --silent --allow-no-sex --allow-extra-chr --bfile $tmp/plink.filter2 --make-bed --out $tmp/plink.filter3

if [ $? -gt 0 ]
then
  echo "plink error 4:" `grep Error tmp/plink.filter3.log | sed 's/Error: //'` >> $out/runStatus.out
  exit 1
fi

retainSNP=`wc -l $tmp/plink.filter1.retain.snp | awk '{print $1}'`
flipSNP=`wc -l $tmp/plink.filter1.flip.snp | awk '{print $1}'`

echo "$(date +"%m-%d-%Y-%T"):progress: $retainSNP autosomal variants remain after filtering duplicates, ambiguous strand, non-matching allele variants."
echo "$(date +"%m-%d-%Y-%T"):progress: of these, $flipSNP variants were flipped to the forward strand."

# 1e) convert to VCF and set reference allele, split vcf by chromosome

cut -f 1,7 $tmp/plink.filter1.alleles.retain > $tmp/plink.filter3.ref.allele

snpNumber=()

for chr in `seq 1 1 18`
do
	snpNumber[$chr]=0
done

cut -f 1 $tmp/plink.filter3.bim | sort -k1,1g | uniq -c | awk '{print $2"\t"$1}' > $tmp/plink.filter3.chr.snp.number


while read line
do
	chr=`echo $line | awk '{print $1}'`
	nSNPs=`echo $line | awk '{print $2}'`
	snpNumber[$chr]=$nSNPs
	if [ $nSNPs -gt 20 ]
	then
    echo "$(date +"%m-%d-%Y-%T"):progress: making VCF for chromosome $chr with $nSNPs variants."
    $PLINK --silent --allow-no-sex --allow-extra-chr --bfile $tmp/plink.filter3 --chr $chr --recode vcf --a2-allele $tmp/plink.filter3.ref.allele 2 1 --out $tmp/plink.final.chr"$chr"
    
    if [ $? -gt 0 ]
    then
      echo "plink error 5 for chr $chr:" `grep Error tmp/plink.final.chr"$chr" | sed 's/Error: //'` >> $out/runStatus.out
      exit 1
    fi
    
    bgzip $tmp/plink.final.chr"$chr".vcf
    $BCFTOOLS index $tmp/plink.final.chr"$chr".vcf.gz
	fi
done < $tmp/plink.filter3.chr.snp.number

echo "$(date +"%m-%d-%Y-%T"):progress: finished preparing VCFs."

# 2. phase and impute
# ============================================================

chrSuccess=()

if [ $nind2 -gt 10000 ]
then
  
	for chr in `seq 1 1 18`
	do
	  
		if [ ${snpNumber[$chr]} -gt 20 ]
		then
      echo "$(date +"%m-%d-%Y-%T"):progress: phasing chromosome $chr."
      $SHAPEIT --input $tmp/plink.final.chr"$chr".vcf.gz --thread 10 --map $SHAPEITMAP/genMap_1cMperMb_"$chr".txt --region $chr --reference $REFHAP/chr"$chr".QC.final.shapeit4.phased.vcf.gz --pbwt-mdr 0.05 --pbwt-mac 2 --effective-size 100 --output $tmp/chr"$chr".shapeit4.phased.vcf.gz > $tmp/plink.final.chr"$chr".vcf.shapeit4.phase.log 2>&1
  
      if [ $? -gt 0 ]
      then
        echo "shapeit error for chr $chr." >> $out/runStatus.out
        exit 1
      fi
  
	    $BCFTOOLS index $tmp/chr"$chr".shapeit4.phased.vcf.gz
	    echo "$(date +"%m-%d-%Y-%T"):progress: finished phasing chromosome $chr."
		
		  # loop through regions for each chromosome
      while read line
      do
        chrInt=`echo $line | awk '{print $1":"$2"-"$3}'`
        chrIntName=`echo $line | awk '{print $4}'`
				chrIntStart=`echo $line | awk '{print $2}'`
				chrIntEnd=`echo $line | awk '{print $3}'`
   	    
				chrIntnSNP=`awk '$1 == '$chr' && $4 >= '$chrIntStart' && $4 <= '$chrIntEnd'' $tmp/plink.filter3.bim | wc -l | awk '{print $1}'`
				
				echo "$(date +"%m-%d-%Y-%T"):progress: imputing chromosome $chr region $chrInt."
				
				if [ $chrIntnSNP -gt 10 ]
				then
					$IMPUTE --h $REFHAP/chr"$chr".QC.final.shapeit4.phased.vcf.gz --g $tmp/chr"$chr".shapeit4.phased.vcf.gz --r $chrInt --b 300 --ne 100 --o $tmp/chr"$chr".swim.imputed.$chrIntName.vcf.gz --out-gp-field --out-ap-field --threads 10 > $tmp/chr"$chr".impute5.$chrIntName.log 2>&1
	        if [ $? -gt 0 ]
	        then
	          echo "impute error for chr $chr region $chrInt." >> $out/runStatus.out
	          exit 1
	        fi
					echo $tmp/chr"$chr".swim.imputed.$chrIntName.vcf.gz >> $tmp/chr"$chr".impute.file.list
				fi
				
      done < <(awk '$1 == '$chr'{print $0"\t"NR}' $resource/all.chr.int.bed)
		
 	    echo "$(date +"%m-%d-%Y-%T"):progress: merging imputed genomes for chromosome $chr."
    
      $BCFTOOLS concat -f $tmp/chr"$chr".impute.file.list -O z -o $out/chr"$chr".swim.imputed.vcf.gz > $tmp/concat.chr"$chr".log 2>&1
    
      if [ $? -gt 0 ]
      then
        echo "bcftools concat error for chr $chr." >> $out/runStatus.out
        exit 1
      fi
      chrSuccess+=("$chr")
			
		fi
    
	done

else

  for chr in `seq 1 1 18`
  do
	  if [ ${snpNumber[$chr]} -gt 20 ]
	  then
      echo "$(date +"%m-%d-%Y-%T"):progress: phasing chromosome $chr."
      $SHAPEIT --input $tmp/plink.final.chr"$chr".vcf.gz --thread 10 --map $SHAPEITMAP/genMap_1cMperMb_"$chr".txt --region $chr --reference $REFHAP/chr"$chr".QC.final.shapeit4.phased.vcf.gz --pbwt-mdr 0.05 --pbwt-mac 2 --effective-size 100 --output $tmp/chr"$chr".shapeit4.phased.vcf.gz > $tmp/plink.final.chr"$chr".vcf.shapeit4.phase.log 2>&1
    
      if [ $? -gt 0 ]
      then
        echo "shapeit error for chr $chr." >> $out/runStatus.out
        exit 1
      fi
    
		  $BCFTOOLS index $tmp/chr"$chr".shapeit4.phased.vcf.gz
		  echo "$(date +"%m-%d-%Y-%T"):progress: finished phasing chromosome $chr."
      echo "$(date +"%m-%d-%Y-%T"):progress: imputing chromosome $chr."
		  $IMPUTE --h $REFHAP/chr"$chr".QC.final.shapeit4.phased.vcf.gz --g $tmp/chr"$chr".shapeit4.phased.vcf.gz --r $chr --b 300 --ne 100 --o $out/chr"$chr".swim.imputed.vcf.gz --out-gp-field --out-ap-field --threads 10 > $tmp/chr"$chr".impute5.log 2>&1
    
      if [ $? -gt 0 ]
      then
        echo "impute error for chr $chr." >> $out/runStatus.out
        exit 1
      fi
    
      bothSNP=`grep -A 1 "Imputation region" "$tmp"/chr"$chr".impute5.log | tail -n 1 | sed 's/.*\[//' | awk -F "/" '{print $2}'`
      refSNP=`grep -A 1 "Imputation region" "$tmp"/chr"$chr".impute5.log | tail -n 1 | sed 's/.*\[//' | awk -F "/" '{print $1}'`
      echo "$(date +"%m-%d-%Y-%T"):progress: $bothSNP variants on chromosome $chr from chip and ref panel with consistent alleles."
      echo "$(date +"%m-%d-%Y-%T"):progress: imputed $refSNP varaints on chromosome $chr."
		  echo "$(date +"%m-%d-%Y-%T"):progress: finished imputing chromosome $chr."
      chrSuccess+=("$chr")
	  fi
  done

fi

echo "success" >> $out/runStatus.out
join=$(IFS=, ; echo "${chrSuccess[*]}")
echo $join >> $out/runStatus.out
echo $nind2 >> $out/runStatus.out

endTime=$(date)

echo "# run summary generated on $endTime" >> $out/runSummary.txt
echo "start: $startTime" >> $out/runSummary.txt
echo "end: $endTime" >> $out/runSummary.txt
echo "n_individuals_imported: $nind1" >> $out/runSummary.txt
echo "n_snps_imported: $nsnp1" >> $out/runSummary.txt
echo "n_individuals_retained: $nind2" >> $out/runSummary.txt
echo "n_snps_retained: $retainSNP" >> $out/runSummary.txt

rm -r $tmp
