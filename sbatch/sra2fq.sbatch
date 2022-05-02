#!/bin/bash

# Job name:
#SBATCH --job-name=sra2fastq
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=2
#
# Memory per node:
#SBATCH --mem=5G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=24:00:00
# ============================================================

source "$env"

for iter in $(echo $srr | tr '+' ' ')
do
  echo $iter
  $PREFETCH -O "$dir" --max-size 200G --force ALL --transport http "$iter" > "$dir"/$iter.prefetch.log 2>&1
  if [[ `grep err: "$dir"/$iter.prefetch.log | wc -l` -gt 0 ]] || [[ `grep failed "$dir"/$iter.prefetch.log | wc -l` -gt 0 ]]
  then
    >&2 echo "$(date +"%m-%d-%Y-%T"):error: failed to download SRA $iter."
    exit 1
  else
  $FQDUMP --split-files --origfmt -O "$dir"/$iter "$dir"/$iter/$iter.sra > "$dir"/$iter.fqdump.log 2>&1
  rm "$dir"/$iter/$iter.sra
  if [[ -e $iter ]]
  then
    rm -r $iter
  fi
  fi
done

for iter in $(echo $srr | tr '+' ' ')
do
  if [[ -e "$dir"/"$iter"/"$iter"_1.fastq ]]
  then
    nf=`head -n 1 "$dir"/"$iter"/"$iter"_1.fastq | awk '{print $1}' | awk -F ":" '{print NF}'`
    
    if [[ -e "$dir"/"$iter"/"$iter"_2.fastq ]]
    then
      layout="pe"
      # when paired, sort and then join to make sure that the files actually pair
      # check for read name pairing
      unpline=`paste <(cat "$dir"/"$iter"/"$iter"_1.fastq | paste - - - - | cut -f 1 | sed 's/\/1//') <(cat "$dir"/"$iter"/"$iter"_2.fastq | paste - - - - | cut -f 1 | sed 's/\/1//') | awk '$1 != $2' | head | wc -l`
      
      if [[ $unpline -gt 0 ]]
      then
        cat "$dir"/"$iter"/"$iter"_1.fastq | paste - - - - | sort -k1,1 > "$dir"/"$iter"/"$iter"_1.sorted.fastq
        cat "$dir"/"$iter"/"$iter"_2.fastq | paste - - - - | sort -k1,1 > "$dir"/"$iter"/"$iter"_2.sorted.fastq
        join -t $'\t' "$dir"/"$iter"/"$iter"_1.sorted.fastq "$dir"/"$iter"/"$iter"_2.sorted.fastq | cut -f 1-4 | tr '\t' '\n' > "$dir"/"$iter"/"$iter"_1.fastq
        join -t $'\t' "$dir"/"$iter"/"$iter"_1.sorted.fastq "$dir"/"$iter"/"$iter"_2.sorted.fastq | cut -f 1,5-7 | tr '\t' '\n' > "$dir"/"$iter"/"$iter"_2.fastq
        rm "$dir"/"$iter"/"$iter"_1.sorted.fastq "$dir"/"$iter"/"$iter"_2.sorted.fastq
      fi
      
    else
      layout="se"
    fi
    
    if [[ $nf -eq 7 ]] # 7 fields
    then
      cat "$dir"/"$iter"/"$iter"_1.fastq | paste - - - - | cut -f 1 | sed 's/^@//' | awk -F ":" '{print $1"."$2"."$3"."$4" "$1" "$2" "$3" "$4" @"$1":"$2":"$3":"$4" '$dir'/'$sample'/'$iter'/'$iter'"}' | uniq | sort | uniq | awk '{print "'$layout'."$0}' >> "$dir"/file.info
    elif [[ $nf -eq 5 ]]
    then
      cat "$dir"/"$iter"/"$iter"_1.fastq | paste - - - - | cut -f 1 | sed 's/^@//' | awk -F ":" '{print "-.-."$1"."$2" - - "$1" "$2" @"$1":"$2" '$dir'/'$sample'/'$iter'/'$iter'"}' | uniq | sort | uniq | awk '{print "'$layout'."$0}' >> "$dir"/file.info
    else
      cat "$dir"/"$iter"/"$iter"_1.fastq | paste - - - - | head | cut -f 1 | sed 's/^@//' | awk -F ":" '{print "-.-.'$iter'.1 - - '$iter' 1 - '$dir'/'$sample'/'$iter'/'$iter'"}' | uniq | sort | uniq | awk '{print "'$layout'."$0}' >> "$dir"/file.info
      cat "$dir"/"$iter"/"$iter"_1.fastq | paste - - - - | sed 's/^@/@'$iter'-/' | awk -F "\t" '{print $1"\t"$2"\t+\t"$4}' | tr '\t' '\n' > "$dir"/"$iter"/"$iter"_1.tmp.fastq
      mv "$dir"/"$iter"/"$iter"_1.tmp.fastq "$dir"/"$iter"/"$iter"_1.fastq
      
      if [[ "$layout" = "pe" ]]
      then
        cat "$dir"/"$iter"/"$iter"_2.fastq | paste - - - - | sed 's/^@/@'$iter'-/' | awk -F "\t" '{print $1"\t"$2"\t+\t"$4}' | tr '\t' '\n' > "$dir"/"$iter"/"$iter"_2.tmp.fastq
        mv "$dir"/"$iter"/"$iter"_2.tmp.fastq "$dir"/"$iter"/"$iter"_2.fastq
      fi
      
    fi
  fi
done

# split fastq
# ============================================================

sed 's/ /\t/g' "$dir"/file.info | sort -k1,1 | $BEDTOOLS groupby -g 1 -c 2,3,4,5,6,7 -o distinct,distinct,distinct,distinct,distinct,distinct | sed 's/\t/ /g' > "$dir"/file.info.tmp
mv "$dir"/file.info.tmp "$dir"/file.info

while read line           
do           
  # output fastq file
	pre=$(echo $line | cut -d " " -f 1)
	flc=$(echo $line | cut -d " " -f 4)
	lane=$(echo $line | cut -d " " -f 5)
  reg=$(echo $line | cut -d " " -f 6)
  file=$(echo $line | cut -d " " -f 7)
  
  if [[ "$pre" =~ ^"pe.".* ]]
  then
    if [[ "$reg" = "-" ]]
    then
      echo $file, | sed 's/,/_1\.fastq /g' | xargs cat | gzip > "$dir"/"$pre"_1.fastq.gz
      echo $file, | sed 's/,/_2\.fastq /g' | xargs cat | gzip > "$dir"/"$pre"_2.fastq.gz
    else   
      echo $file, | sed 's/,/_1\.fastq /g' | xargs cat | grep -A 3 ^"$reg" --no-group-separator | gzip > "$dir"/"$pre"_1.fastq.gz
      echo $file, | sed 's/,/_2\.fastq /g' | xargs cat | grep -A 3 ^"$reg" --no-group-separator | gzip > "$dir"/"$pre"_2.fastq.gz
    fi
  else
    if [[ "$reg" = "-" ]]
    then
      echo $file, | sed 's/,/_1\.fastq /g' | xargs cat | gzip > "$dir"/"$pre"_1.fastq.gz
    else
      echo $file, | sed 's/,/_1\.fastq /g' | xargs cat | grep -A 3 ^"$reg" --no-group-separator | gzip > "$dir"/"$pre"_1.fastq.gz
    fi
  fi
done < "$dir"/file.info

chmod 440 "$dir"/*.fastq.gz

for iter in $(echo $srr | tr '+' ' ')
do
  rm -r "$dir"/$iter
done

echo $(date +"%m-%d-%Y-%T"):info: done.main.process
