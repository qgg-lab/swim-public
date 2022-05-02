# =================
# = combine gvcfs =
# =================

# in /mnt/home/huangw53/scratch2/swim/gvcfs
# ============================================================

ls /mnt/ufs18/scratch/huangw53/pigWGS/ | ~/qgg/software/parallel-20200722/bin/parallel -j 48 "bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/ufs18/scratch/huangw53/pigWGS --tmp /mnt/ufs18/scratch/huangw53/tmp --out /mnt/ufs18/scratch/huangw53/swim/gvcfs/ --sample {} > /mnt/ufs18/scratch/huangw53/swim/gvcfs/log/{}.out 2>&1" &

# when it's done, run a test to see if it is successful
# ============================================================

for sample in $(ls /mnt/ufs18/scratch/huangw53/pigWGS/)
do
  if [ `grep error log/$sample.out | wc -l` -gt 0 ] || [ `grep haplotype log/$sample.out | wc -l` -eq 0 ]
  then
    rm -r $sample
    echo $sample
  fi
done > 20201010.rerun.sample

cat 20201010.rerun.sample | ~/qgg/software/parallel-20200722/bin/parallel -j 48 "bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/ufs18/scratch/huangw53/pigWGS --tmp /mnt/ufs18/scratch/huangw53/tmp --out /mnt/ufs18/scratch/huangw53/swim/gvcfs/ --sample {} > /mnt/ufs18/scratch/huangw53/swim/gvcfs/log/{}.out 2>&1" &

# re-run ambiguous sex after updating the sex check sbatch
# ============================================================

for sample in $(ls log/*.out | sed 's/log\///' | sed 's/\.out//')
do
  if [ `grep amb log/$sample.out | wc -l` -gt 0 ] || [ `grep error log/$sample.out | wc -l` -gt 0 ] || [ `grep "completed haplotype" log/$sample.out | wc -l` -eq 0 ]
  then
    echo $sample
    rm -r $sample
    rm $sample.out
  else
    echo $sample >> gvcf.success.txt
  fi
done

ls /mnt/ufs18/scratch/huangw53/pigWGS/ | grep "yorkshire\|S21_\|S22_" | sort | comm -23 - <(sort gvcf.success.txt) > 20201214.run.sample

cat 20201214.run.sample | ~/qgg/software/parallel-20200722/bin/parallel -j 16 "bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/ufs18/scratch/huangw53/pigWGS --tmp /mnt/ufs18/scratch/huangw53/tmp --out /mnt/ufs18/scratch/huangw53/swim/gvcfs/ --sample {} > /mnt/ufs18/scratch/huangw53/swim/gvcfs/log/{}.out 2>&1" &

for sample in $(cat 20201214.run.sample)
do
  if [ `grep amb log/$sample.out | wc -l` -gt 0 ] || [ `grep error log/$sample.out | wc -l` -gt 0 ] || [ `grep "completed haplotype" log/$sample.out | wc -l` -eq 0 ]
  then
    echo $sample
    rm -r $sample
    rm $sample.out
  else
    echo $sample >> gvcf.success.txt
  fi
done

# 03092021: run SRA samples, get the list for all SRA
# ============================================================

cp /mnt/ufs18/scratch/huangw53/pigWGS/sra03092021success.txt 03092021.run.sample
cat 03092021.run.sample | ~/qgg/software/parallel-20200722/bin/parallel -j 48 "bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/ufs18/scratch/huangw53/pigWGS --tmp /mnt/ufs18/scratch/huangw53/tmp --out /mnt/ufs18/scratch/huangw53/swim/gvcfs/ --sample {} > /mnt/ufs18/scratch/huangw53/swim/gvcfs/log/{}.out 2>&1" &

# 04012021: there was problem with the HPCC in the middle
# ============================================================

for sample in $(ls log/*.out | sed 's/log\///' | sed 's/\.out//')
do
  if [ `grep amb log/$sample.out | wc -l` -gt 0 ] || [ `grep error log/$sample.out | wc -l` -gt 0 ] || [ `grep "completed haplotype" log/$sample.out | wc -l` -eq 0 ]
  then
    echo $sample
    rm -r $sample
    rm log/$sample.out
  else
    echo $sample >> gvcf.04012021.success.txt
  fi
done

comm -23 <(sort 03092021.run.sample) <(sort gvcf.04012021.success.txt) | ~/qgg/software/parallel-20200722/bin/parallel -j 64 "bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/ufs18/scratch/huangw53/pigWGS --tmp /mnt/ufs18/scratch/huangw53/tmp --out /mnt/ufs18/scratch/huangw53/swim/gvcfs/ --sample {} > /mnt/ufs18/scratch/huangw53/swim/gvcfs/log/{}.out 2>&1" &

# 04042021: another problem with the HPCC in the middle
# ============================================================

for sample in $(ls log/*.out | sed 's/log\///' | sed 's/\.out//')
do
  if [ `grep amb log/$sample.out | wc -l` -gt 0 ] || [ `grep error log/$sample.out | wc -l` -gt 0 ] || [ `grep "completed haplotype" log/$sample.out | wc -l` -eq 0 ]
  then
    echo $sample
    rm -r $sample
    rm log/$sample.out
  else
    echo $sample >> gvcf.04042021.success.txt
  fi
done

comm -23 <(sort 03092021.run.sample) <(sort gvcf.04042021.success.txt) | ~/qgg/software/parallel-20200722/bin/parallel -j 64 "bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/ufs18/scratch/huangw53/pigWGS --tmp /mnt/ufs18/scratch/huangw53/tmp --out /mnt/ufs18/scratch/huangw53/swim/gvcfs/ --sample {} > /mnt/ufs18/scratch/huangw53/swim/gvcfs/log/{}.out 2>&1" &


# 04072021: too many failures, update scripts
# to increase memory for bwa and merge
# ============================================================

for sample in $(ls log/*.out | sed 's/log\///' | sed 's/\.out//')
do
  if [ `grep amb log/$sample.out | wc -l` -gt 0 ] || [ `grep error log/$sample.out | wc -l` -gt 0 ] || [ `grep "completed haplotype" log/$sample.out | wc -l` -eq 0 ]
  then
    echo $sample
    rm -r $sample
    rm log/$sample.out
  else
    echo $sample >> gvcf.04072021.success.txt
  fi
done

comm -23 <(sort 03092021.run.sample) <(sort gvcf.04072021.success.txt) | ~/qgg/software/parallel-20200722/bin/parallel -j 64 "bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/ufs18/scratch/huangw53/pigWGS --tmp /mnt/ufs18/scratch/huangw53/tmp --out /mnt/ufs18/scratch/huangw53/swim/gvcfs/ --sample {} > /mnt/ufs18/scratch/huangw53/swim/gvcfs/log/{}.out 2>&1" &

# 05062021: another batch of gvcf run
# ============================================================

for sample in $(ls log/*.out | sed 's/log\///' | sed 's/\.out//')
do
  if [ `grep amb log/$sample.out | wc -l` -gt 0 ] || [ `grep error log/$sample.out | wc -l` -gt 0 ] || [ `grep "completed haplotype" log/$sample.out | wc -l` -eq 0 ]
  then
    echo $sample
    rm -r $sample
    rm log/$sample.out
  else
    echo $sample >> gvcf.05062021.success.txt
  fi
done

comm -23 <(sort 03092021.run.sample) <(sort gvcf.05062021.success.txt) | ~/qgg/software/parallel-20200722/bin/parallel -j 64 "bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/ufs18/scratch/huangw53/pigWGS --tmp /mnt/ufs18/scratch/huangw53/tmp --out /mnt/ufs18/scratch/huangw53/swim/gvcfs/ --sample {} > /mnt/ufs18/scratch/huangw53/swim/gvcfs/log/{}.out 2>&1" &

# 05282021: try to finish all SRA samples
# ============================================================

for sample in $(ls log/*.out | sed 's/log\///' | sed 's/\.out//')
do
  if [ `grep amb log/$sample.out | wc -l` -gt 0 ] || [ `grep error log/$sample.out | wc -l` -gt 0 ] || [ `grep "completed haplotype" log/$sample.out | wc -l` -eq 0 ]
  then
    echo $sample
    rm -r $sample
    rm log/$sample.out
  else
    echo $sample >> gvcf.05282021.success.txt
  fi
done

comm -23 <(sort 03092021.run.sample) <(sort gvcf.05282021.success.txt) | comm -23 - <(cut -f 1 bad.sample | sort) | ~/qgg/software/parallel-20200722/bin/parallel -j 64 "bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/ufs18/scratch/huangw53/pigWGS --tmp /mnt/ufs18/scratch/huangw53/tmp --out /mnt/ufs18/scratch/huangw53/swim/gvcfs/ --sample {} > /mnt/ufs18/scratch/huangw53/swim/gvcfs/log/{}.out 2>&1" &

# 06062021: run gvcfs for the new SRA sequences
# ============================================================

for sample in $(ls log/*.out | sed 's/log\///' | sed 's/\.out//')
do
  if [ `grep amb log/$sample.out | wc -l` -gt 0 ] || [ `grep error log/$sample.out | wc -l` -gt 0 ] || [ `grep "completed haplotype" log/$sample.out | wc -l` -eq 0 ]
  then
    echo $sample
    rm -r $sample
    rm log/$sample.out
  else
    echo $sample >> gvcf.06062021.success.txt
  fi
done

comm -23 <(sort 03092021.run.sample) <(sort gvcf.06062021.success.txt) | comm -23 - <(cut -f 1 bad.sample | sort) | ~/qgg/software/parallel-20200722/bin/parallel -j 6 "bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/ufs18/scratch/huangw53/pigWGS --tmp /mnt/ufs18/scratch/huangw53/tmp --out /mnt/ufs18/scratch/huangw53/swim/gvcfs/ --sample {} > /mnt/ufs18/scratch/huangw53/swim/gvcfs/log/{}.out 2>&1" &

cp /mnt/ufs18/scratch/huangw53/pigWGS/sra.srs.05032021 06062021.run.sample

comm -23 <(sort 06062021.run.sample) <(sort gvcf.06062021.success.txt) | comm -23 - <(cut -f 1 bad.sample | sort) | ~/qgg/software/parallel-20200722/bin/parallel -j 74 "bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/ufs18/scratch/huangw53/pigWGS --tmp /mnt/ufs18/scratch/huangw53/tmp --out /mnt/ufs18/scratch/huangw53/swim/gvcfs/ --sample {} > /mnt/ufs18/scratch/huangw53/swim/gvcfs/log/{}.out 2>&1" &

# 06092021: run gvcfs for new seq
# ============================================================

ls /mnt/gs18/scratch/users/longnany/wen/pigWGS2/ | grep -v log | grep -v split > 06092021.run.sample
sort 06092021.run.sample | comm -23 - <(cut -f 1 bad.sample | sort) | ~/qgg/software/parallel-20200722/bin/parallel -j 96 "bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/gs18/scratch/users/longnany/wen/pigWGS2 --tmp /mnt/ufs18/scratch/huangw53/tmp --out /mnt/ufs18/scratch/huangw53/swim/gvcfs/ --sample {} > /mnt/ufs18/scratch/huangw53/swim/gvcfs/log/{}.out 2>&1" &

# delete success SRA samples
# ============================================================

for sample in $(grep "SRS\|ERS" gvcf.06062021.success.txt)
do
  echo $sample
  rm -rf /mnt/ufs18/scratch/huangw53/pigWGS/$sample
done

# get sucessful runs
# ============================================================

for sample in $(ls log/*.out | sed 's/log\///' | sed 's/\.out//')
do
  if [ `grep amb log/$sample.out | wc -l` -gt 0 ] || [ `grep error log/$sample.out | wc -l` -gt 0 ] || [ `grep "completed haplotype" log/$sample.out | wc -l` -eq 0 ]
  then
    echo $sample
    rm -r $sample
    rm log/$sample.out
  else
    echo $sample >> gvcf.06182021.success.txt
  fi
done

for sample in $(comm -12 <(sort 06092021.run.sample) <(sort gvcf.06182021.success.txt))
do
  echo $sample
  rm -rf /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
done

# continue the previous iteration on 06092021
# ============================================================

sort 06092021.run.sample | comm -23 - <(cut -f 1 bad.sample | cat - gvcf.06182021.success.txt | sort) | ~/qgg/software/parallel-20200722/bin/parallel -j 96 "bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/gs18/scratch/users/longnany/wen/pigWGS2 --tmp /mnt/ufs18/scratch/huangw53/tmp --out /mnt/ufs18/scratch/huangw53/swim/gvcfs/ --sample {} > /mnt/ufs18/scratch/huangw53/swim/gvcfs/log/{}.out 2>&1" &

# 06/30/2021: some problem with HPCC, all jobs stopped
# ============================================================

for sample in $(ls log/*.out | sed 's/log\///' | sed 's/\.out//')
do
  if [ `grep amb log/$sample.out | wc -l` -gt 0 ] || [ `grep error log/$sample.out | wc -l` -gt 0 ] || [ `grep "completed haplotype" log/$sample.out | wc -l` -eq 0 ]
  then
    echo $sample
    rm -r $sample
    rm log/$sample.out
  else
    echo $sample >> gvcf.06302021.success.txt
  fi
done

sort 06092021.run.sample | comm -23 - <(cut -f 1 bad.sample | cat - gvcf.06302021.success.txt | sort) | ~/qgg/software/parallel-20200722/bin/parallel -j 96 "bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/gs18/scratch/users/longnany/wen/pigWGS2 --tmp /mnt/ufs18/scratch/huangw53/tmp --out /mnt/ufs18/scratch/huangw53/swim/gvcfs/ --sample {} > /mnt/ufs18/scratch/huangw53/swim/gvcfs/log/{}.out 2>&1" &

# 07/02/2021: hopefully final run of re-running
# ============================================================

for sample in $(ls log/*.out | sed 's/log\///' | sed 's/\.out//')
do
  if [ `grep amb log/$sample.out | wc -l` -gt 0 ] || [ `grep error log/$sample.out | wc -l` -gt 0 ] || [ `grep "completed haplotype" log/$sample.out | wc -l` -eq 0 ]
  then
    echo $sample
    rm -r $sample
    rm log/$sample.out
  else
    echo $sample >> gvcf.07022021.success.txt
  fi
done

sort pigWGS.sample | comm -23 - <(cut -f 1 bad.sample | cat - gvcf.07022021.success.txt | sort) | ~/qgg/software/parallel-20200722/bin/parallel -j 96 "bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/ufs18/scratch/huangw53/pigWGS --tmp /mnt/ufs18/scratch/huangw53/tmp --out /mnt/ufs18/scratch/huangw53/swim/gvcfs/ --sample {} > /mnt/ufs18/scratch/huangw53/swim/gvcfs/log/{}.out 2>&1" &

# 07/21/2021: get the ones that have duplicates
# ============================================================

for sample in $(ls log/*.out | sed 's/log\///' | sed 's/\.out//')
do
  if [ `grep amb log/$sample.out | wc -l` -gt 0 ] || [ `grep error log/$sample.out | wc -l` -gt 0 ] || [ `grep "completed haplotype" log/$sample.out | wc -l` -eq 0 ]
  then 
    echo $sample
  fi
done

for sample in $(ls log/*.out | sed 's/log\///' | sed 's/\.out//')
do
  if [ `grep amb log/$sample.out | wc -l` -gt 0 ] || [ `grep error log/$sample.out | wc -l` -gt 0 ] || [ `grep "completed haplotype" log/$sample.out | wc -l` -eq 0 ] || [ `grep $sample bad.sample | wc -l` -gt 0 ]
  then
    echo $sample
    rm -r $sample
    rm log/$sample.out
  fi
done

# get a complete list of successful gvcf runs and their inferred sex

for sample in $(ls log/*.out | sed 's/log\///' | sed 's/\.out//')
do
  if [ `grep amb log/$sample.out | wc -l` -gt 0 ] || [ `grep error log/$sample.out | wc -l` -gt 0 ] || [ `grep "completed haplotype" log/$sample.out | wc -l` -eq 0 ] || [ `grep $sample bad.sample | wc -l` -gt 0 ]
  then
    echo $sample error
  else
    awk '{print "'$sample'\t"$1}' $sample/sex.info >> gvcf.complete.list.07252021
  fi
done



