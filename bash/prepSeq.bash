# ==========================
# = prepare sequence files =
# ==========================

# store fastq files in directories named with unique IDs
# uniform fastq format fastq.gz
# and change all to mode r--r----- = chmod 440

# note that this process cannot be fully automated
# and must be done with a certain level of human supervision
# future additions will be more structured

# in /mnt/ufs18/scratch/huangw53/duroc
# American Durocs S21
# these files are in the format of A12_1_clean.fq.gz
# and A12_2_clean.fq.gz for the two ends for the paired-end
# not sure what the clean means
# each file contains pooled flow cell/lane from multiple runs
# need to separate them into different fastq files and 
# get read information
# ============================================================

for sample in `ls *_1_clean.fq.gz | sed 's/_1_clean\.fq\.gz//'`
do
  mkdir /mnt/ufs18/scratch/huangw53/pigWGS/S21_"$sample"/
  sbatch --export=sample="$sample" --output=/mnt/ufs18/scratch/huangw53/pigWGS/S21_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/S21_"$sample"/fq.split.err splitFastQS21.sbatch
done

# check for success
for sample in `ls *_1_clean.fq.gz | sed 's/_1_clean\.fq\.gz//'`
do
  if [ `wc -l /mnt/ufs18/scratch/huangw53/pigWGS/S21_"$sample"/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "Account" /mnt/ufs18/scratch/huangw53/pigWGS/S21_"$sample"/fq.split.out | wc -l | awk '{print $1}'` -eq 0 ]
  then
    echo "$sample"
    sbatch --export=sample="$sample" --output=/mnt/ufs18/scratch/huangw53/pigWGS/S21_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/S21_"$sample"/fq.split.err splitFastQS21.sbatch
  fi
done

# one sample A29 needs re-run, somehow the execution stopped without outputting error
# done and check to make sure it's OK

# in /mnt/gs18/scratch/users/huangw53/seq/s22/
# Canadian Durocs S22
# some of the files do not have the string "clean" in file names
# ============================================================

for sample in `ls *_1_clean.fq.gz | sed 's/_1_clean\.fq\.gz//'`
do
	mkdir /mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/
	sbatch --export=sample="$sample",add="_clean" --output=/mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.err splitFastQS22.sbatch
done

for sample in `ls *_1.fq.gz | sed 's/_1\.fq\.gz//'`
do
  mkdir /mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/
  sbatch --export=sample="$sample",add="" --output=/mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.err splitFastQS22.sbatch
done

# check for success
for sample in `ls *_1_clean.fq.gz | sed 's/_1_clean\.fq\.gz//'`
do
  if [ `wc -l /mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "Account" /mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.out | wc -l | awk '{print $1}'` -eq 0 ]
  then
    echo "$sample"
    sbatch --export=sample="$sample",add="_clean" --output=/mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.err splitFastQS22.sbatch
	fi
done
# all good!

for sample in `ls *_1.fq.gz | sed 's/_1\.fq\.gz//'`
do
  if [ `wc -l /mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "Account" /mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.out | wc -l | awk '{print $1}'` -eq 0 ]
  then
    echo "$sample"
    sbatch --export=sample="$sample",add="" --output=/mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.err splitFastQS22.sbatch
	fi
done
# all good!

# S22 in /mnt/gs18/scratch/users/huangw53/seq/01205-1.9/X101SC19092770-Z01-J002_1.7T/cleandata
# ============================================================

for sample in `ls *_1_clean.fq.gz | sed 's/_1_clean\.fq\.gz//'`
do
  mkdir /mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/
  sbatch --export=sample="$sample",add="_clean" --output=/mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.err splitFastQS22_2.sbatch
done

for sample in `ls *_1_clean.fq.gz | sed 's/_1_clean\.fq\.gz//'`
do
  if [ `wc -l /mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "Account" /mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.out | wc -l | awk '{print $1}'` -eq 0 ]
  then
    echo "$sample"
    sbatch --export=sample="$sample",add="_clean" --output=/mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.err splitFastQS22_2.sbatch
  fi
done
# all good!

# S22 in /mnt/gs18/scratch/users/huangw53/seq/01315-1.9/X101SC19092770-Z01-J001_2T/cleandata
# ============================================================

for sample in `ls *_1_clean.fq.gz | sed 's/_1_clean\.fq\.gz//'`
do
  mkdir /mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/
  sbatch --export=sample="$sample",add="_clean" --output=/mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.err splitFastQS22_3.sbatch
done

for sample in `ls *_1_clean.fq.gz | sed 's/_1_clean\.fq\.gz//'`
do
  if [ `wc -l /mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "Account" /mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.out | wc -l | awk '{print $1}'` -eq 0 ]
  then
    echo "$sample"
    sbatch --export=sample="$sample",add="_clean" --output=/mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/S22_"$sample"/fq.split.err splitFastQS22_3.sbatch
  fi
done
# S0725, S0741, S0787 are problematic and confirmed fixed after re-run

# yorkshire in /mnt/ufs18/scratch/huangw53/yorkshire
# ============================================================

for sample in `ls *_1_clean.fq.gz | sed 's/_1_clean\.fq\.gz//'`
do
	mkdir /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/
	sbatch --export=sample="$sample",add="_clean" --output=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.err splitFastQyorkshire1.sbatch
done

for sample in `ls *_1.fq.gz | sed 's/_1\.fq\.gz//'`
do
  mkdir /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/
  sbatch --export=sample="$sample",add="" --output=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.err splitFastQyorkshire1.sbatch
done

for sample in `ls *_1_clean.fq.gz | sed 's/_1_clean\.fq\.gz//'`
do
  if [ `wc -l /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "Account" /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.out | wc -l | awk '{print $1}'` -eq 0 ]
  then
    echo "$sample"
		# sbatch --export=sample="$sample",add="_clean" --output=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.err splitFastQyorkshire1.sbatch
  fi
done
# two samples B10, B20 have unexpected end of file errors, don't deal with these for the moment

for sample in `ls *_1.fq.gz | sed 's/_1\.fq\.gz//'`
do
  if [ `wc -l /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "Account" /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.out | wc -l | awk '{print $1}'` -eq 0 ]
  then
    echo "$sample"
    sbatch --export=sample="$sample",add="" --output=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.err splitFastQyorkshire1.sbatch
done
# all good!

# deal with the B10 and B20 samples, both samples have problem with the first read
# find the reads that pair and the remaining reads convert to single end.
gunzip -c B10_1_clean.fq.gz > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/B10_1_clean.fq 2> /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/B10_1_clean.fq.gunzip.log &
gunzip -c B10_2_clean.fq.gz > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/B10_2_clean.fq 2> /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/B10_2_clean.fq.gunzip.log &

gunzip -c B20_1_clean.fq.gz > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B20/B20_1_clean.fq 2> /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B20/B20_1_clean.fq.gunzip.log &
gunzip -c B20_2_clean.fq.gz > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B20/B20_2_clean.fq 2> /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B20/B20_2_clean.fq.gunzip.log &

wc -l /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/B10_1_clean.fq > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/B10_1_clean.wc.out &
wc -l /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/B10_2_clean.fq > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/B10_2_clean.wc.out &
# @A00133:48:H2N53DMXX:1: line 1-265821808
# @A00133:48:H2N53DMXX:2: line 265821809-end
head -n 265821808 /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/B10_1_clean.fq | gzip > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/pe.A00133.48.H2N53DMXX.1_1.fastq.gz &
head -n 265821808 /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/B10_2_clean.fq | gzip > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/pe.A00133.48.H2N53DMXX.1_2.fastq.gz &
head -n 482865236 /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/B10_1_clean.fq | tail -n+265821809 | gzip > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/pe.A00133.48.H2N53DMXX.2_1.fastq.gz &
head -n 482865236 /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/B10_2_clean.fq | tail -n+265821809 | gzip > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/pe.A00133.48.H2N53DMXX.2_2.fastq.gz &
tail -n+482865237 /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/B10_2_clean.fq | gzip > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/se.A00133.48.H2N53DMXX.3_1.fastq.gz &
echo se.A00133.48.H2N53DMXX.2 A00133 48 H2N53DMXX 2 >> /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B10/file.info

wc -l /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B20/B20_1_clean.fq > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B20/B20_1_clean.wc.out &
wc -l /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B20/B20_2_clean.fq > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B20/B20_2_clean.wc.out &
# @A00159:35:H2JHFDMXX:1: line 1-310898772
# @A00159:35:H2JHFDMXX:2: line 310898772-end
# _2 1-268118912
head -n 268118912 /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B20/B20_1_clean.fq | gzip > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B20/pe.A00159.35.H2JHFDMXX.1_1.fastq.gz &
head -n 268118912 /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B20/B20_2_clean.fq | gzip > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B20/pe.A00159.35.H2JHFDMXX.1_2.fastq.gz &
head -n 310898772 /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B20/B20_1_clean.fq | tail -n+268118913 | gzip > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B20/se.A00159.35.H2JHFDMXX.1_1.fastq.gz &
tail -n+310898773 /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B20/B20_1_clean.fq | gzip > /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B20/se.A00159.35.H2JHFDMXX.2_1.fastq.gz &
echo se.A00159.35.H2JHFDMXX.1 A00159 35 H2JHFDMXX 1 >> /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B20/file.info
sed -i 's/pe.A00159.35.H2JHFDMXX.2/se.A00159.35.H2JHFDMXX.2/' /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_B20/file.info

# yorkshire in /mnt/gs18/scratch/users/huangw53/seq/Clean/
# all are paired-end, and all have only 1 pair of files
# ============================================================

for sample in `ls -d */ | sed 's/\///'`
do
  file=`ls "$sample" | grep _1 | sed 's/_1\.fq\.gz//'`
	mkdir /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/
	sbatch --export=sample="$sample",add="",file="$file" --output=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.err splitFastQyorkshire2.sbatch	
done

for sample in `ls -d */ | sed 's/\///'`
do
  if [ `wc -l /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "Account" /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.out | wc -l | awk '{print $1}'` -eq 0 ]
  then
    echo "$sample"
		file=`ls "$sample" | grep _1 | sed 's/_1\.fq\.gz//'`
    sbatch --export=sample="$sample",add="",file="$file" --output=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.err splitFastQyorkshire2.sbatch	
done
# all good!

# yorkshire in /mnt/gs18/scratch/users/huangw53/seq/F18FTSSCKF1561_PIGpkwR_1/Clean
# ============================================================

for sample in `ls -d */ | sed 's/\///'`
do
  file=`ls "$sample" | grep _1 | sed 's/_1\.fq\.gz//'`
	mkdir /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/
	sbatch --export=sample="$sample",add="",file="$file" --output=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.err splitFastQyorkshire3.sbatch	
done

for sample in `ls -d */ | sed 's/\///'`
do
  if [ `wc -l /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "Account" /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.out | wc -l | awk '{print $1}'` -eq 0 ]
  then
    echo "$sample"
		file=`ls "$sample" | grep _1 | sed 's/_1\.fq\.gz//'`
		sbatch --export=sample="$sample",add="",file="$file" --output=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.err splitFastQyorkshire3.sbatch	
done
# all good!

# yorkshire in /mnt/gs18/scratch/users/huangw53/seq/F18FTSSCKF1561_PIGpkwR_3/Clean
# ============================================================

for sample in `ls -d */ | sed 's/\///'`
do
  file=`ls "$sample" | grep _1 | sed 's/_1\.fq\.gz//'`
	mkdir /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/
	sbatch --export=sample="$sample",add="",file="$file" --output=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.err splitFastQyorkshire4.sbatch	
done

for sample in `ls -d */ | sed 's/\///'`
do
  if [ `wc -l /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "Account" /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.out | wc -l | awk '{print $1}'` -eq 0 ]
  then
    echo "$sample"
		file=`ls "$sample" | grep _1 | sed 's/_1\.fq\.gz//'`
		sbatch --export=sample="$sample",add="",file="$file" --output=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.err splitFastQyorkshire4.sbatch	
done
# all good!

# yorkshire newly shared in /mnt/research/qgg/share/
# ============================================================

for sample in B60 B61 B62 B63 B64 B65
do
  file=$sample
	mkdir /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/
	sbatch --export=sample="$sample",add="_clean",file="$file" --output=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.err splitFastQyorkshire5.sbatch	
done

for sample in B60 B61 B62 B63 B64 B65
do
  if [ `wc -l /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "Account" /mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.out | wc -l | awk '{print $1}'` -eq 0 ]
  then
    echo "$sample"
    file=$sample
  	sbatch --export=sample="$sample",add="_clean",file="$file" --output=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.out --error=/mnt/ufs18/scratch/huangw53/pigWGS/yorkshire_"$sample"/fq.split.err splitFastQyorkshire5.sbatch	
done
# all good!

# get sequences from SRA
# ============================================================

# get SRR for SRSs
# ============================================================

sort /mnt/ufs18/scratch/huangw53/pigWGS/sra.srs.01212021 | join -t $'\t' - /mnt/research/qgg/resource/SRAdb/SRA_Accessions.tab.latest.sample.sorted > /mnt/ufs18/scratch/huangw53/pigWGS/sra.srs.srr.01212021

# prep SRA sequence
# ============================================================

cat /mnt/ufs18/scratch/huangw53/pigWGS/sra.srs.srr.01212021 | ~/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2 -o collapse | tr '\t' '\n' | ~/qgg/software/parallel-20200722/bin/parallel -N 2 -j 16 "bash /mnt/research/qgg/resource/swim/sra2fqwait.bash --dir /mnt/ufs18/scratch/huangw53/pigWGS/ --sample {1} --srr {2} --resource /mnt/research/qgg/resource/swim/ > /mnt/ufs18/scratch/huangw53/pigWGS/log/sra.{1}.log 2>&1" > /mnt/ufs18/scratch/huangw53/pigWGS/log/sra10022020.log 2>&1 &

# many has not completed, due to slurm problem
# ============================================================

grep completed log/sra.*.log | awk -F ":" '{print $1}' | sed 's/log\/sra\.//' | sed 's/\.log//' | sort | uniq > sra02242021success.txt

cat /mnt/ufs18/scratch/huangw53/pigWGS/sra.srs.srr.01212021 | ~/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2 -o collapse | join -a 1 -v 1 -t $'\t' - sra02242021success.txt | cut -f 1 | xargs rm -r
cat /mnt/ufs18/scratch/huangw53/pigWGS/sra.srs.srr.01212021 | ~/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2 -o collapse | join -a 1 -v 1 -t $'\t' - sra02242021success.txt | tr '\t' '\n' | ~/qgg/software/parallel-20200722/bin/parallel -N 2 -j 24 "bash /mnt/research/qgg/resource/swim/sra2fqwait.bash --dir /mnt/ufs18/scratch/huangw53/pigWGS/ --sample {1} --srr {2} --resource /mnt/research/qgg/resource/swim/ > /mnt/ufs18/scratch/huangw53/pigWGS/log/sra.{1}.log 2>&1" > /mnt/ufs18/scratch/huangw53/pigWGS/log/sra02242021.log 2>&1 &

# check again
# ============================================================

grep completed log/sra.*.log | awk -F ":" '{print $1}' | sed 's/log\/sra\.//' | sed 's/\.log//' | sort | uniq > sra02272021success.txt
cat /mnt/ufs18/scratch/huangw53/pigWGS/sra.srs.srr.01212021 | ~/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2 -o collapse | join -a 1 -v 1 -t $'\t' - sra02272021success.txt | cut -f 1 | xargs rm -r
cat /mnt/ufs18/scratch/huangw53/pigWGS/sra.srs.srr.01212021 | ~/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2 -o collapse | join -a 1 -v 1 -t $'\t' - sra02272021success.txt | tr '\t' '\n' | ~/qgg/software/parallel-20200722/bin/parallel -N 2 -j 24 "bash /mnt/research/qgg/resource/swim/sra2fqwait.bash --dir /mnt/ufs18/scratch/huangw53/pigWGS/ --sample {1} --srr {2} --resource /mnt/research/qgg/resource/swim/ > /mnt/ufs18/scratch/huangw53/pigWGS/log/sra.{1}.log 2>&1" > /mnt/ufs18/scratch/huangw53/pigWGS/log/sra02272021.log 2>&1 &

# check another time
# ============================================================

grep completed log/sra.*.log | awk -F ":" '{print $1}' | sed 's/log\/sra\.//' | sed 's/\.log//' | sort | uniq > sra03012021success.txt
cat /mnt/ufs18/scratch/huangw53/pigWGS/sra.srs.srr.01212021 | ~/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2 -o collapse | join -a 1 -v 1 -t $'\t' - sra03012021success.txt | cut -f 1 | xargs rm -r
# check all of those and they all have problem prefetch
# one did not finish, restart this one "SRS703324"
bash /mnt/research/qgg/resource/swim/sra2fqwait.bash --dir /mnt/ufs18/scratch/huangw53/pigWGS/ --sample SRS703324 --srr SRR1581126,SRR1581127,SRR1581128 --resource /mnt/research/qgg/resource/swim/ > /mnt/ufs18/scratch/huangw53/pigWGS/log/sra.SRS703324.log 2>&1 &

# 03072021: add two Iberian pigs to the SRA download list
# ============================================================

echo -e "SRS875335\nSRS2170012" > sra.srs.03072021
sort /mnt/ufs18/scratch/huangw53/pigWGS/sra.srs.03072021 | join -t $'\t' - /mnt/research/qgg/resource/SRAdb/SRA_Accessions.tab.latest.sample.sorted > /mnt/ufs18/scratch/huangw53/pigWGS/sra.srs.srr.03072021

cat /mnt/ufs18/scratch/huangw53/pigWGS/sra.srs.srr.03072021 | ~/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2 -o collapse | tr '\t' '\n' | ~/qgg/software/parallel-20200722/bin/parallel -N 2 -j 16 "bash /mnt/research/qgg/resource/swim/sra2fqwait.bash --dir /mnt/ufs18/scratch/huangw53/pigWGS/ --sample {1} --srr {2} --resource /mnt/research/qgg/resource/swim/ > /mnt/ufs18/scratch/huangw53/pigWGS/log/sra.{1}.log 2>&1" > /mnt/ufs18/scratch/huangw53/pigWGS/log/sra03072021.log 2>&1 &

# 03092021: get all SRA success sequences
# ============================================================

grep completed log/sra.*.log | awk -F ":" '{print $1}' | sed 's/log\/sra\.//' | sed 's/\.log//' | sort | uniq > sra03092021success.txt

# 04/10/2021: a set of new sequences from SCAU
# unknow, unknow2, unknow3, unknow4
# store sequences in /mnt/gs18/scratch/users/longnany/wen/pigWGS2
# ============================================================

# unknow
for sample in `ls /mnt/gs18/scratch/users/huangw53/newseq/unknow/*_1_clean.fq.gz | sed 's/_1_clean\.fq\.gz//' | sed 's/\/mnt\/gs18\/scratch\/users\/huangw53\/newseq\/unknow\///'`
do
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
  file1="/mnt/gs18/scratch/users/huangw53/newseq/unknow/"$sample"_1_clean.fq.gz"
  file2="/mnt/gs18/scratch/users/huangw53/newseq/unknow/"$sample"_2_clean.fq.gz"
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1="$file1",file2="$file2" --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQunknow.sbatch
done

for sample in `ls /mnt/gs18/scratch/users/huangw53/newseq/unknow/*_1_clean.fq.gz | sed 's/_1_clean\.fq\.gz//' | sed 's/\/mnt\/gs18\/scratch\/users\/huangw53\/newseq\/unknow\///'`
do
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
	then
    echo $sample
	fi
done
# all good

# unknow2
for sample in `ls /mnt/gs18/scratch/users/huangw53/newseq/unknow2/*_1_clean.fq.gz | sed 's/_1_clean\.fq\.gz//' | sed 's/\/mnt\/gs18\/scratch\/users\/huangw53\/newseq\/unknow2\///'`
do
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
  file1="/mnt/gs18/scratch/users/huangw53/newseq/unknow2/"$sample"_1_clean.fq.gz"
  file2="/mnt/gs18/scratch/users/huangw53/newseq/unknow2/"$sample"_2_clean.fq.gz"
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1="$file1",file2="$file2" --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQunknow.sbatch
done

for sample in `ls /mnt/gs18/scratch/users/huangw53/newseq/unknow2/*_1_clean.fq.gz | sed 's/_1_clean\.fq\.gz//' | sed 's/\/mnt\/gs18\/scratch\/users\/huangw53\/newseq\/unknow2\///'`
do
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
  then
    echo $sample
  fi
done
# all good

# unknow3
for sample in `ls /mnt/gs18/scratch/users/huangw53/newseq/unknow3/*_1_clean.fq.gz | sed 's/_1_clean\.fq\.gz//' | sed 's/\/mnt\/gs18\/scratch\/users\/huangw53\/newseq\/unknow3\///'`
do
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
  file1="/mnt/gs18/scratch/users/huangw53/newseq/unknow3/"$sample"_1_clean.fq.gz"
  file2="/mnt/gs18/scratch/users/huangw53/newseq/unknow3/"$sample"_2_clean.fq.gz"
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1="$file1",file2="$file2" --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQunknow.sbatch
done

for sample in `ls /mnt/gs18/scratch/users/huangw53/newseq/unknow3/*_1_clean.fq.gz | sed 's/_1_clean\.fq\.gz//' | sed 's/\/mnt\/gs18\/scratch\/users\/huangw53\/newseq\/unknow3\///'`
do
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
	then
    echo $sample
	fi
done
# all good
# there are four samples that have been processed in unknow and unknow2
# those have larger numbers of reads

# unknow4
for sample in `ls /mnt/gs18/scratch/users/huangw53/newseq/unknow4/*_1_clean.fq.gz | sed 's/_1_clean\.fq\.gz//' | sed 's/\/mnt\/gs18\/scratch\/users\/huangw53\/newseq\/unknow4\///'`
do
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
  file1="/mnt/gs18/scratch/users/huangw53/newseq/unknow4/"$sample"_1_clean.fq.gz"
  file2="/mnt/gs18/scratch/users/huangw53/newseq/unknow4/"$sample"_2_clean.fq.gz"
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1="$file1",file2="$file2" --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQunknow.sbatch
done
# two samples have been processed, checked to make sure that the previous ones are larger
for sample in `ls /mnt/gs18/scratch/users/huangw53/newseq/unknow4/*_1_clean.fq.gz | sed 's/_1_clean\.fq\.gz//' | sed 's/\/mnt\/gs18\/scratch\/users\/huangw53\/newseq\/unknow4\///'`
do
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
	then
    echo $sample
	fi
done
# all good

# 04/23/2021 Korean pigs
# in /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/D10,DK,KL,KNPs,L5,LK
# ============================================================

for file1 in `cat <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/D10/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/DK/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/KL/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/KNPs/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/L5/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/LK/ -name "*_1.fastq.gz")`
do
  sample=korea_`echo $file1 | sed 's/.*\///' | sed 's/_filtered_1\.fastq\.gz//'`
  file2=`echo $file1 | sed 's/filtered_1\.fastq\.gz/filtered_2\.fastq\.gz/'`
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
	njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
	while [ $njobs -ge 16 ]
	do
    sleep 2m
    njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  done
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1="$file1",file2="$file2" --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQunknow.sbatch
  sleep 1m
done > /mnt/gs18/scratch/users/longnany/wen/pigWGS2/log/04232021.split.korea.fastq.log 2>&1 &

for file1 in `cat <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/D10/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/DK/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/KL/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/KNPs/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/L5/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/LK/ -name "*_1.fastq.gz")`
do
  sample=korea_`echo $file1 | sed 's/.*\///' | sed 's/_filtered_1\.fastq\.gz//'`
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
  then
    echo $sample
    rm -r /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
    file2=`echo $file1 | sed 's/filtered_1\.fastq\.gz/filtered_2\.fastq\.gz/'`
    mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
    sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1="$file1",file2="$file2" --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQunknow.sbatch
  fi
done > /mnt/gs18/scratch/users/longnany/wen/pigWGS2/log/04292021.split.korea.fastq.log 2>&1 &

# still three did not finish, increase time
for file1 in `cat <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/D10/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/DK/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/KL/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/KNPs/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/L5/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/LK/ -name "*_1.fastq.gz")`
do
  sample=korea_`echo $file1 | sed 's/.*\///' | sed 's/_filtered_1\.fastq\.gz//'`
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
  then
    echo $sample
    rm -r /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
    file2=`echo $file1 | sed 's/filtered_1\.fastq\.gz/filtered_2\.fastq\.gz/'`
    mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
    sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1="$file1",file2="$file2" --time=24:00:00 --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQunknow.sbatch
  fi
done > /mnt/gs18/scratch/users/longnany/wen/pigWGS2/log/05032021.split.korea.fastq.log 2>&1 &

for file1 in `cat <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/D10/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/DK/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/KL/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/KNPs/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/L5/ -name "*_1.fastq.gz") <(find /mnt/ufs18/scratch/pelicion/rawData/pigKoreanSequence/LK/ -name "*_1.fastq.gz")`
do
  sample=korea_`echo $file1 | sed 's/.*\///' | sed 's/_filtered_1\.fastq\.gz//'`
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
  then
    echo $sample
  fi
done
# all good!

# 05032021: check SRA samples, add new samples
# ============================================================

sort /mnt/ufs18/scratch/huangw53/pigWGS/sra.srs.05032021 | join -t $'\t' - /mnt/research/qgg/resource/SRAdb/SRA_Accessions.tab.latest.sample.sorted > /mnt/ufs18/scratch/huangw53/pigWGS/sra.srs.srr.05032021

cat /mnt/ufs18/scratch/huangw53/pigWGS/sra.srs.srr.05032021 | ~/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2 -o collapse | tr '\t' '\n' | ~/qgg/software/parallel-20200722/bin/parallel -N 2 -j 80 "bash /mnt/research/qgg/resource/swim/sra2fqwait.bash --dir /mnt/ufs18/scratch/huangw53/pigWGS/ --sample {1} --srr {2} --resource /mnt/research/qgg/resource/swim/ > /mnt/ufs18/scratch/huangw53/pigWGS/log/sra.{1}.log 2>&1" > /mnt/ufs18/scratch/huangw53/pigWGS/log/sra05032021.log 2>&1 &

grep completed log/sra.*.log | awk -F ":" '{print $1}' | sed 's/log\/sra\.//' | sed 's/\.log//' | sort | uniq > sra05062021success.txt
cat /mnt/ufs18/scratch/huangw53/pigWGS/sra.srs.srr.05032021 | ~/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2 -o collapse | join -a 1 -v 1 -t $'\t' - sra05062021success.txt
# no unsuccessful downloads

# 05032021: more scau sequences
# there is one animal that has been genotyped by both illumina and bgi
# the illumina sequences
# ============================================================

mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/CC3_C5_illumina/
sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/CC3_C5_illumina",file1="/mnt/home/huangw53/scratch/newseq/illumina/CC3_C5_FDSW202560868-1r_1.clean.fq.gz",file2="/mnt/home/huangw53/scratch/newseq/illumina/CC3_C5_FDSW202560868-1r_2.clean.fq.gz" --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/CC3_C5_illumina/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/CC3_C5_illumina/fq.split.err splitFastQunknow.sbatch

# 05142021: 36 BGI files in directories
# ============================================================

for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/CC3etc/`
do
  sample=`echo $sdir | awk -F "_" '{print $1}'`
  file1=/mnt/gs18/scratch/users/huangw53/newseq/CC3etc/$sdir/$sdir"_1.clean.fq.gz"
  file2=/mnt/gs18/scratch/users/huangw53/newseq/CC3etc/$sdir/$sdir"_2.clean.fq.gz"
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
  njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  while [ $njobs -ge 16 ]
  do
    sleep 2m
    njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  done
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1=$file1,file2=$file2,tmp=/mnt/gs18/scratch/users/longnany/wen/tmp,sample=$sample --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQBGI.sbatch
  sleep 1m
done &

# many failed because of std disk quota issue, due to pipe I think
# and many jobs running on the same node
for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/CC3etc/`
do
  sample=`echo $sdir | awk -F "_" '{print $1}'`
  file1=/mnt/gs18/scratch/users/huangw53/newseq/CC3etc/$sdir/$sdir"_1.clean.fq.gz"
  file2=/mnt/gs18/scratch/users/huangw53/newseq/CC3etc/$sdir/$sdir"_2.clean.fq.gz"
  njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
  then
    echo $sample
    rm -rf /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
    mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
    while [ $njobs -ge 8 ]
    do
      sleep 5m
      njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
    done
    sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1=$file1,file2=$file2,tmp=/mnt/gs18/scratch/users/longnany/wen/tmp,sample=$sample --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQBGI.sbatch
    sleep 5m
  fi
done &

for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/CC3etc/`
do
  sample=`echo $sdir | awk -F "_" '{print $1}'`
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
  then
    echo $sample
  fi
done
# all good

# clean_0120
for file1 in `ls /mnt/gs18/scratch/users/huangw53/newseq/clean_0120/*_1.clean.fq.gz`
do
  sample=`echo $file1 | sed 's/.*clean_0120\///' | sed 's/\_1\.clean\.fq\.gz//'`
  file2=`echo $file1 | sed 's/1\.clean\.fq\.gz/2\.clean\.fq\.gz/'`
  njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
  while [ $njobs -ge 16 ]
  do
    sleep 1m
    njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  done
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1=$file1,file2=$file2,tmp=/mnt/gs18/scratch/users/longnany/wen/tmp,sample=$sample --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQBGI.sbatch
  sleep 1m
done &

for file1 in `ls /mnt/gs18/scratch/users/huangw53/newseq/clean_0120/*_1.clean.fq.gz`
do
  sample=`echo $file1 | sed 's/.*clean_0120\///' | sed 's/\_1\.clean\.fq\.gz//'`
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
  then
    echo $sample
  fi
done
# all good

# clean_0120/2.cleandata in directory
for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/clean_0120/2.cleandata/`
do
  sample="$sdir"
  file1=/mnt/gs18/scratch/users/huangw53/newseq/clean_0120/2.cleandata/$sdir/$sdir"_1.clean.fq.gz"
  file2=/mnt/gs18/scratch/users/huangw53/newseq/clean_0120/2.cleandata/$sdir/$sdir"_2.clean.fq.gz"
  njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
  while [ $njobs -ge 16 ]
  do
    sleep 1m
    njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  done
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1=$file1,file2=$file2,tmp=/mnt/gs18/scratch/users/longnany/wen/tmp,sample=$sample --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQBGI.sbatch
  sleep 1m
done &

for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/clean_0120/2.cleandata/`
do
  sample="$sdir"
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
  then
    echo $sample
  fi
done 

# 05/19/2021: Landrace_82
# ============================================================

for file1 in `ls /mnt/gs18/scratch/users/huangw53/newseq/Landrace_82/*_1_clean.fq.gz`
do
  sample=`echo $file1 | sed 's/.*Landrace_82\///' | sed 's/\_1_clean\.fq\.gz//'`
  file2=`echo $file1 | sed 's/1_clean\.fq\.gz/2_clean\.fq\.gz/'`
  njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
  while [ $njobs -ge 16 ]
  do
    sleep 1m
    njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  done
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1="$file1",file2="$file2" --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQunknow.sbatch
  sleep 1m
done &

for file1 in `ls /mnt/gs18/scratch/users/huangw53/newseq/Landrace_82/*_1_clean.fq.gz`
do
  sample=`echo $file1 | sed 's/.*Landrace_82\///' | sed 's/\_1_clean\.fq\.gz//'`
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
  then
    echo $sample
  fi
done

# 05/20/2021: more
# ============================================================

for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/more/`
do
  sample=`echo $sdir | awk -F "_" '{print $1}'`
  file1=/mnt/gs18/scratch/users/huangw53/newseq/more/$sdir/$sdir"_1.clean.fq.gz"
  file2=/mnt/gs18/scratch/users/huangw53/newseq/more/$sdir/$sdir"_2.clean.fq.gz"
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
  njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  while [ $njobs -ge 16 ]
  do
    sleep 1m
    njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  done
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1=$file1,file2=$file2,tmp=/mnt/gs18/scratch/users/longnany/wen/tmp,sample=$sample --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQBGI.sbatch
  sleep 1m
done &

for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/more/`
do
  sample=`echo $sdir | awk -F "_" '{print $1}'`
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
  then
    echo $sample
  fi
done

# 05/20/2021: songyingshen0104035
# ============================================================

for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/songyingshen0104035/1.cleandata/`
do
  sample=`echo $sdir`
file1=/mnt/gs18/scratch/users/huangw53/newseq/songyingshen0104035/1.cleandata/$sdir/$sdir"_1.clean.fq.gz"
file2=/mnt/gs18/scratch/users/huangw53/newseq/songyingshen0104035/1.cleandata/$sdir/$sdir"_2.clean.fq.gz"
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
  njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  while [ $njobs -ge 16 ]
  do
    sleep 1m
    njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  done
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1=$file1,file2=$file2,tmp=/mnt/ufs18/scratch/huangw53/tmp,sample=$sample --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQBGI.sbatch
  sleep 1m
done > log/05022021.songyingshen0104035.fastq.log 2>&1 &

for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/songyingshen0104035/1.cleandata/`
do
  sample=`echo $sdir`
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
  then
    echo $sample
  fi
done

# 05/21/2021: songyingshen1229030
# ============================================================

for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/songyingshen1229030/2.cleandata/`
do
  sample=`echo $sdir | sed 's/_.*//'`
  file1=`ls /mnt/gs18/scratch/users/huangw53/newseq/songyingshen1229030/2.cleandata/$sdir/*.fq.gz | head -n 1`
  file2=`echo $file1 | sed 's/1\.clean\.fq\.gz/2\.clean\.fq\.gz/'`
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
  njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  while [ $njobs -ge 16 ]
  do
    sleep 1m
    njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  done
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1=$file1,file2=$file2,tmp=/mnt/ufs18/scratch/huangw53/tmp,sample=$sample --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQBGI.sbatch
  sleep 1m
done > log/05212021.songyingshen0104035.fastq.log 2>&1 &

for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/songyingshen1229030/2.cleandata/`
do
  sample=`echo $sdir | sed 's/_.*//'`
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
  then
    echo $sample
  fi
done

# 05/22/2021: songyingshen1229033
# ============================================================

for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/songyingshen1229033/2.cleandata/`
do
  sample=`echo $sdir`
  file1="/mnt/gs18/scratch/users/huangw53/newseq/songyingshen1229033/2.cleandata/$sdir/$sdir"_1.clean.fq.gz
  file2="/mnt/gs18/scratch/users/huangw53/newseq/songyingshen1229033/2.cleandata/$sdir/$sdir"_2.clean.fq.gz
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
  njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  while [ $njobs -ge 16 ]
  do
    sleep 1m
    njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  done
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1=$file1,file2=$file2,tmp=/mnt/ufs18/scratch/huangw53/tmp,sample=$sample --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQBGI.sbatch
  sleep 1m
done > log/05222021.songyingshen1229033.fastq.log 2>&1 &

for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/songyingshen1229033/2.cleandata/`
do
  sample=`echo $sdir`
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
  then
    echo $sample
  fi
done

# 05/23/2021: western_pig
# ============================================================

for file1 in `ls /mnt/gs18/scratch/users/huangw53/newseq/western_pig/*_1_clean.fq.gz`
do
  sample=`echo $file1 | sed 's/.*western_pig\///' | awk -F "_" '{print $2}'`
  file2=`echo $file1 | sed 's/1_clean\.fq\.gz/2_clean\.fq\.gz/'`
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
  njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  while [ $njobs -ge 16 ]
  do
    sleep 1m
    njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  done
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1="$file1",file2="$file2" --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQunknow.sbatch
  sleep 1m
done > log/05232021.western_pig.fastq.log 2>&1 &

for file1 in `ls /mnt/gs18/scratch/users/huangw53/newseq/western_pig/*_1_clean.fq.gz`
do
  sample=`echo $file1 | sed 's/.*western_pig\///' | awk -F "_" '{print $2}'`
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
  then
    echo $sample
  fi
done

# 05/23/2021: western_pig_B10_B20_clean_data
# ============================================================

for file1 in `ls /mnt/gs18/scratch/users/huangw53/newseq/western_pig_B10_B20_clean_data/*_1_clean.fq.gz`
do
  sample=`echo $file1 | sed 's/.*western_pig_B10_B20_clean_data\///' | sed 's/\_1_clean\.fq\.gz//'`
  file2=`echo $file1 | sed 's/1_clean\.fq\.gz/2_clean\.fq\.gz/'`
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1="$file1",file2="$file2" --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQunknow.sbatch
  sleep 1m
done > log/05232021.western_pig_B10_B20.fastq.log 2>&1 &
  
# 05/24/2021: western_pig_MGISEQ_2000
# ============================================================

while read line
do
  sample=`echo $line | awk '{print $1}'`
  file1=`echo $line | awk '{print $2}'`
  file2=`echo $file1 | sed 's/_1\.fq\.gz/_2\.fq\.gz/g'`
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
  njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  while [ $njobs -ge 16 ]
  do
    sleep 1m
    njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  done
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1=$file1,file2=$file2,tmp=/mnt/ufs18/scratch/huangw53/tmp,sample=$sample --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQBGI.sbatch
  sleep 1m
done < <(ls /mnt/gs18/scratch/users/huangw53/newseq/western_pig_MGISEQ_2000/*_1.fq.gz | awk '{print $1"\t"$1}' | sed 's/\/mnt\/gs18\/scratch\/users\/huangw53\/newseq\/western_pig_MGISEQ_2000\///' | cut -d "_" -f 2,3- | sed 's/_1\.fq\.gz//' | sed 's/[ab]\t/\t/' | sort -k1,1 | ~/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2 -o collapse -delim "|") > log/05242021.western_pig_B10_B20.fastq.log 2>&1 &

while read line
do
  sample=`echo $line | awk '{print $1}'`
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
  then
    echo $sample
  fi
done < <(ls /mnt/gs18/scratch/users/huangw53/newseq/western_pig_MGISEQ_2000/*_1.fq.gz | awk '{print $1"\t"$1}' | sed 's/\/mnt\/gs18\/scratch\/users\/huangw53\/newseq\/western_pig_MGISEQ_2000\///' | cut -d "_" -f 2,3- | sed 's/_1\.fq\.gz//' | sed 's/[ab]\t/\t/' | sort -k1,1 | ~/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2 -o collapse -delim "|")

# 05/25/2021: western_pig_songyingshen
# ============================================================

for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/western_pig_songyingshen`
do
  sample=$sdir
  file1=`ls /mnt/gs18/scratch/users/huangw53/newseq/western_pig_songyingshen/$sdir/*_1.clean.fq.gz`
  file2=`echo $file1 | sed 's/_1\.clean\.fq\.gz/_2\.clean\.fq\.gz/'`
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
  njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  while [ $njobs -ge 16 ]
  do
    sleep 1m
    njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  done
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1=$file1,file2=$file2,tmp=/mnt/ufs18/scratch/huangw53/tmp,sample=$sample --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQBGI.sbatch
  sleep 1m
done > log/05252021.western_pig_songyingshen.fastq.log &

for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/western_pig_songyingshen`
do
  sample=$sdir
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
  then
    echo $sample
  fi
done

# 05/26/2021: western_pig_songyingshen1216014
# ============================================================

for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/western_pig_songyingshen1216014/00.CleanData`
do
  sample=$sdir
  file1=`ls /mnt/gs18/scratch/users/huangw53/newseq/western_pig_songyingshen1216014/00.CleanData/$sdir/*_1.clean.fq.gz`
  file2=`echo $file1 | sed 's/_1\.clean\.fq\.gz/_2\.clean\.fq\.gz/'`
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
  njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  while [ $njobs -ge 24 ]
  do
    sleep 1m
    njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  done
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1=$file1,file2=$file2,tmp=/mnt/ufs18/scratch/huangw53/tmp,sample=$sample --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQBGI.sbatch
  sleep 1m
done > log/05262021.western_pig_songyingshen1216014.fastq.log &

for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/western_pig_songyingshen1216014/00.CleanData`
do
  sample=$sdir
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
  then
    echo $sample
  fi
done

# 05/27/2021: western_pig_songyingshen1216015
# ============================================================

for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/western_pig_songyingshen1216015/X101SC20101498-Z01-J012/00.CleanData`
do
  sample=$sdir
  file1=`ls /mnt/gs18/scratch/users/huangw53/newseq/western_pig_songyingshen1216015/X101SC20101498-Z01-J012/00.CleanData/$sdir/*_1.clean.fq.gz`
  file2=`echo $file1 | sed 's/_1\.clean\.fq\.gz/_2\.clean\.fq\.gz/'`
  mkdir /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample
  njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  while [ $njobs -ge 24 ]
  do
    sleep 1m
    njobs=$(squeue -u huangw53 -O state | egrep '(RUNNING|PENDING)' | wc -l)
  done
  sbatch --export=dir="/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample",file1=$file1,file2=$file2,tmp=/mnt/ufs18/scratch/huangw53/tmp,sample=$sample --output=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out --error=/mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err splitFastQBGI.sbatch
  sleep 1m
done > log/05272021.western_pig_songyingshen1216015.fastq.log &

for sdir in `ls /mnt/gs18/scratch/users/huangw53/newseq/western_pig_songyingshen1216015/X101SC20101498-Z01-J012/00.CleanData`
do
  sample=$sdir
  if [ `wc -l /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.err | awk '{print $1}'` -ne 0 -o `grep "done.main.process" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -eq 0  -o `grep "error" /mnt/gs18/scratch/users/longnany/wen/pigWGS2/$sample/fq.split.out | wc -l | awk '{print $1}'` -gt 0 ]
  then
    echo $sample
  fi
done

