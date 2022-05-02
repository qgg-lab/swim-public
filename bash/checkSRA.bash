# ==========================
# = check and fix SRA runs =
# ==========================

# the SRA accessions information is present at /mnt/research/qgg/SRAdb/SRA_Accessions.tab.latest
# downloaded using codes in a different project directory sradb/
# where the run and sample information for the "live" samples has been extracted and sorted by sample accession
# SRA_Accessions.tab.latest.sample.sorted
# and also by run accession
# SRA_Accessions.tab.latest.run.sorted
# some run may not have sample association and is indicated by "-"
# ============================================================

# 1. extract run and sample information from the list
# ============================================================

tail -n+2 pig_SRA_RUN-2020-12-29.txt | awk '{print $2"\t"$1}' | sort -k1,1 -k2,2 | ~/qgg/software/bedtools-2.27.1/bin/bedtools groupby -g 1 -c 2 -o collapse > pig.sra.20201229.sample.txt

# 2. check if there are the same runs that appear in multiple samples
# ============================================================

tail -n+2 pig_SRA_RUN-2020-12-29.txt | awk '{print $1}' | sort | uniq -c | awk '$1 > 1'
# for example there are a bunch of runs that belong to two samples SRS965765, and SRS965762

# 3. check if all the list captures all SRA runs for a given SRA sample
# ============================================================

# there are 726 SRA samples in total
join -t $'\t' pig.sra.20101229.sample.txt SRA_Accessions.tab.latest.sample.sorted | cut -f 1 | sort | uniq | wc -l
# get 725, this means one of the samples (SRS2040573_1) is not in the SRA accession
# this is OK this one has been checked to be same sample ID but multiple
# runs belonging to two different animals

join -t $'\t' pig.sra.20201229.sample.txt SRA_Accessions.tab.latest.sample.sorted | sort -k1,1 -k2,2 -k3,3 | ~/qgg/software/bedtools-2.27.1/bin/bedtools groupby -g 1 -c 2,3 -o distinct,collapse | awk '$2 != $3' > pig.sra.20201229.discordant.runs
