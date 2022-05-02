# get database files
# gff and fa
# modify config
# add the following lines
# data.dir = /mnt/research/qgg/software/snpEff-5.0e/data/
# ssc11.1.genome : Sus scrofa
# ============================================================

cp /mnt/research/qgg/resource/sscrofa11.1/genome/ /mnt/research/qgg/software/snpEff-5.0e/data/ssc11.1/genes.gff
cp /mnt/research/qgg/resource/sscrofa11.1/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa /mnt/research/qgg/software/snpEff-5.0e/data/ssc11.1/sequences.fa

java -jar /mnt/research/qgg/software/snpEff-5.0e/snpEff.jar build -gff3 -v ssc11.1 > /mnt/research/qgg/software/snpEff-5.0e/data/ssc11.1/build.log 2>&1 &
