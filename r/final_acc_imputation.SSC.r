args <- commandArgs(TRUE)
server <- args[1]

final_chr_acc = data.frame(SSC=1:18,N=1:18,CR=1:18,r=1:18)
for(chr in 1:18){
chr_acc =read.table(paste0( "/mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/imputation_result/summary/accuracy_maf.",server,".chr", chr,".txt"),head=T,stringsAsFactors = F)
final_chr_acc[chr,"CR"] = mean(chr_acc$concRate,na.rm=T)
final_chr_acc[chr,"r"] = mean(chr_acc$r_calc,na.rm=T)
final_chr_acc[chr,"N"] = nrow(chr_acc)
}

write.table(final_chr_acc,paste0( "/mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/imputation_result/summary/accuracy_SSC.",server,".txt"),quote = FALSE, row.names=FALSE, col.names=TRUE)

