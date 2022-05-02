library(data.table)
library(dplyr)

args <- commandArgs(TRUE)
inputname <- args[1]
outputname <- args[2]

maf_result=read.table(paste0(inputname,".txt"),head=F)
colnames(maf_result)=c("maf","concRate","r_calc")

maf_result$mafBin <- cut(maf_result$maf,
                         seq(0, 0.5, 0.02),
                         include.lowest=TRUE, right=T,
                         labels=seq(0.02, 0.5, 0.02)
)

maf_result_mean <- maf_result %>% group_by(mafBin) %>% summarise(N=n(),mean.CR=mean(concRate,na.rm=T) ,mean.r=mean(r_calc,na.rm=T)) %>% as.data.frame()

write.table(maf_result_mean,paste0(outputname, ".txt"),quote = FALSE, row.names=FALSE, col.names=TRUE)

