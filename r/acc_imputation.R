args <- commandArgs(TRUE)
chr <- as.numeric(args[1])
server <- args[2]

setwd("/mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/imputation_result/012")

#---------- ref maf ----------#
maf=read.table(paste0("maf.",server,".chr", chr,".target.txt"),head=F)
colnames(maf) <- c("maf")
numColObs <- as.numeric(dim(maf)[1])
numColPred <- as.numeric(dim(maf)[1])

#---------- Observed genotypes ----------#
print("loading observed genotypes (012): Validation population")
geno_obs <- matrix(scan(paste0("phasedObs012",".chr", chr,".raw"),what = "character"), nrow = (100+1), ncol = (numColObs + 6), byrow=T)
geno_obs <- geno_obs[-c(1), -c(1:6)]
geno_obs <- t(geno_obs)
dim(geno_obs)
geno_obs <- as.data.frame(geno_obs)

#---------- Predicted genotypes ----------#
print("loading predicted genotypes (012): Imputed population")
geno_pred <- matrix(scan(paste0("phasedPred012.",server,".chr", chr,".raw"),what = "character"), nrow = (100+1), ncol = (numColPred + 6), byrow=T)
geno_pred <- geno_pred[-c(1), -c(1:6)]
geno_pred <- t(geno_pred)
dim(geno_pred)
geno_pred <- as.data.frame(geno_pred)


#---------- Concordance Rate ----------#
compare <- geno_obs == geno_pred

CR_snp <- rowMeans(compare, na.rm = TRUE)
CR_snp <- matrix(CR_snp, ncol=1)
colnames(CR_snp) <- c("concRate")

#---------- Correlation ----------#

geno_obs <- as.matrix(geno_obs)
mode(geno_obs) <- "numeric"

geno_pred <- as.matrix(geno_pred)
mode(geno_pred) <- "numeric"

COR_snp <- sapply(1:nrow(geno_obs), function(i) cor(geno_obs[i,], geno_pred[i,] ,use="pairwise.complete.obs"))
COR_snp <- matrix(COR_snp, ncol=1)
colnames(COR_snp) <- c("r_calc")

#---------- accuracy by maf ----------#
accuracy_maf <- cbind(maf,CR_snp,COR_snp)
accuracy_maf=accuracy_maf[!is.na(accuracy_maf$r_calc),]
dir.create("/mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/imputation_result/summary", recursive = T)

write.table(accuracy_maf,paste0( "/mnt/ls15/scratch/users/dingrong/swim/imputationQCvcf/imputation_result/summary/accuracy_maf.",server,".chr", chr,".txt"),quote = FALSE, row.names=FALSE, col.names=TRUE)
