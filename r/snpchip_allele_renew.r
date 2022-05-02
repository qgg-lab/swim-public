args <- commandArgs(TRUE)

mapping_info <- read.table(args[1],head =F )
bimfile <- read.table(args[2],head =F )
outname <- as.character(args[3])

colnames(mapping_info) <- c("CHR", "POS", "SNP_ID", "REF", "ALT", "QUAL1", "QUAL2", "SNPCHIP")
colnames(bimfile) <- c("chr", "SNP_ID", "map", "pos", "alleleA", "alleleB")
info1=merge(mapping_info,bimfile,by = "SNP_ID", all.x = TRUE)

base.tr <- c("A", "T", "C", "G")
names(base.tr) <- c("T", "A", "G", "C")

alleleA_new <- apply(cbind(as.character(info1$REF), as.character(info1$alleleA), as.character(info1$alleleB)),
                 1, function(x) { if (length(setdiff(x[2:3], x[1])) == 2) { x[2:3] <- base.tr[x[2:3]]; }
					y <- x[2];
						return(y) })
					
alleleB_new <- apply(cbind(as.character(info1$REF), as.character(info1$alleleA), as.character(info1$alleleB)),
                 1, function(x) { if (length(setdiff(x[2:3], x[1])) == 2) { x[2:3] <- base.tr[x[2:3]]; }
				 y<- x[3]; 
				 return(y) })
				 
info1$alleleA_new <- alleleA_new 
info1$alleleB_new <- alleleB_new 
info2=subset(info1,select=c(SNP_ID, alleleA, alleleB, alleleA_new, alleleB_new))

write.table(info2,paste(outname,"_UpdateAlleles.txt",sep=""), row.names = FALSE, col.names =FALSE, quote =FALSE)
