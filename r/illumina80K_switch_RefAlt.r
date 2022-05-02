info1=read.table("80K_info.txt",head=T)
info2=read.table("pig_80k_61371_info2.txt",head=T)
info3=merge(info2,info1,by = "Name", all.x = TRUE)

info3$REF <- toupper(info3$REF)
snp.alleles <- matrix(unlist(strsplit(gsub("\\]", "", gsub("\\[", "", info3$SNP)), split = "/")),
                      byrow = 2, ncol = 2)
base.tr <- c("A", "T", "C", "G")
names(base.tr) <- c("T", "A", "G", "C")
snp.alt <- apply(cbind(snp.alleles, as.character(info3$Strand), as.character(info3$REF)),
                 1, function(x) { if (x[3] == "-") { x[1:2] <- base.tr[x[1:2]]; }
                                                     y <- setdiff(x[1:2], x[4]);
                                                     y <- ifelse(length(y) == 1, y, NA);
                                                     return(y) })
 
info3$alt <- snp.alt                                             

info3$QUAL="."
info3$INFO="PorcineSNP80K"

info4= info3[!is.na(info3$alt),]
info5=subset(info4,select=c(chr,pos,Name,REF,alt,QUAL,QUAL,INFO))
nrow(info5)
write.table(info5,"PorcineSNP80K_60818.info", row.names = FALSE, col.names =FALSE, quote =FALSE)
