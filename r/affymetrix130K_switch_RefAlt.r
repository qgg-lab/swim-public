info1=read.table("130K_info.txt",head=T)
info1= info1[!duplicated(info1$Name),]
info2_1=read.table("pig_130k_A111753_info2.txt",head=T)
info2_2=read.table("pig_130k_B5401_info2.txt",head=T)
info2 = rbind(info2_1,info2_2)
info3=merge(info2,info1,by = "Name", all.x = TRUE)

info3$REF <- toupper(info3$REF)

base.tr <- c("A", "T", "C", "G")
names(base.tr) <- c("T", "A", "G", "C")
snp.alt <- apply(cbind(as.character(info3$AlleleA),as.character(info3$AlleleB), as.character(info3$Strand), as.character(info3$REF)),
                 1, function(x) { if (x[3] == "-") { x[1:2] <- base.tr[x[1:2]]; }
                                                     y <- setdiff(x[1:2], x[4]);
                                                     y <- ifelse(length(y) == 1, y, NA);
                                                     return(y) })

info3$alt <- snp.alt                                             

info3$QUAL="."
info3$INFO="Affymetrix_PigSNP130K"

info4= info3[!is.na(info3$alt),]

#info5 = info4[order(info4$chr,info4$pos,info4$Name),]
#去重重复位置的SNP
#info5$chrpos=paste(info5$chr,info5$pos,sep=":")
#info6 = info5[!duplicated(info5$chrpos),1:11]

#info7=subset(info6,select=c(chr,pos,Name,REF,alt,QUAL,QUAL,INFO))
info7=subset(info4,select=c(chr,pos,Name,REF,alt,QUAL,QUAL,INFO))
dim(info7)
write.table(info7,"Affymetrix_PigSNP130K_111522.info", row.names = FALSE, col.names =FALSE, quote =FALSE)



#info1=read.table("130K_info.txt",head=T)
#info2_1=read.table("pig_130k_A111753_info2.txt",head=T)
#info2_2=read.table("pig_130k_B5401_info2.txt",head=T)
#info2 = rbind(info2_1,info2_2)
#info3=merge(info2,info1,by = "Name", all.x = TRUE)

#info3$REF <- toupper(info3$REF)

#base.tr <- c("A", "T", "C", "G")
#names(base.tr) <- c("T", "A", "G", "C")
#snp.alt <- apply(cbind(as.character(info3$AlleleA),as.character(info3$AlleleB), as.character(info3$Strand), as.character(info3$REF)),
                 1, function(x) { if (x[3] == "-") { x[1:2] <- base.tr[x[1:2]]; }
                                                     y <- setdiff(x[1:2], x[4]);
                                                     y <- ifelse(length(y) == 1, y, NA);
                                                     return(y) })

#info3$alt <- snp.alt                                             

#info3$QUAL="."
#info3$INFO="Affymetrix_PigSNP130K"

#info4= info3[!is.na(info3$alt),]

#info5 = info4[order(info4$chr,info4$pos,info4$Name),]
#去重重复位置的SNP
#info5$chrpos=paste(info5$chr,info5$pos,sep=":")
#info6 = info5[!duplicated(info5$chrpos),1:11]

#info7=subset(info6,select=c(chr,pos,Name,REF,alt,QUAL,QUAL,INFO))
#dim(info6)
#write.table(info7,"Affymetrix_PigSNP130K_99018.info", row.names = FALSE, col.names =FALSE, quote =FALSE)
