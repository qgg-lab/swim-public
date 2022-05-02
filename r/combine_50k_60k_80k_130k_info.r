#combine 50k 60k 80k 130k snpchip info
k50=read.table("PorcineSNP50K_45164.info",head=F)
k60=read.table("PorcineSNP60K_52866.info",head=F)
k80=read.table("PorcineSNP80K_60818.info",head=F)
k130=read.table("Affymetrix_PigSNP130K_99018.info",head=F)
comb =rbind(k50,k60,k80,k130)
comb1 = comb[order(comb$V1,comb$V2,comb$V3),]
#去重重复位置的SNP
comb1$chrpos=paste(comb1$V1,comb1$V2,sep=":")
comb2=comb1[!duplicated(comb1$chrpos),1:8]
dim(comb2)
write.table(comb2,"PorcineSNPcomb_138162.info", row.names = FALSE, col.names =FALSE, quote =FALSE)
