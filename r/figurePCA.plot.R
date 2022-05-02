library("RColorBrewer")
#brewer.pal.all()
par(las = 1, tcl = -0.2, mar = c(2, 2.5, 1, 0.9), ps = 7, lwd = 0.5)
date.col <- c( brewer.pal(9, "Paired")[8],brewer.pal(8, "Purples")[8],
			   brewer.pal(8, "Paired")[6], brewer.pal(9, "Reds")[9], 
			   brewer.pal(8, "Set1")[8],  brewer.pal(8, "Dark2")[1],
			   brewer.pal(8, "Greens")[8], brewer.pal(9, "Blues")[8], brewer.pal(9, "BuPu")[8])

info1=read.table("pca.p20.g80.het.filter.pruned.eigenvec",head=F)
info2=read.csv("SWIM_sampleInfo_V2.csv",head=T)
colnames(info1)=c("ID","ID1","pca1","pca2","pca3","pca4","pca5")
pca=merge(info1,info2,"ID",all.x=T)
pca_select <- subset(pca,select=c(pca1,pca2,Breed6))

png("pca.png",width=1800,height=1800,res=350)
par(tcl = -0.2, mar = c(2.7, 2.8, 1, 0.9), ps = 7, lwd = 1.1)

plot(c(-0.016, 0.06), c(-0.02, 0.045), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)
points(pca_select[pca_select$Breed6=="Landrace","pca1"], pca_select[pca_select$Breed6=="Landrace","pca2"],
col = brewer.pal(9, "BuPu")[8] , pch = 1, cex = 1.1)
points(pca_select[pca_select$Breed6=="Large_White","pca1"], pca_select[pca_select$Breed6=="Large_White","pca2"], 
col = brewer.pal(9, "Blues")[8], pch = 1, cex = 1.1)
points(pca_select[pca_select$Breed6=="European_breed","pca1"], pca_select[pca_select$Breed6=="European_breed","pca2"], 
col = brewer.pal(8, "Set1")[8], pch = 1, cex = 1.1)
points(pca_select[pca_select$Breed6=="Hybrid_pigs","pca1"], pca_select[pca_select$Breed6=="Hybrid_pigs","pca2"], 
col = brewer.pal(8, "Dark2")[1], pch = 1, cex = 1.1)
points(pca_select[pca_select$Breed6=="European_Wild_Boars","pca1"], pca_select[pca_select$Breed6=="European_Wild_Boars","pca2"], 
col = brewer.pal(9, "Reds")[9], pch = 1, cex = 1.1)
points(pca_select[pca_select$Breed6=="Duroc","pca1"], pca_select[pca_select$Breed6=="Duroc","pca2"], 
col = brewer.pal(8, "Greens")[8], pch = 1, cex = 1.1)
points(pca_select[pca_select$Breed6=="China_breed","pca1"], pca_select[pca_select$Breed6=="China_breed","pca2"], 
col = brewer.pal(8, "Paired")[6], pch = 1, cex = 1.1)
points(pca_select[pca_select$Breed6=="Asian_Wild_Boars","pca1"], pca_select[pca_select$Breed6=="Asian_Wild_Boars","pca2"], 
col = brewer.pal(8, "Purples")[8], pch = 1, cex = 1.1)
points(pca_select[pca_select$Breed6=="Asian_other_pigs","pca1"], pca_select[pca_select$Breed6=="Asian_other_pigs","pca2"], 
col = brewer.pal(9, "Paired")[8], pch = 1, cex = 1.1)


legend("topright", 
bty = "n", 
pch = c(1,1,1,1,1,1,1) , 
col = date.col,
 legend = c("Asian other pigs","Asian Wild Boars","China breed","European Wild Boars","European breed","Hybrid pigs","Duroc","Large White","Landrace"),
 cex =7/par("ps")/par("cex"), 
 text.col = date.col)   
		
axis(side = 1, mgp = c(1.2, 0.3, 0), cex.axis = 10/par("ps")/par("cex"), lwd = 0.7, at = seq(-0.01, 0.05, 0.01),las = 1)
axis(side = 2, lwd = 0.7, cex.axis = 10/par("ps")/par("cex"), mgp = c(1.2, 0.4, 0), at = seq(-0.01, 0.04, 0.01),las = 1)
title(xlab = "PC1", cex.lab = 10/par("ps")/par("cex"), mgp = c(1.5, 0, 0))
title(ylab = "PC2", cex.lab = 10/par("ps")/par("cex"), mgp = c(1.8, 0, 0))

box(bty = "l")

dev.off()

save(pca,pca_select, file="figurePCA.RData")
