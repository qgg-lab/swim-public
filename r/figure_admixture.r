##R代码做structure图

library(RColorBrewer)

rm(list=ls());options(stringsAsFactors=FALSE)

breedsOrder = c("AWB","MS","EHL","WNB","JQH","WNS","JH","AQ","MIN","HT","LWH","LPS","TC","DWZ","RC","WJ","NJ","TT","XP","DNP","WZS","BMX","LUC","DHB","XEH","LTZ","YDH","KP","KJBP","LCH","TWH","WD","USMARC","DLY","DP","YL","DRC","LW","LR","PIT","HP","BK","IBR","EWB")

colors <- c( brewer.pal(9, "Reds")[9], brewer.pal(8, "Dark2")[1], brewer.pal(8, "Greens")[8], brewer.pal(9, "Purples")[8], brewer.pal(8, "Paired")[8], brewer.pal(8, "Paired")[6],brewer.pal(9, "Blues")[8], brewer.pal(9, "BuPu")[8],brewer.pal(9, "Set1")[6],brewer.pal(9, "PiYG")[2])

#colors = c("blue2","green","mediumturquoise","orange","chartreuse4","chocolate4","coral1","chartreuse4","darkred","lightgoldenrod3","burlywood1","darkolivegreen2","aquamarine1","darkslategray1","lightslateblue","blueviolet","orange",
#"red","blue","grey","olivedrab1","hotpink4","lightcoral","chocolate4","springgreen","lavenderblush","lightgoldenrod","lightpink4","lightskyblue3","lightyellow3","limegreen","mediumvioletred","navyblue","orchid1","paleturquoise3","peachpuff", "rosybrown1","seashell2","tan3","tomato3","turquoise2","violetred"
#,"yellow1")

for(n in 2:12){
	fam = read.table(file="select.all.pruned.fam",header=FALSE,colClasses="character")
	dat = read.table(paste("select.all.pruned.",n,".Q",sep=""),header=FALSE) 
	fam2=cbind(fam,dat)
	colnames(fam2)=c("ID","ID1","un1","un2","un3","un4",paste("un",5:(4+n),sep=""))
	info = read.table(file="SWIM_sampleInfo_V3.txt",header=TRUE,colClasses="character")
	fam3 = merge(fam2,info,"ID",all.x=T)
	famID =fam3[,7:(7+n)]
	newDat = NULL
	mkr = 0
	tmp = NULL
for(i in 1:length(breedsOrder)){
	tmp = famID [famID$BreedAbbr== breedsOrder[i],-((1+n):(n+2))]
	newDat = rbind(newDat,tmp)
	N_ind = nrow(tmp)
	mkr = c(mkr,mkr[length(mkr)]+N_ind)#下标
}
	dim(newDat)
	mkr = mkr[-1]
	at = numeric(length(mkr))
	at[1] = mkr[1]/2
for(m in 2:length(at)){
	at[m] = (mkr[m]+mkr[m-1])/2
}	
png(paste("barplot.",n,".Q",".png",sep=""),width=8000,height=1200,res=350)
	par(mar = c(5.5, 2, 1.1, 0.9))
		barplot(t(as.matrix(newDat)), col=colors[1:ncol(newDat)], ylim=c(-0.01,1.14),
		xlab = "", ylab = "", main = "", axes = FALSE, border=NA, space=0, cex.names=NULL, xaxt="n", xaxt = "n")
		axis(side = 2, cex.axis =1.2, font=1, las = 1, gap.axis=0.25, line=(-3))
	
for(k in 1:length(mkr)){
		segments(mkr[k], 0, mkr[k], 1, lty=2,lwd = 1)
		axis(side=1,at=at,labels=breedsOrder,cex.axis=1.2,font=1,las=2)
}

####if necessary
h1=1.04
h2=1.14
lwd1=5.5

segments(0, h1, 8, h1, lwd = lwd1, col = brewer.pal(8, "Purples")[8])
text(4.5, h2, "AWB", xpd = TRUE, cex = 1.2, col = brewer.pal(8, "Purples")[8])

segments(8, h1, 26, h1, lwd = lwd1, col = brewer.pal(9, "Reds")[9])
text(17, h2, "ECN", xpd = TRUE, cex = 1.2, col = brewer.pal(9, "Reds")[9])

segments(26, h1, 42, h1, lwd = lwd1, col = brewer.pal(8, "Dark2")[1])
text(34, h2, "NCN", xpd = TRUE, cex = 1.2, col = brewer.pal(8, "Dark2")[1])

segments(42, h1, 49, h1, lwd = lwd1, col = brewer.pal(8, "Greens")[8])
text(45.5, h2, "CCN", xpd = TRUE, cex = 1.2, col = brewer.pal(8, "Greens")[8])

segments(49, h1, 62, h1, lwd = lwd1, col = brewer.pal(8, "Set1")[8])
text(55.5, h2, "SWCN", xpd = TRUE, cex = 1.2, col = brewer.pal(8, "Set1")[8])

segments(62, h1, 94, h1, lwd = lwd1, col = brewer.pal(8, "Paired")[6])
text(78, h2, "SCN", xpd = TRUE, cex = 1.2, col = brewer.pal(8, "Paired")[6])

segments(94, h1, 102, h1, lwd = lwd1, col = brewer.pal(9, "Paired")[8])
text(98, h2, "AND", xpd = TRUE, cex = 1.2, col = brewer.pal(9, "Paired")[8])

segments(102, h1, 125, h1, lwd = lwd1, col = brewer.pal(9, "Blues")[8])
text(113.5, h2, "Hybrid", xpd = TRUE, cex = 1.2, col = brewer.pal(9, "Blues")[8])

segments(125, h1, 177, h1, lwd = lwd1, col = brewer.pal(9, "BuPu")[8])
text(151, h2, "EUD", xpd = TRUE, cex = 1.2, col = brewer.pal(9, "BuPu")[8])

segments(177, h1, 185, h1, lwd = lwd1, col = brewer.pal(9, "Reds")[9])
text(181, h2 , "EWB", xpd = TRUE, cex = 1.2, col = brewer.pal(9, "Reds")[9])

text(-2.5, h2 , "K=6", xpd = TRUE, cex = 1.3, col = "black")
	
	dev.off()
}
