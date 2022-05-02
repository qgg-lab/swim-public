# ===========================================
# = figure for PCA, genetic structure, etc. =
# ===========================================

args <- commandArgs(TRUE) # args = c("../reportData/all.ids.random.snp.count", "../reportData/Fig_muti.Landrace", "../reportData/Fig_muti.Yorkshire", "../reportData/Fig_muti.Duroc", "../reportData/Fig_muti.Meishan", "../reportData/Fig_muti.Tibetan", "../reportData/Fig_muti.WB_Asia", "../reportData/Fig_muti.WB_European", "../reportData/figurePCA.RData", "../reportData/SWIM_BreedInfo.txt", "../reportData/select.all.pruned.fam", "../reportData/select.all.pruned.2.Q", "../reportData/select.all.pruned.4.Q", "../reportData/select.all.pruned.6.Q", "../reportData/breed.code.csv", "../report/figure_genetics.pdf", "Myriad Pro")
library("RColorBrewer")

# prepare file to plot
# ============================================================

file.width = 180
cairo_pdf(file = args[16], width = file.width/25.4, height = file.width*0.45/25.4, family = args[17])
layout(matrix(c(1, 2, 2, 4, 3, 3, 5, 5, 3, 3, 6, 6, 3, 3, 7, 7), byrow = T, ncol = 4, nrow = 4), widths = c(0.18, 0.22, 0.05, 0.55), heights = c(1, 0.5, 0.5, 0.5))

# variants plot
# ============================================================

snp.count <- read.table(args[1], header = FALSE, as.is = TRUE)
par(las = 1, tcl = -0.2, mai = c(0.075, 0.1, 0.05, 0.01)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = T)
plot(c(0, 2500), c(0, max(snp.count[, 2])/1000000), type = "n", xlab = "", ylab = "", axes = FALSE)
points(c(0, snp.count[, 1]), c(0, snp.count[, 2]/1000000), type = "l", lwd = 1, col = "grey30")
axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), at = c(0, 1000, 2259), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(0, 40, 10), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "Number of animals", mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Number of variants (M)", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))

points(c(1000, 2259), c(snp.count[snp.count[, 1] == 1000, 2]/1000000, snp.count[snp.count[, 1] == 2259, 2]/1000000), pch = 4, lwd = 1, col = "grey30", cex = 0.8)
text(1000, snp.count[snp.count[, 1] == 1000, 2]/1000000 - 6, paste(formatC(snp.count[snp.count[, 1] == 1000, 2]/1000000, digits = 2, format = "f"), " M", sep = ""), xpd = TRUE, cex = 6/par("ps")/par("cex"))
text(2259, snp.count[snp.count[, 1] == 2259, 2]/1000000 - 4, paste(formatC(snp.count[snp.count[, 1] == 2259, 2]/1000000, digits = 2, format = "f"), " M", sep = ""), xpd = TRUE, cex = 6/par("ps")/par("cex"))

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# LD plot
# ============================================================

landrace.ld <- read.table(args[2], header = FALSE, as.is = TRUE)
yorkshire.ld <- read.table(args[3], header = FALSE, as.is = TRUE)
duroc.ld <- read.table(args[4], header = FALSE, as.is = TRUE)
meishan.ld <- read.table(args[5], header = FALSE, as.is = TRUE)
tibet.ld <- read.table(args[6], header = FALSE, as.is = TRUE)
asianWB.ld <- read.table(args[7], header = FALSE, as.is = TRUE)
euWB.ld <- read.table(args[8], header = FALSE, as.is = TRUE)


point.col <- c(brewer.pal(9, "Blues")[9], brewer.pal(8, "Dark2"), brewer.pal(9, "Reds")[9])
names(point.col) <- c("Asian_other_pigs", "Asian_Wild_Boars", "China_breed", "Duroc", "European_breed", "European_Wild_Boars", "Hybrid_pigs", "Landrace", "Large_White", "tibet")


par(las = 1, tcl = -0.2, mai = c(0.075, 0.1, 0.05, 0.05)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = T)
plot(c(0, 160000), c(0, 1), type = "n", xlab = "", ylab = "", axes = FALSE)
points(landrace.ld[landrace.ld[, 1] <= 100000, 1:2], type = "l", lwd = 1, col = point.col["Landrace"])
points(duroc.ld[duroc.ld[, 1] <= 100000, 1:2], type = "l", lwd = 1, col = point.col["Duroc"])
points(yorkshire.ld[yorkshire.ld[, 1] <= 100000, 1:2], type = "l", lwd = 1, col = point.col["Large_White"])
points(euWB.ld[euWB.ld[, 1] <= 100000, 1:2], type = "l", lwd = 1, col = point.col["European_Wild_Boars"])
# points(tibet.ld[tibet.ld[, 1] <= 100000, 1:2], type = "l", lwd = 1, col = point.col["tibet"])
points(asianWB.ld[asianWB.ld[, 1] <= 100000, 1:2], type = "l", lwd = 1, col = point.col["Asian_Wild_Boars"])
points(meishan.ld[meishan.ld[, 1] <= 100000, 1:2], type = "l", lwd = 1, col = point.col["China_breed"])

axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), at = c(0, 50000, 100000), label = c("0", "50,000", "100,000"), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "Distance (bp)", mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste("Average ", italic(r), {}^2)), mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))



text(110000, 0.54, "Duroc", col = point.col["Duroc"], pos = 4, cex = 6/par("ps")/par("cex"))
segments(102000, duroc.ld[duroc.ld[, 1] == 100000, 2], 106000, duroc.ld[duroc.ld[, 1] == 100000, 2], col = point.col["Duroc"])
segments(106000, duroc.ld[duroc.ld[, 1] == 100000, 2], 110000, 0.54, col = point.col["Duroc"])
segments(110000, 0.54, 114000, 0.54, col = point.col["Duroc"])

text(110000, 0.46, "Landrace", col = point.col["Landrace"], pos = 4, cex = 6/par("ps")/par("cex"))
segments(102000, landrace.ld[landrace.ld[, 1] == 100000, 2], 106000, landrace.ld[landrace.ld[, 1] == 100000, 2], col = point.col["Landrace"])
segments(106000, landrace.ld[landrace.ld[, 1] == 100000, 2], 110000, 0.46, col = point.col["Landrace"])
segments(110000, 0.46, 114000, 0.46, col = point.col["Landrace"])


text(110000, 0.38, "Yorkshire", col = point.col["Large_White"], pos = 4, cex = 6/par("ps")/par("cex"))
segments(102000, yorkshire.ld[yorkshire.ld[, 1] == 100000, 2], 106000, yorkshire.ld[yorkshire.ld[, 1] == 100000, 2], col = point.col["Large_White"])
segments(106000, yorkshire.ld[yorkshire.ld[, 1] == 100000, 2], 110000, 0.38, col = point.col["Large_White"])
segments(110000, 0.38, 114000, 0.38, col = point.col["Large_White"])


text(110000, 0.3, "Meishan", col = point.col["China_breed"], pos = 4, cex = 6/par("ps")/par("cex"))
segments(102000, meishan.ld[meishan.ld[, 1] == 100000, 2], 106000, meishan.ld[meishan.ld[, 1] == 100000, 2], col = point.col["China_breed"])
segments(106000, meishan.ld[meishan.ld[, 1] == 100000, 2], 110000, 0.30, col = point.col["China_breed"])
segments(110000, 0.30, 114000, 0.30, col = point.col["China_breed"])


text(110000, 0.16, "European wild", col = point.col["European_Wild_Boars"], pos = 4, cex = 6/par("ps")/par("cex"), xpd = T)
segments(102000, euWB.ld[euWB.ld[, 1] == 100000, 2], 106000, euWB.ld[euWB.ld[, 1] == 100000, 2], col = point.col["European_Wild_Boars"])
segments(106000, euWB.ld[euWB.ld[, 1] == 100000, 2], 110000, 0.16, col = point.col["European_Wild_Boars"])
segments(110000, 0.16, 114000, 0.16, col = point.col["European_Wild_Boars"])



text(110000, 0.08, "Asian wild", col = point.col["Asian_Wild_Boars"], pos = 4, cex = 6/par("ps")/par("cex"), xpd = T)
segments(102000, asianWB.ld[asianWB.ld[, 1] == 100000, 2], 106000, asianWB.ld[asianWB.ld[, 1] == 100000, 2], col = point.col["Asian_Wild_Boars"])
segments(106000, asianWB.ld[asianWB.ld[, 1] == 100000, 2], 110000, 0.08, col = point.col["Asian_Wild_Boars"])
segments(110000, 0.08, 114000, 0.08, col = point.col["Asian_Wild_Boars"])
#
#
# text(110000, 0, "Tibetan", col = point.col["tibet"], pos = 4, cex = 6/par("ps")/par("cex"), xpd = T)
# segments(102000, tibet.ld[tibet.ld[, 1] == 100000, 2], 106000, tibet.ld[tibet.ld[, 1] == 100000, 2], col = point.col["tibet"])
# segments(106000, tibet.ld[tibet.ld[, 1] == 100000, 2], 110000, 0, col = point.col["tibet"])
# segments(110000, 0, 114000, 0, col = point.col["tibet"])


text(grconvertX(0.05 + file.width/25.4*0.18, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)
text(grconvertX(file.width/25.4*0.45 - 0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

#  load data for PCA
# ============================================================

load(args[9])
breed.info <- read.table(args[10], header = T, as.is = T)
breed.info[breed.info[, 3] == "TT", 2] <- "tibet"
pca_combine <- merge(pca_raw, breed.info, "ID", all.x = T)
pca <- subset(pca_combine,select=c(pca1,pca2,Breed))

point.col <- c(brewer.pal(9, "Blues")[9], brewer.pal(8, "Dark2"), brewer.pal(9, "Reds")[9])
breed.count <- table(pca$Breed)
pca <- pca[order(factor(pca$Breed, levels = names(sort(breed.count, decreasing = T)))), ]
names(point.col) <- c("Asian_other_pigs", "Asian_Wild_Boars", "China_breed", "Duroc", "European_breed", "European_Wild_Boars", "Hybrid_pigs", "Landrace", "Large_White", "tibet")
point.col["tibet"] <- point.col["China_breed"]

#sapply(split(pca$pca1, pca$Breed), mean)
#sapply(split(pca$pca2, pca$Breed), mean)

par(las = 1, tcl = -0.2, mai = c(0.075, 0.1, 0.02, 0.01)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = T)

plot(c(-0.02, 0.06), c(-0.03, 0.05), type = "n", axes = FALSE, xlab = "", ylab = "")
points(pca$pca1, pca$pca2, col = point.col[pca$Breed], cex = 0.5)
axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(-0.02, 0.04, 0.02), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "PC1", mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "PC2", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))

text(mean(pca$pca1[pca$Breed == "Asian_other_pigs"]), mean(pca$pca2[pca$Breed == "Asian_other_pigs"]) + 0.004, "Korean breeds", col = point.col[1], cex = 7/par("ps")/par("cex"))
text(mean(pca$pca1[pca$Breed == "Asian_Wild_Boars"]) - 0.002, mean(pca$pca2[pca$Breed == "Asian_Wild_Boars"]) + 0.005, "Asian wild boars", col = point.col[2], cex = 7/par("ps")/par("cex"))
text(mean(pca$pca1[pca$Breed == "China_breed"]) - 0.004, mean(pca$pca2[pca$Breed == "China_breed"]) - 0.008, "Chinese breeds", col = point.col[3], cex = 7/par("ps")/par("cex"))
text(mean(pca$pca1[pca$Breed == "Duroc"]) + 0.008, mean(pca$pca2[pca$Breed == "Duroc"]) + 0.005, "Duroc", col = point.col[4], cex = 7/par("ps")/par("cex"))
text(mean(pca$pca1[pca$Breed == "European_breed"]) + 0.020, mean(pca$pca2[pca$Breed == "European_breed"]) - 0.004, "Other European breeds", col = point.col[5], cex = 7/par("ps")/par("cex"))
text(mean(pca$pca1[pca$Breed == "European_Wild_Boars"]) + 0.012, mean(pca$pca2[pca$Breed == "European_Wild_Boars"]) - 0.003, "European wild boars", col = point.col[6], cex = 7/par("ps")/par("cex"))
text(mean(pca$pca1[pca$Breed == "Hybrid_pigs"]), mean(pca$pca2[pca$Breed == "Hybrid_pigs"]) + 0.015, "Composite and hybrid", col = point.col[7], cex = 7/par("ps")/par("cex"))
text(mean(pca$pca1[pca$Breed == "Landrace"]), mean(pca$pca2[pca$Breed == "Landrace"]) - 0.007, "Landrace", col = point.col[8], cex = 7/par("ps")/par("cex"))
text(mean(pca$pca1[pca$Breed == "Large_White"]) + 0.010, mean(pca$pca2[pca$Breed == "Large_White"]), "Yorkshire", col = point.col[9], cex = 7/par("ps")/par("cex"))
# text(mean(pca$pca1[pca$Breed == "tibet"]), mean(pca$pca2[pca$Breed == "tibet"]) - 0.004, "Tibetan wild boars", col = point.col[10], cex = 7/par("ps")/par("cex"))

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# admixture plots
# ============================================================

breed.code <- read.csv(args[15], header = FALSE, as.is = TRUE)
rownames(breed.code) <- breed.code[, 2]
admix.id <- read.table(args[11], header = FALSE, as.is = TRUE)
admix.breed <- breed.info[match(admix.id[, 1], breed.info[, 1]), 3]
breed.order <- c("AWB","TT","RC","WJ",
				"NJ","XP","DNP","WZS","BMX","LUC","DHB",
				"XEH","LTZ","YDH", "LPS","TC","DWZ", "MS","EHL","WNB","JQH","WNS","JH","AQ",
				"MIN","HT","LWH","KP","KJBP","LCH","TWH",
				"WD","USMARC","DLY","DP","YL","DRC","LW",
				"LR","PIT","HP","BK","IBR","EWB")

breed.code["LW", 1] = "Yorkshire"
breed.code["TT", 1] = "Tibetan"
breed.code["AWB", 1] = "Asian wild boars"
breed.code["EWB", 1] = "European wild boars"


admix.breed.order <- admix.breed[order(as.numeric(factor(admix.breed, levels = breed.order)))]
admix.k2 <- read.table(args[12], header = FALSE, as.is = TRUE)[order(as.numeric(factor(admix.breed, levels = breed.order))), ]
admix.k4 <- read.table(args[13], header = FALSE, as.is = TRUE)[order(as.numeric(factor(admix.breed, levels = breed.order))), ]
admix.k6 <- read.table(args[14], header = FALSE, as.is = TRUE)[order(as.numeric(factor(admix.breed, levels = breed.order))), ]

par(las = 1, tcl = -0.2, mai = c(0, 0, 0.02, 0.02)*file.width/25.4/4, ps = 7, lwd = 0.5, xpd = T)

plot(c(0, 185), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(6, 179))
breed.mark <- numeric(length(breed.order) + 1)
for (i in 1:length(breed.order)) {
	breed.mark[i + 1] <- max(which(admix.breed.order == breed.order[i]))
	segments((breed.mark[i + 1] + breed.mark[i])/2, par('usr')[3], (breed.mark[i + 1] + breed.mark[i])/2, par('usr')[3]+ 0.02)
	segments((breed.mark[i + 1] + breed.mark[i])/2, par('usr')[3] + 0.02, i*4.1, par('usr')[3]+ 0.05)

}

text((1:length(breed.order))*4.1 - 3, rep(0.03, length(breed.order)), label = breed.code[breed.order, 1], srt = 90, pos = 4, cex = 6/par("ps")/par("cex"))

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), at = breed.mark, label = rep("", length(breed.mark)), cex.axis = 7/par("ps")/par("cex"))



# set color for ancestral populations
anc.col <- brewer.pal(8, "Dark2")[c(1:4, 6, 8)]


# SW china
segments(which(breed.order == "TT")*4.1, 0.45, which(breed.order == "NJ")*4.1, 0.45, col = anc.col[2], lwd = 1)
text((which(breed.order == "TT")*4.1 + which(breed.order == "NJ")*4.1)/2 - 4, 0.5, "SW China\n(SWCN)", col = anc.col[2], srt = 45, pos = 4, cex = 6/par("ps")/par("cex"))

# s China
segments(which(breed.order == "XP")*4.1, 0.5, which(breed.order == "YDH")*4.1, 0.5, col = anc.col[4], lwd = 1)
text((which(breed.order == "XP")*4.1 + which(breed.order == "YDH")*4.1)/2 - 4, 0.55, "South China\n(SCN)", col = anc.col[4], srt = 45, pos = 4, cex = 6/par("ps")/par("cex"))


# c China
segments(which(breed.order == "LPS")*4.1, 0.52, which(breed.order == "DWZ")*4.1, 0.52, col = anc.col[5], lwd = 1)
text((which(breed.order == "LPS")*4.1 + which(breed.order == "DWZ")*4.1)/2 - 4, 0.57, "Central China\n(CCN)", col = anc.col[5], srt = 45, pos = 4, cex = 6/par("ps")/par("cex"))

# e China
segments(which(breed.order == "MS")*4.1, 0.71, which(breed.order == "AQ")*4.1, 0.71, col = anc.col[5], lwd = 1)
text((which(breed.order == "MS")*4.1 + which(breed.order == "AQ")*4.1)/2 - 4, 0.76, "East China\n(ECN)", col = anc.col[5], srt = 45, pos = 4, cex = 6/par("ps")/par("cex"))

# n China
segments(which(breed.order == "MIN")*4.1, 0.28, which(breed.order == "LWH")*4.1, 0.28, col = anc.col[5], lwd = 1)
text((which(breed.order == "MIN")*4.1 + which(breed.order == "LWH")*4.1)/2 - 2, 0.33, "North China\n(NCN)", col = anc.col[5], srt = 90, pos = 4, cex = 6/par("ps")/par("cex"))

# korea
segments(which(breed.order == "KP")*4.1, 0.6, which(breed.order == "KJBP")*4.1, 0.6, col = anc.col[6], lwd = 1)
text((which(breed.order == "KP")*4.1 + which(breed.order == "KJBP")*4.1)/2 - 4, 0.65, "Korea", col = anc.col[6], srt = 45, pos = 4, cex = 6/par("ps")/par("cex"))

# composite and hybrid
segments(which(breed.order == "LCH")*4.1, 0.84, which(breed.order == "YL")*4.1, 0.84, col = anc.col[3], lwd = 1)
text((which(breed.order == "LCH")*4.1 + which(breed.order == "YL")*4.1)/2, 0.81, "Composite and\nhybrid (CH)", col = anc.col[3], pos = 3, cex = 6/par("ps")/par("cex"))

# composite and hybrid
segments(which(breed.order == "DRC")*4.1, 0.37, which(breed.order == "IBR")*4.1, 0.37, col = anc.col[1], lwd = 1)
text((which(breed.order == "DRC")*4.1 + which(breed.order == "IBR")*4.1)/2, 0.34, "Europe", col = anc.col[1], pos = 3, cex = 6/par("ps")/par("cex"))

par(las = 1, tcl = -0.2, mai = c(0.02, 0.2, 0, 0.02)*file.width/25.4/4, ps = 7, lwd = 0.5, xpd = F)


plot(c(0, 185), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(6, 179))
barplot(as.matrix(t(admix.k6)), space = 0, border = NA, col = anc.col[c(3, 4, 5, 2, 6, 1)], add = TRUE, axes = F, names.arg = rep("", length(admix.breed.order)))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(0, 1, 0.2), label = c(seq(0, 100, 20)), cex.axis = 7/par("ps")/par("cex"))
text(-2, 0.85, "K = 6", col = "white", cex = 7/par("ps")/par("cex"), pos = 4)

segments(min(which(admix.breed.order == "RC")) - 0.5, -0.02, max(which(admix.breed.order == "NJ")) - 0.5, -0.02, lwd = 1, xpd = TRUE, col = anc.col[2])
text((min(which(admix.breed.order == "RC")) + max(which(admix.breed.order == "NJ")))/2 - 0.5, -0.07, "SWCN", xpd = TRUE, cex = 4/par("ps")/par("cex"), col = anc.col[2])

segments(min(which(admix.breed.order == "XP")) - 0.5, -0.02, max(which(admix.breed.order == "YDH")) - 0.5, -0.02, lwd = 1, xpd = TRUE, col = anc.col[4])
text((min(which(admix.breed.order == "XP")) + max(which(admix.breed.order == "YDH")))/2 - 0.5, -0.07, "SCN", xpd = TRUE, cex = 4/par("ps")/par("cex"), col = anc.col[4])

segments(min(which(admix.breed.order == "LPS")) - 0.5, -0.02, max(which(admix.breed.order == "DWZ")) - 0.5, -0.02, lwd = 1, xpd = TRUE, col = anc.col[5])
text((min(which(admix.breed.order == "LPS")) + max(which(admix.breed.order == "DWZ")))/2 - 0.5, -0.07, "CCN", xpd = TRUE, cex = 4/par("ps")/par("cex"), col = anc.col[5])

segments(min(which(admix.breed.order == "MS")) - 0.5, -0.02, max(which(admix.breed.order == "AQ")) - 0.5, -0.02, lwd = 1, xpd = TRUE, col = anc.col[5])
text((min(which(admix.breed.order == "MS")) + max(which(admix.breed.order == "AQ")))/2 - 0.5, -0.07, "ECN", xpd = TRUE, cex = 4/par("ps")/par("cex"), col = anc.col[5])

segments(min(which(admix.breed.order == "MIN")) - 0.5, -0.02, max(which(admix.breed.order == "LWH")) - 0.5, -0.02, lwd = 1, xpd = TRUE, col = anc.col[5])
text((min(which(admix.breed.order == "MIN")) + max(which(admix.breed.order == "LWH")))/2 - 0.5, -0.07, "NCN", xpd = TRUE, cex = 4/par("ps")/par("cex"), col = anc.col[5])

segments(min(which(admix.breed.order == "KP")) - 0.5, -0.02, max(which(admix.breed.order == "KJBP")) - 0.5, -0.02, lwd = 1, xpd = TRUE, col = anc.col[6])
text((min(which(admix.breed.order == "KP")) + max(which(admix.breed.order == "KJBP")))/2 - 0.5, -0.07, "Korea", xpd = TRUE, cex = 4/par("ps")/par("cex"), col = anc.col[6])

segments(min(which(admix.breed.order == "LCH")) - 0.5, -0.02, max(which(admix.breed.order == "YL")) - 0.5, -0.02, lwd = 1, xpd = TRUE, col = anc.col[3])
text((min(which(admix.breed.order == "LCH")) + max(which(admix.breed.order == "YL")))/2 - 0.5, -0.07, "CH", xpd = TRUE, cex = 4/par("ps")/par("cex"), col = anc.col[3])

segments(min(which(admix.breed.order == "DRC")) - 0.5, -0.02, max(which(admix.breed.order == "IBR")) - 0.5, -0.02, lwd = 1, xpd = TRUE, col = anc.col[1])
text((min(which(admix.breed.order == "DRC")) + max(which(admix.breed.order == "IBR")))/2 - 0.5, -0.07, "Europe", xpd = TRUE, cex = 4/par("ps")/par("cex"), col = anc.col[1])

text((min(which(admix.breed.order == "DRC")) + max(which(admix.breed.order == "DRC")))/2 - 0.5, 0.5, "Duroc", cex = 5/par("ps")/par("cex"), col = "white")

plot(c(0, 185), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(6, 179))
barplot(as.matrix(t(admix.k4)), space = 0, border = NA, col = anc.col[c(4, 1, 2, 3)], add = TRUE, axes = F, names.arg = rep("", length(admix.breed.order)))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(0, 1, 0.2), label = c(seq(0, 100, 20)), cex.axis = 7/par("ps")/par("cex"))
title(ylab = "Ancestry (%)", mgp = c(1.5, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(-2, 0.85, "K = 4", col = "white", cex = 7/par("ps")/par("cex"), pos = 4)

plot(c(0, 185), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(6, 179))
barplot(as.matrix(t(admix.k2)), space = 0, border = NA, col = anc.col[1:2], add = TRUE, axes = F, names.arg = rep("", length(admix.breed.order)))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(0, 1, 0.2), label = c(seq(0, 100, 20)), cex.axis = 7/par("ps")/par("cex"))
text(-2, 0.85, "K = 2", col = "white", cex = 7/par("ps")/par("cex"), pos = 4)




dev.off()


