# ================================
# = figure for backfat thickness =
# ================================

args <- commandArgs(TRUE) # args = c("../report/figure_bf.pdf", "Myriad Pro", "../reportData/bf.s21.pval.RData", "../reportData/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.fai", "../reportData/bf.man.png", "../reportData/bf.s21.impute.pval.r2.RData", "../reportData/chr1.mc4r.gene.bound")
library("RColorBrewer")
library("png")

# prepare file to plot
# ============================================================

file.width = 90
cairo_pdf(file = args[1], width = file.width/25.4, height = file.width/25.4*1.4, family = args[2])
layout(matrix(c(1, 2, 3), byrow = T, ncol = 1, nrow = 3), heights = c(0.5, 0.60, 0.30))

# manhattan plot
# ============================================================

par(las = 1, tcl = -0.2, mai = c(0.12, 0.15, 0.1, 0.05)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = T)

load(args[3])

chr.len <- as.numeric(read.table(args[4], header = FALSE, as.is = TRUE)[1:18, 2])
names(chr.len) <- as.character(1:18)
chr.cum <- c(0, cumsum(chr.len)[1:17]) + (0:17)*1000000
names(chr.cum) <- as.character(1:18)
chr.col <- rep(c(brewer.pal(9, "Blues")[4], brewer.pal(9, "Blues")[8]), 9)

plot(c(0, chr.cum[18] + chr.len[18]), c(0, 14), type = "n", axes = FALSE, xaxs = "i", yaxs = "i", xlim = c(-2e7, chr.cum[18] + chr.len[18] + 2e7), ylim = c(-0.5, 14))
png.plot <- readPNG(args[5])
rasterImage(png.plot, -2e7, -0.5, chr.cum[18] + chr.len[18] + 2e7, 14)
points(chr.cum[as.character(chip.pval[, 1])] + chip.pval[, 3], -log10(chip.pval[, 9]), bg = chr.col[chip.pval[, 1]], pch = 21, cex = 0.1, lwd = 0.2, col = NULL)
axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), at = chr.cum + chr.len/2, label = 1:18, cex.axis = 5/par("ps")/par("cex"))
axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), at = chr.cum[c(11, 17)] + chr.len[c(11, 17)]/2, label = c(11, 17), cex.axis = 5/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "Chromosome", mgp = c(0.7, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste("-log", {}[10], "(", italic(P), "-value)")), mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))

points(c(1.8e9, 1.85e9), c(13.5, 13.5), pch = 21, cex = 0.5, bg = c(brewer.pal(9, "Blues")[4], brewer.pal(9, "Blues")[8]), col = NULL)
text(1.84e9, 13.5, "GGP Porcine 50K", pos = 4, cex = 6/par("ps")/par("cex"))

points(c(1.8e9, 1.85e9), c(12.5, 12.5), pch = 21, cex = 0.5, bg = c("grey50", "grey80"), col = NULL)
text(1.84e9, 12.5, "Imputed", pos = 4, cex = 6/par("ps")/par("cex"))

points(chr.cum[1] + chip.pval[chip.pval[, 2] == "12784636", 3], -log10(chip.pval[chip.pval[, 2] == "12784636", 9]), cex = 0.3)
points(chr.cum[1] + 161511936, -log10(2.98122e-13), cex = 0.3)

arrows(2.5e8, 11.5, chr.cum[1] + chip.pval[chip.pval[, 2] == "12784636", 3] + 10000000, -log10(chip.pval[chip.pval[, 2] == "12784636", 9]), length = 0.03)
text(2.2e8, 11.2, expression(paste("Array ", italic("MC4R"), " SNP")), pos = 4, cex = 6/par("ps")/par("cex"))
arrows(2.5e8, 13, chr.cum[1] + 161511936 + 10000000, -log10(2.98122e-13), length = 0.03)
text(2.2e8, 13, expression(paste("Imputed chr1:161511936:T>C")), pos = 4, cex = 6/par("ps")/par("cex"))

text(chr.cum[8], 13, expression(paste(italic(h), {}["array"]^2, " = 0.27 (0.02)")), cex = 6/par("ps")/par("cex"), pos = 4)
text(chr.cum[8], 11.5, expression(paste(italic(h), {}["imputed"]^2, " = 0.28 (0.02)")), cex = 6/par("ps")/par("cex"), pos = 4)

box(bty = "l")

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# local plot
# ============================================================

load(args[6])

par(las = 1, tcl = -0.2, mai = c(0.01, 0.15, 0.1, 0.05)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = T)

impute.pval <- impute.pval[nrow(impute.pval):1, ]

plot(c(158500000, 162500000), c(0, 14), type = "n", axes = FALSE, xaxs = "i", yaxs = "i", xlim = c(158450000, 162550000), ylim = c(-0.5, 14))

r2.col <- brewer.pal(9, "Blues")[c(2, 3, 5, 7, 9)]

rect(xleft = rep(158500000, length(r2.col)),
     ybottom = seq(length = length(r2.col), from = 10.5, to = 13.5)
               - (13.5-10.5)/(length(r2.col) - 1),
     xright = rep(158700000, length(r2.col)),
     ytop = seq(length = length(r2.col), from = 10.5, to = 13.5),
     col = r2.col)
     
text(158600000, 14, expression(paste(italic(r), {}^2)), cex = 6/par("ps")/par("cex"))
text(rep(158650000, 6),
     c(seq(length = length(r2.col), from = 10.5, to = 13.5) - (13.5-10.5)/(length(r2.col) - 1),
       seq(length = length(r2.col), from = 10.5, to = 13.5)[5]),
       formatC(seq(0, 1, 0.2), format = "f", digits = 1), pos = 4)

points(impute.pval[, 3], -log10(impute.pval[, 9]), pch = 21, bg = r2.col[as.numeric(cut(r2[match(impute.pval[, 2], r2[, 6]), 7], breaks = seq(0, 1, 0.2), include.lowest = TRUE))], col = NULL, cex = 0.8)
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
title(ylab = expression(paste("-log", {}[10], "(", italic(P), "-value)")), mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))

points(impute.pval[impute.pval[, 3] == 160773437, 3], -log10(impute.pval[impute.pval[, 3] == 160773437, 9]), cex = 1.2)
points(impute.pval[impute.pval[, 3] == 161511936, 3], -log10(impute.pval[impute.pval[, 3] == 161511936, 9]), cex = 1.2)

arrows(impute.pval[impute.pval[, 3] == 161511936, 3], 13.5, impute.pval[impute.pval[, 3] == 161511936, 3], -log10(impute.pval[impute.pval[, 3] == 161511936, 9]) + 0.2, length = 0.03)
text(impute.pval[impute.pval[, 3] == 161511936, 3], 13.2, expression(paste("Imputed chr1:161511936:T>C")), pos = 3, cex = 6/par("ps")/par("cex"))

arrows(impute.pval[impute.pval[, 3] == 160773437, 3], 12.5, impute.pval[impute.pval[, 3] == 160773437, 3], -log10(impute.pval[impute.pval[, 3] == 160773437, 9]) + 0.2, length = 0.03)
text(impute.pval[impute.pval[, 3] == 160773437, 3], 12.2, expression(paste("Imputed ", italic("MC4R"), " SNP")), pos = 3, cex = 6/par("ps")/par("cex"))

box(bty = "l")
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# mark chip SNPs

chip.impute.pval <- impute.pval[na.omit(match(paste(chip.pval[, 1], chip.pval[, 3], sep = "_"), paste(impute.pval[, 1], impute.pval[, 3], sep = "_"))), ]
points(chip.impute.pval[, 3], -log10(chip.impute.pval[, 9]), cex = 0.8, pch = 4)

legend(159000000, 14.8, pch = c(21, 4), pt.bg = c(brewer.pal(9, "Blues")[9], "black"), legend = c("Imputed", "Array"), bty = "n", x.intersp = 0.5, y.intersp = 0.8, cex = 6/par("ps")/par("cex"))

# genes
# ============================================================

genes <- read.table(args[7], header = FALSE, as.is = TRUE)

par(las = 1, tcl = -0.2, mai = c(0.15, 0.15, 0.01, 0.05)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = T)

plot(c(158500000, 162500000), c(0, 2), type = "n", axes = FALSE, xaxs = "i", yaxs = "i", xlim = c(158450000, 162550000), ylim = c(0, 2))
rect(xleft = genes[, 2],
     ybottom = rep(0.9, nrow(genes)),
     xright = genes[, 3],
     ytop = rep(1.1, nrow(genes)),
     col = ifelse(genes[, 4] == "+", brewer.pal(9, "Reds")[9], brewer.pal(9, "Blues")[9]), border = NA)
		 
		 


# segments((genes[genes[, 4] == "+" & genes[, 5] != ".", 2] + genes[genes[, 4] == "+" & genes[, 5] != ".", 3])/2, rep(1.1, nrow(genes[genes[, 4] == "+" & genes[, 5] != ".", ])), (genes[genes[, 4] == "+" & genes[, 5] != ".", 2] + genes[genes[, 4] == "+" & genes[, 5] != ".", 3])/2, rep(1.2, nrow(genes[genes[, 4] == "+" & genes[, 5] != ".", ])), lwd = 0.5, col = brewer.pal(9, "Reds")[9])

# genes that can go directly above
# horizontal

segments((genes[genes[, 5] == "MC4R", 2] + genes[genes[, 5] == "MC4R", 3])/2, 1.1, (genes[genes[, 5] == "MC4R", 2] + genes[genes[, 5] == "MC4R", 3])/2, 1.5, lwd = 0.5, col = brewer.pal(9, "Reds")[9])
text((genes[genes[, 5] == "MC4R", 2] + genes[genes[, 5] == "MC4R", 3])/2 - 80000, 1.5, expression(italic(">MC4R")), col = brewer.pal(9, "Reds")[9], pos = 4)

segments((genes[genes[, 5] == "PIGN", 2] + genes[genes[, 5] == "PIGN", 3])/2, 1.1, (genes[genes[, 5] == "PIGN", 2] + genes[genes[, 5] == "PIGN", 3])/2, 1.6, lwd = 0.5, col = brewer.pal(9, "Reds")[9])
text((genes[genes[, 5] == "PIGN", 2] + genes[genes[, 5] == "PIGN", 3])/2 - 80000, 1.6, expression(italic(">PIGN")), col = brewer.pal(9, "Reds")[9], pos = 4)

segments((genes[genes[, 5] == "RNF152", 2] + genes[genes[, 5] == "RNF152", 3])/2, 1.1, (genes[genes[, 5] == "RNF152", 2] + genes[genes[, 5] == "RNF152", 3])/2, 1.45, lwd = 0.5, col = brewer.pal(9, "Reds")[9])
text((genes[genes[, 5] == "RNF152", 2] + genes[genes[, 5] == "RNF152", 3])/2 - 80000, 1.45, expression(italic(">RNF152")), col = brewer.pal(9, "Reds")[9], pos = 4)

segments((genes[genes[, 5] == "CCBE1", 2] + genes[genes[, 5] == "CCBE1", 3])/2, 1.1, (genes[genes[, 5] == "CCBE1", 2] + genes[genes[, 5] == "CCBE1", 3])/2, 1.8, lwd = 0.5, col = brewer.pal(9, "Reds")[9])
text((genes[genes[, 5] == "CCBE1", 2] + genes[genes[, 5] == "CCBE1", 3])/2 - 80000, 1.8, expression(italic(">CCBE1")), col = brewer.pal(9, "Reds")[9], pos = 4)

segments((genes[genes[, 5] == "LMAN1", 2] + genes[genes[, 5] == "LMAN1", 3])/2, 1.1, (genes[genes[, 5] == "LMAN1", 2] + genes[genes[, 5] == "LMAN1", 3])/2, 1.65, lwd = 0.5, col = brewer.pal(9, "Reds")[9])
text((genes[genes[, 5] == "LMAN1", 2] + genes[genes[, 5] == "LMAN1", 3])/2 - 80000, 1.65, expression(italic(">LMAN1")), col = brewer.pal(9, "Reds")[9], pos = 4)

segments((genes[genes[, 5] == "CPLX4", 2] + genes[genes[, 5] == "CPLX4", 3])/2, 1.1, (genes[genes[, 5] == "CPLX4", 2] + genes[genes[, 5] == "CPLX4", 3])/2, 1.5, lwd = 0.5, col = brewer.pal(9, "Reds")[9])
text((genes[genes[, 5] == "CPLX4", 2] + genes[genes[, 5] == "CPLX4", 3])/2 - 80000, 1.5, expression(italic(">CPLX4")), col = brewer.pal(9, "Reds")[9], pos = 4)

segments((genes[genes[, 5] == "RAX", 2] + genes[genes[, 5] == "RAX", 3])/2, 1.1, (genes[genes[, 5] == "RAX", 2] + genes[genes[, 5] == "RAX", 3])/2, 1.35, lwd = 0.5, col = brewer.pal(9, "Reds")[9])
text((genes[genes[, 5] == "RAX", 2] + genes[genes[, 5] == "RAX", 3])/2 - 80000, 1.35, expression(italic(">RAX")), col = brewer.pal(9, "Reds")[9], pos = 4)

segments((genes[genes[, 5] == "ALPK2", 2] + genes[genes[, 5] == "ALPK2", 3])/2, 1.1, (genes[genes[, 5] == "ALPK2", 2] + genes[genes[, 5] == "ALPK2", 3])/2, 1.5, lwd = 0.5, col = brewer.pal(9, "Reds")[9])
text((genes[genes[, 5] == "ALPK2", 2] + genes[genes[, 5] == "ALPK2", 3])/2 - 80000, 1.5, expression(italic(">ALPK2")), col = brewer.pal(9, "Reds")[9], pos = 4)

segments((genes[genes[, 5] == "PHLPP1", 2] + genes[genes[, 5] == "PHLPP1", 3])/2, 0.9, (genes[genes[, 5] == "PHLPP1", 2] + genes[genes[, 5] == "PHLPP1", 3])/2, 0.2, lwd = 0.5, col = brewer.pal(9, "Blues")[9])
text((genes[genes[, 5] == "PHLPP1", 2] + genes[genes[, 5] == "PHLPP1", 3])/2 - 80000, 0.2, expression(italic("<PHLPP1")), col = brewer.pal(9, "Blues")[9], pos = 4)

segments((genes[genes[, 5] == "ZCCHC2", 2] + genes[genes[, 5] == "ZCCHC2", 3])/2, 0.9, (genes[genes[, 5] == "ZCCHC2", 2] + genes[genes[, 5] == "ZCCHC2", 3])/2, 0.35, lwd = 0.5, col = brewer.pal(9, "Blues")[9])
text((genes[genes[, 5] == "ZCCHC2", 2] + genes[genes[, 5] == "ZCCHC2", 3])/2 - 80000, 0.35, expression(italic("<ZCCHC2")), col = brewer.pal(9, "Blues")[9], pos = 4)

segments((genes[genes[, 5] == "TNFRSF11A", 2] + genes[genes[, 5] == "TNFRSF11A", 3])/2, 0.9, (genes[genes[, 5] == "TNFRSF11A", 2] + genes[genes[, 5] == "TNFRSF11A", 3])/2, 0.5, lwd = 0.5, col = brewer.pal(9, "Blues")[9])
text((genes[genes[, 5] == "TNFRSF11A", 2] + genes[genes[, 5] == "TNFRSF11A", 3])/2 - 80000, 0.5, expression(italic("<TNFRSF11A")), col = brewer.pal(9, "Blues")[9], pos = 4)

segments((genes[genes[, 5] == "RELCH", 2] + genes[genes[, 5] == "RELCH", 3])/2, 0.9, (genes[genes[, 5] == "RELCH", 2] + genes[genes[, 5] == "RELCH", 3])/2, 0.65, lwd = 0.5, col = brewer.pal(9, "Blues")[9])
text((genes[genes[, 5] == "RELCH", 2] + genes[genes[, 5] == "RELCH", 3])/2 - 80000, 0.65, expression(italic("<RELCH")), col = brewer.pal(9, "Blues")[9], pos = 4)


segments((genes[genes[, 5] == "CDH20", 2] + genes[genes[, 5] == "CDH20", 3])/2, 0.9, (genes[genes[, 5] == "CDH20", 2] + genes[genes[, 5] == "CDH20", 3])/2, 0.65, lwd = 0.5, col = brewer.pal(9, "Blues")[9])
text((genes[genes[, 5] == "CDH20", 2] + genes[genes[, 5] == "CDH20", 3])/2 - 80000, 0.65, expression(italic("<CDH20")), col = brewer.pal(9, "Blues")[9], pos = 4)


segments((genes[genes[, 5] == "GRP", 2] + genes[genes[, 5] == "GRP", 3])/2, 0.9, (genes[genes[, 5] == "GRP", 2] + genes[genes[, 5] == "GRP", 3])/2, 0.2, lwd = 0.5, col = brewer.pal(9, "Blues")[9])
text((genes[genes[, 5] == "GRP", 2] + genes[genes[, 5] == "GRP", 3])/2 - 80000, 0.2, expression(italic("<GRP")), col = brewer.pal(9, "Blues")[9], pos = 4)



segments((genes[genes[, 5] == "SEC11C", 2] + genes[genes[, 5] == "SEC11C", 3])/2, 0.9, (genes[genes[, 5] == "SEC11C", 2] + genes[genes[, 5] == "SEC11C", 3])/2, 0.35, lwd = 0.5, col = brewer.pal(9, "Blues")[9])
text((genes[genes[, 5] == "SEC11C", 2] + genes[genes[, 5] == "SEC11C", 3])/2 - 80000, 0.35, expression(italic("<SEC11C")), col = brewer.pal(9, "Blues")[9], pos = 4)


segments((genes[genes[, 5] == "ZNF532", 2] + genes[genes[, 5] == "ZNF532", 3])/2, 0.9, (genes[genes[, 5] == "ZNF532", 2] + genes[genes[, 5] == "ZNF532", 3])/2, 0.5, lwd = 0.5, col = brewer.pal(9, "Blues")[9])
text((genes[genes[, 5] == "ZNF532", 2] + genes[genes[, 5] == "ZNF532", 3])/2 - 80000, 0.5, expression(italic("<ZNF532")), col = brewer.pal(9, "Blues")[9], pos = 4)

segments((genes[genes[, 5] == "MALT1", 2] + genes[genes[, 5] == "MALT1", 3])/2, 0.9, (genes[genes[, 5] == "MALT1", 2] + genes[genes[, 5] == "MALT1", 3])/2, 0.65, lwd = 0.5, col = brewer.pal(9, "Blues")[9])
text((genes[genes[, 5] == "MALT1", 2] + genes[genes[, 5] == "MALT1", 3])/2 - 80000, 0.65, expression(italic("<MALT1")), col = brewer.pal(9, "Blues")[9], pos = 4)

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(159000000, 162000000, 1000000), label = seq(159, 162, 1), cex.axis = 6/par("ps")/par("cex"))
title(xlab = "Position on Chromosome 1 (Mb)", mgp = c(0.7, 0, 0), cex.lab = 7/par("ps")/par("cex"))

box(bty = "o")

dev.off()





