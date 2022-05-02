# ================================
# = figure for backfat thickness =
# ================================

args <- commandArgs(TRUE) # args = c("../report/figure_bf.pdf", "Myriad Pro", "../reportData/bl.w64.pval.RData", "../reportData/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.fai", "../reportData/bl.man.png", "../reportData/bl.w64.impute.pval.r2.RData", "../reportData/chr17.bmp2.gene.bound", "../reportData/BL_topSNP_phe.csv")
library("RColorBrewer")
library("png")

# prepare file to plot
# ============================================================

file.width = 90
cairo_pdf(file = args[1], width = file.width/25.4, height = file.width/25.4*1.25, family = args[2])
layout(matrix(c(1, 1, 1, 2, 4, 4, 2, 5, 6, 2, 7, 8, 3, 7, 8), byrow = T, ncol = 3, nrow = 5), heights = c(0.5, 0.375, 0.1875, 0.0375, 0.15), widths = c(0.6, 0.2, 0.2))

# layout(matrix(c(1, 1, 2, 4, 2, 5, 3, 5), byrow = T, ncol = 2, nrow = 4), heights = c(0.5, 0.375, 0.225, 0.15), widths = c(0.6, 0.4))


# manhattan plot
# ============================================================

par(las = 1, tcl = -0.2, mai = c(0.12, 0.15, 0.1, 0.05)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = T)

load(args[3])

chr.len <- as.numeric(read.table(args[4], header = FALSE, as.is = TRUE)[1:18, 2])
names(chr.len) <- as.character(1:18)
chr.cum <- c(0, cumsum(chr.len)[1:17]) + (0:17)*1000000
names(chr.cum) <- as.character(1:18)
chr.col <- rep(c(brewer.pal(9, "Blues")[4], brewer.pal(9, "Blues")[8]), 9)

plot(c(0, chr.cum[18] + chr.len[18]), c(0, 40), type = "n", axes = FALSE, xaxs = "i", yaxs = "i", xlim = c(-2e7, chr.cum[18] + chr.len[18] + 2e7), ylim = c(-1.4, 40))
png.plot <- readPNG(args[5])
rasterImage(png.plot, -2e7, -1.4, chr.cum[18] + chr.len[18] + 2e7, 40)
points(chr.cum[as.character(chip.pval[, 1])] + chip.pval[, 3], -log10(chip.pval[, 9]), bg = chr.col[chip.pval[, 1]], pch = 21, cex = 0.1, lwd = 0.2, col = NULL)
axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), at = chr.cum + chr.len/2, label = 1:18, cex.axis = 5/par("ps")/par("cex"))
axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), at = chr.cum[c(11, 17)] + chr.len[c(11, 17)]/2, label = c(11, 17), cex.axis = 5/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "Chromosome", mgp = c(0.7, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste("-log", {}[10], "(", italic(P), "-value)")), mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))

points(c(0.1e9, 0.15e9), c(37, 37), pch = 21, cex = 0.5, bg = c(brewer.pal(9, "Blues")[4], brewer.pal(9, "Blues")[8]), col = NULL)
text(0.14e9, 37, "Affymetrix 55K", pos = 4, cex = 6/par("ps")/par("cex"))

points(c(0.1e9, 0.15e9), c(40, 40), pch = 21, cex = 0.5, bg = c("grey50", "grey80"), col = NULL)
text(0.14e9, 40, "Imputed", pos = 4, cex = 6/par("ps")/par("cex"))

# points(chr.cum[1] + chip.pval[chip.pval[, 2] == "12784636", 3], -log10(chip.pval[chip.pval[, 2] == "12784636", 9]), cex = 0.3)
# points(chr.cum[1] + 161511936, -log10(2.98122e-13), cex = 0.3)

# arrows(2.5e8, 11.5, chr.cum[1] + chip.pval[chip.pval[, 2] == "12784636", 3] + 10000000, -log10(chip.pval[chip.pval[, 2] == "12784636", 9]), length = 0.03)
# text(2.2e8, 11.2, expression(paste("Array ", italic(MC4R), " SNP")), pos = 4, cex = 6/par("ps")/par("cex"))
# arrows(2.5e8, 13, chr.cum[1] + 161511936 + 10000000, -log10(2.98122e-13), length = 0.03)
# text(2.2e8, 13, expression(paste("Imputed SNP chr1:161511936:T>C")), pos = 4, cex = 6/par("ps")/par("cex"))


text(chr.cum[8], 40, expression(paste(italic(h), {}["array"]^2, " = 0.32 (0.04)")), cex = 6/par("ps")/par("cex"), pos = 4)
text(chr.cum[8], 36, expression(paste(italic(h), {}["imputed"]^2, " = 0.34 (0.04)")), cex = 6/par("ps")/par("cex"), pos = 4)

box(bty = "l")

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# local plot
# ============================================================

load(args[6])

par(las = 1, tcl = -0.2, mai = c(0.01, 0.15, 0.1, 0.05)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = T)

impute.pval <- impute.pval[nrow(impute.pval):1, ]

plot(c(15300000, 16300000), c(-1.4, 40), type = "n", axes = FALSE, xaxs = "i", yaxs = "i", xlim = c(15250000, 16350000), ylim = c(-1.4, 40))
r2.col <- brewer.pal(9, "Blues")[c(2, 3, 5, 7, 9)]

rect(xleft = rep(15300000, length(r2.col)),
     ybottom = seq(length = length(r2.col), from = 30, to = 38)
               - (38-30)/(length(r2.col) - 1),
     xright = rep(15400000, length(r2.col)),
     ytop = seq(length = length(r2.col), from = 30, to = 38),
     col = r2.col)
text(15350000, 39.5, expression(paste(italic(r), {}^2)), cex = 6/par("ps")/par("cex"))
text(rep(15375000, 6),
     c(seq(length = length(r2.col), from = 30, to = 38) - (38-30)/(length(r2.col) - 1),
       seq(length = length(r2.col), from = 30, to = 38)[5]),
       formatC(seq(0, 1, 0.2), format = "f", digits = 1), pos = 4)

points(impute.pval[, 3], -log10(impute.pval[, 9]), pch = 21, bg = r2.col[as.numeric(cut(r2[match(impute.pval[, 2], r2[, 6]), 7], breaks = seq(0, 1, 0.2), include.lowest = TRUE))], col = NULL, cex = 0.8)
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
title(ylab = expression(paste("-log", {}[10], "(", italic(P), "-value)")), mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))

points(impute.pval[impute.pval[, 3] == 15643342, 3], -log10(impute.pval[impute.pval[, 3] == 15643342, 9]), cex = 1.2)
points(impute.pval[impute.pval[, 3] == 15626425, 3], -log10(impute.pval[impute.pval[, 3] == 15626425, 9]), cex = 1.2)
arrows(impute.pval[impute.pval[, 3] == 15643342, 3] + 80000, 39, impute.pval[impute.pval[, 3] == 15643342, 3] + 20000, -log10(impute.pval[impute.pval[, 3] == 15643342, 9]), length = 0.03)
text(impute.pval[impute.pval[, 3] == 15643342, 3] + 50000, 39, expression(paste("Imputed chr17:15643342:C>T")), pos = 4, cex = 6/par("ps")/par("cex"))

arrows(impute.pval[impute.pval[, 3] == 15626425, 3] + 80000, 36, impute.pval[impute.pval[, 3] == 15626425, 3] + 20000, -log10(impute.pval[impute.pval[, 3] == 15626425, 9]), length = 0.03)
text(impute.pval[impute.pval[, 3] == 15626425, 3] + 50000, 36, expression(paste("Imputed chr17:15626425:C>T")), pos = 4, cex = 6/par("ps")/par("cex"))

# mark chip SNPs

chip.impute.pval <- impute.pval[na.omit(match(paste(chip.pval[, 1], chip.pval[, 3], sep = "_"), paste(impute.pval[, 1], impute.pval[, 3], sep = "_"))), ]
points(chip.impute.pval[, 3], -log10(chip.impute.pval[, 9]), cex = 0.8, pch = 4)

legend("right", pch = c(21, 4), pt.bg = c(brewer.pal(9, "Blues")[9], "black"), legend = c("Imputed", "Array"), bty = "n", x.intersp = 0.5, y.intersp = 0.8, cex = 6/par("ps")/par("cex"))


box(bty = "l")
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# genes
# ============================================================

genes <- read.table(args[7], header = FALSE, as.is = TRUE)

par(las = 1, tcl = -0.2, mai = c(0.15, 0.15, 0.01, 0.05)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = T)

plot(c(15300000, 16300000), c(0, 2), type = "n", axes = FALSE, xaxs = "i", yaxs = "i", xlim = c(15250000, 16350000), ylim = c(0, 2))
rect(xleft = genes[, 2],
     ybottom = rep(1.2, nrow(genes)),
     xright = genes[, 3],
     ytop = rep(1.9, nrow(genes)),
     col = ifelse(genes[, 4] == "+", brewer.pal(9, "Reds")[9], brewer.pal(9, "Blues")[9]), border = NA)

text((genes[genes[, 5] == "BMP2", 2] + genes[genes[, 5] == "BMP2", 3])/2, 0.5, ">BMP2", cex = 6/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])
axis(side = 1, lwd = 0.5, at = seq(15400000, 16200000, 400000), label = seq(15.4, 16.2, 0.4), mgp = c(0.8, 0, 0), cex.axis = 6/par("ps")/par("cex"))
title(xlab = "Position on Chromosome 17 (Mb)", mgp = c(0.7, 0, 0), cex.lab = 7/par("ps")/par("cex"))
#
box(bty = "o")

# phenotypes
# ============================================================

par(las = 1, tcl = -0.2, mai = c(0.15, 0.2, 0.1, 0.05)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = T)

bl.pheno <- read.csv(args[8], header = T, as.is = T)

plot(c(0.5, 3), c(105, 145), type = "n", axes = FALSE, xaxs = "i", yaxs = "i", xlim = c(0.5, 3.5), ylim = c(105, 145), xlab = "", ylab = "")
points(as.numeric(factor(bl.pheno[, 3])) + runif(nrow(bl.pheno), -0.4, 0.4), bl.pheno[, 2], cex = 0.2, pch = 3, lwd = 0.1)
this.box <- boxplot(bl.pheno[, 2] ~ bl.pheno[, 3], range = 0, plot = F)
bxp(this.box, add = T, axes = FALSE, boxfill = NA, medlwd = 2, boxcol = brewer.pal(9, "Blues")[7], whiskcol = brewer.pal(9, "Blues")[7], staplecol = brewer.pal(9, "Blues")[7], medcol = brewer.pal(9, "Blues")[7], whisklty = 1, lwd = 0.75)
axis(side = 1, lwd = 0.5, at = 1:3, label = c("C/C", "C/T", "T/T"), mgp = c(0.8, 0, 0), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))

title(xlab = "Genotype", mgp = c(0.7, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Body length (cm)", mgp = c(1.4, 0, 0), cex.lab = 7/par("ps")/par("cex"))

box(bty = "l")

text(grconvertX(0.05 + file.width/25.4*0.6, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# allele frequency
# ============================================================

par(las = 1, tcl = -0.2, mai = c(0.02, 0, 0, 0)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = T)

pie(c(0.92,0.08), col = c(brewer.pal(9, "Blues")[8], brewer.pal(8, "Dark2")[1]), border = NA, radius = 0.8, labels = "")

text(-0.2, 0.2, "C:92%", col = "white", cex = 6/par("ps")/par("cex"))
text(0.5, -0.1, "T:8%", col = "white", cex = 6/par("ps")/par("cex"))

text(grconvertX(0.05 + file.width/25.4*0.6, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)
text(0, -0.95, "Landrace", cex = 7/par("ps")/par("cex"))

pie(c(0.22,0.78), col = c(brewer.pal(9, "Blues")[8], brewer.pal(8, "Dark2")[1]), border = NA, radius = 0.8, labels = "")

text(0.37, 0.2, "C:22%", col = "white", cex = 6/par("ps")/par("cex"))
text(-0.2, -0.2, "T:78%", col = "white", cex = 6/par("ps")/par("cex"))
text(0, -0.95, "Duroc", cex = 7/par("ps")/par("cex"))

pie(c(0.25,0.75), col = c(brewer.pal(9, "Blues")[8], brewer.pal(8, "Dark2")[1]), border = NA, radius = 0.8, labels = "")

text(0.37, 0.2, "C:25%", col = "white", cex = 6/par("ps")/par("cex"))
text(-0.2, -0.2, "T:75%", col = "white", cex = 6/par("ps")/par("cex"))
text(0, -0.95, "Yorkshire", cex = 7/par("ps")/par("cex"))

pie(c(0.23,0.77), col = c(brewer.pal(9, "Blues")[8], brewer.pal(8, "Dark2")[1]), border = NA, radius = 0.8, labels = "")

text(0.37, 0.2, "C:23%", col = "white", cex = 6/par("ps")/par("cex"))
text(-0.2, -0.2, "T:77%", col = "white", cex = 6/par("ps")/par("cex"))
text(0, -0.95, "Chinese breeds", cex = 7/par("ps")/par("cex"))

dev.off()





