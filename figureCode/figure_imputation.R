# ===========================================
# = figure for PCA, genetic structure, etc. =
# ===========================================

args <- commandArgs(TRUE) # args = c("../report/figure_imputation.pdf", "Myriad Pro", "../reportData/all.maf.hist.RData", "../reportData/chip.maf.hist.RData", "impute5.all.RData", "../reportData/impute5.lr550.RData", "../reportData/impute5.dly550.RData",  "../reportData/impute5.nodly550.RData", "pharp.RData")
library("RColorBrewer")

# prepare file to plot
# ============================================================

file.width = 135
cairo_pdf(file = args[1], width = file.width/25.4, height = file.width/25.4*0.6, family = args[2])
layout(matrix(c(1, 2, 2, 3, 4, 5), byrow = T, ncol = 3, nrow = 2))

# histogram of maf for panel and snp chips
# ============================================================

load(args[3])
all.maf.prop <- maf.counts/sum(maf.counts)
load(args[4])
chip.maf.prop <- maf.counts/sum(maf.counts)

par(las = 1, tcl = -0.2, mai = c(0.06, 0.07, 0.02, 0.02)*file.width/25.4, ps = 7, lwd = 0.5, xpd = T)
plot(c(0, 0.5), c(0, 12), type = "n", xlab = "", ylab = "", axes = FALSE)

polygon(rep(seq(0, 0.5, 0.005), each = 2), c(0, rep(all.maf.prop*100, each = 2), 0), col = "#0000FF20")
polygon(rep(seq(0, 0.5, 0.005), each = 2), c(0, rep(chip.maf.prop*100, each = 2), 0), col = "#FF000020")


axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "Minor allele frequency (MAF)", mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Proportion of variants (%)", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))

text(0, 10, "All variants (34 M)", pos = 4, col = "#0000FFA0", cex = 7/par("ps")/par("cex"))
text(0.08, 2.5, "SNP array variants (40 K)", pos = 4, col = "#FF0000A0", cex = 7/par("ps")/par("cex"))

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# schematic diagram of the expeirmental design
# ============================================================

par(las = 1, tcl = -0.2, mai = c(0.05, 0.05, 0.02, 0.02)*file.width/25.4, ps = 7, lwd = 0.5, xpd = T)
plot(c(-200, 2700), c(-0.2, 6.5), type = "n", xlab = "", ylab = "", axes = FALSE)

# in reference - 100, 485 duroc, 551 landrace, 543 yorkshire, 580 others

polygon(c(0, 1006, 1006, 0), c(0, 0, 0.8, 0.8), lty = 2)
polygon(c(0, 1006, 1006, 0), c(0, 0, 0.8, 0.8), border = NA, col = brewer.pal(8, "Dark2")[6])
polygon(c(2800, 2900, 2900, 2800) - 300, c(0, 0, 0.8, 0.8))
polygon(c(2800, 2900, 2900, 2800) - 300, c(0, 0, 0.8, 0.8), border = NA, col = brewer.pal(8, "Dark2")[7])


polygon(c(0, 550, 550, 0), c(0, 0, 0.8, 0.8) + 1)
polygon(c(0, 250, 250, 0), c(0, 0, 0.8, 0.8) + 1, border = NA, col = brewer.pal(8, "Dark2")[7])
polygon(c(250, 550, 550, 250), c(0, 0, 0.8, 0.8) + 1, border = NA, col = brewer.pal(8, "Dark2")[1])
polygon(c(2800, 2900, 2900, 2800) - 300, c(0, 0, 0.8, 0.8)+1)
polygon(c(2800, 2900, 2900, 2800) - 300, c(0, 0, 0.8, 0.8)+1, border = NA, col = brewer.pal(8, "Dark2")[7])

polygon(c(0, 550, 550, 0), c(0, 0, 0.8, 0.8) + 2)
polygon(c(0, 250, 250, 0), c(0, 0, 0.8, 0.8) + 2, border = NA, col = brewer.pal(8, "Dark2")[7])
polygon(c(250, 400, 400, 250), c(0, 0, 0.8, 0.8) + 2, border = NA, col = brewer.pal(8, "Dark2")[3])
polygon(c(400, 550, 550, 400), c(0, 0, 0.8, 0.8) + 2, border = NA, col = brewer.pal(8, "Dark2")[8])
polygon(c(2800, 2900, 2900, 2800) - 300, c(0, 0, 0.8, 0.8)+2)
polygon(c(2800, 2900, 2900, 2800) - 300, c(0, 0, 0.8, 0.8)+2, border = NA, col = brewer.pal(8, "Dark2")[7])

polygon(c(0, 550, 550, 0), c(0, 0, 0.8, 0.8) + 3)
polygon(c(0, 550, 550, 0), c(0, 0, 0.8, 0.8) + 3, border = NA, col = brewer.pal(8, "Dark2")[7])
polygon(c(2800, 2900, 2900, 2800) - 300, c(0, 0, 0.8, 0.8)+3)
polygon(c(2800, 2900, 2900, 2800) - 300, c(0, 0, 0.8, 0.8)+3, border = NA, col = brewer.pal(8, "Dark2")[7])


polygon(c(0, 2159, 2159, 0), c(0, 0, 0.8, 0.8) + 4)
polygon(c(0, 551, 551, 0), c(0, 0, 0.8, 0.8) + 4, border = NA, col = brewer.pal(8, "Dark2")[7])
polygon(c(551, 1036, 1036, 551), c(0, 0, 0.8, 0.8) + 4, border = NA, col = brewer.pal(8, "Dark2")[3])
polygon(c(1036, 1579, 1579, 1036), c(0, 0, 0.8, 0.8) + 4, border = NA, col = brewer.pal(8, "Dark2")[8])
polygon(c(1579, 2159, 2159, 1579), c(0, 0, 0.8, 0.8) + 4, border = NA, col = brewer.pal(8, "Dark2")[1])

polygon(c(2800, 2900, 2900, 2800) - 300, c(0, 0, 0.8, 0.8)+4)
polygon(c(2800, 2900, 2900, 2800) - 300, c(0, 0, 0.8, 0.8)+4, border = NA, col = brewer.pal(8, "Dark2")[7])


text(276, 4.7, "Landrace (L)", pos = 3, cex = 6/par("ps")/par("cex"), col = brewer.pal(8, "Dark2")[7])
text(800, 4.7, "Duroc (D)", pos = 3, cex = 6/par("ps")/par("cex"), col = brewer.pal(8, "Dark2")[3])
text(1320, 4.7, "Yorkshire (Y)", pos = 3, cex = 6/par("ps")/par("cex"), col = brewer.pal(8, "Dark2")[8])
text(1860, 4.7, "Other (O)", pos = 3, cex = 6/par("ps")/par("cex"), col = brewer.pal(8, "Dark2")[1])

text(551/2, 4.4, "551", col = "white", cex = 5/par("ps")/par("cex"))
text(551 + 485/2, 4.4, "485", col = "white", cex = 5/par("ps")/par("cex"))
text(1036 + 543/2, 4.4, "543", col = "white", cex = 5/par("ps")/par("cex"))
text(1579 + 580/2, 4.4, "580", col = "white", cex = 5/par("ps")/par("cex"))

text(550/2, 3.4, "550", col = "white", cex = 5/par("ps")/par("cex"))

text(250/2, 2.4, "250", col = "white", cex = 5/par("ps")/par("cex"))
text(250 + 150/2, 2.4, "150", col = "white", cex = 5/par("ps")/par("cex"))
text(400 + 150/2, 2.4, "150", col = "white", cex = 5/par("ps")/par("cex"))

text(250/2, 1.4, "250", col = "white", cex = 5/par("ps")/par("cex"))
text(250 + 300/2, 1.4, "300", col = "white", cex = 5/par("ps")/par("cex"))

text(1006/2, 0.4, "1006", col = "white", cex = 5/par("ps")/par("cex"))

text(0, 4.4, "All", pos = 2, cex = 6/par("ps")/par("cex"), col = brewer.pal(8, "Dark2")[1])
text(0, 3.4, "L", pos = 2, cex = 6/par("ps")/par("cex"), col = brewer.pal(8, "Dark2")[2])
text(0, 2.4, "DLY", pos = 2, cex = 6/par("ps")/par("cex"), col = brewer.pal(8, "Dark2")[3])
text(0, 1.4, "LO", pos = 2, cex = 6/par("ps")/par("cex"), col = brewer.pal(8, "Dark2")[4])
text(0, 0.4, "PHARP", pos = 2, cex = 6/par("ps")/par("cex"), col = brewer.pal(8, "Dark2")[6])

text(2550, 0.4, "100", col = "white", cex = 5/par("ps")/par("cex"))
text(2550, 1.4, "100", col = "white", cex = 5/par("ps")/par("cex"))
text(2550, 2.4, "100", col = "white", cex = 5/par("ps")/par("cex"))
text(2550, 3.4, "100", col = "white", cex = 5/par("ps")/par("cex"))
text(2550, 4.4, "100", col = "white", cex = 5/par("ps")/par("cex"))


axis(side = 1, at = c(0, 2159), label = c("", ""), lwd = 0.5, mgp = c(0, 0, -0.1), cex.axis = 7/par("ps")/par("cex"), tcl = 0.2)

axis(side = 1, at = c(2500, 2600), label = c("", ""), lwd = 0.5, mgp = c(0, 0, -0.1), cex.axis = 7/par("ps")/par("cex"), tcl = 0.2)

text(2550, 4.7, "Landrace (L)", pos = 3, cex = 6/par("ps")/par("cex"), col = brewer.pal(8, "Dark2")[7])
text(2550, -1.05, "Target", cex = 6/par("ps")/par("cex"))
text(2159/2, -1.05, "Haplotype reference panel", cex = 6/par("ps")/par("cex"))

text(grconvertX(0.05 + file.width/3/25.4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


# 3. concordance rate
# ============================================================
par(las = 1, tcl = -0.2, mai = c(0.06, 0.07, 0.02, 0.02)*file.width/25.4, ps = 7, lwd = 0.5, xpd = T)

# read accuracy data
load(args[5])
all <- accu; all.mean.accu <- mean.accu
load(args[6])
landrace <- accu; landrace.mean.accu <- mean.accu
load(args[7])
dly <- accu; dly.mean.accu <- mean.accu
load(args[8])
lo <- accu; lo.mean.accu <- mean.accu
load(args[9])
pharp <- accu; pharp.mean.accu <- mean.accu


plot(c(0, 0.5), c(0, 100), axes = FALSE, type = "n", xlab = "", ylab = "")
points(seq(0.005, 0.495, 0.01), all[, 1] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[1])
points(seq(0.005, 0.495, 0.01), landrace[, 1] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[2])
points(seq(0.005, 0.495, 0.01), dly[, 1] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[3])
points(seq(0.005, 0.495, 0.01), lo[, 1] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[4])
points(seq(0.005, 0.495, 0.01), pharp[, 1] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[6])

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), at = seq(0, 0.5, 0.1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "Minor allele frequency (MAF)", mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Concordance rate (%)", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))

legend("bottomright", bty = "n", lwd = 1, y.intersp = 0.5, x.intersp = 0.5, legend = c(paste("All: ", formatC(all.mean.accu[1]*100, format = "f", digits = 2), "%", sep = ""), paste("L: ", formatC(landrace.mean.accu[1]*100, format = "f", digits = 2), "%", sep = ""), paste("DLY: ", formatC(dly.mean.accu[1]*100, format = "f", digits = 2), "%", sep = ""), paste("LO: ", formatC(lo.mean.accu[1]*100, format = "f", digits = 2), "%", sep = ""), paste("PHARP: ", formatC(pharp.mean.accu[1]*100, format = "f", digits = 2), "%", sep = "")), col = brewer.pal(8, "Dark2")[c(1:4, 6)], text.col = brewer.pal(8, "Dark2")[c(1:4, 6)], cex = 6/par("ps")/par("cex"))

text(0.3, 50, "Mean concordance rate", cex = 6/par("ps")/par("cex"))

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# 4. NRC
# ============================================================
plot(c(0, 0.5), c(0, 100), axes = FALSE, type = "n", xlab = "", ylab = "")
points(seq(0.005, 0.495, 0.01), all[, 2] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[1])
points(seq(0.005, 0.495, 0.01), landrace[, 2] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[2])
points(seq(0.005, 0.495, 0.01), dly[, 2] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[3])
points(seq(0.005, 0.495, 0.01), lo[, 2] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[4])
points(seq(0.005, 0.495, 0.01), pharp[, 2] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[6])

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), at = seq(0, 0.5, 0.1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "Minor allele frequency (MAF)", mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Non-ref concordance rate (%)", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))

legend("bottomright", bty = "n", lwd = 1, y.intersp = 0.5, x.intersp = 0.5, legend = c(paste("All: ", formatC(all.mean.accu[2]*100, format = "f", digits = 2), "%", sep = ""), paste("L: ", formatC(landrace.mean.accu[2]*100, format = "f", digits = 2), "%", sep = ""), paste("DLY: ", formatC(dly.mean.accu[2]*100, format = "f", digits = 2), "%", sep = ""), paste("LO: ", formatC(lo.mean.accu[2]*100, format = "f", digits = 2), "%", sep = ""), paste("PHARP: ", formatC(pharp.mean.accu[2]*100, format = "f", digits = 2), "%", sep = "")), col = brewer.pal(8, "Dark2")[c(1:4, 6)], text.col = brewer.pal(8, "Dark2")[c(1:4, 6)], cex = 6/par("ps")/par("cex"))

text(0.3, 55, "Mean non-reference\nconcordance rate", cex = 6/par("ps")/par("cex"))

text(grconvertX(0.05 + file.width/25.4/3, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# 4. r2
# ============================================================
plot(c(0, 0.5), c(0, 1), axes = FALSE, type = "n", xlab = "", ylab = "")
points(seq(0.005, 0.495, 0.01), all[, 3], type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[1])
points(seq(0.005, 0.495, 0.01), landrace[, 3], type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[2])
points(seq(0.005, 0.495, 0.01), dly[, 3], type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[3])
points(seq(0.005, 0.495, 0.01), lo[, 3], type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[4])
points(seq(0.005, 0.495, 0.01), pharp[, 3], type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[6])

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), at = seq(0, 0.5, 0.1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "Minor allele frequency (MAF)", mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic(r^2), " (imputed versus observed)")), mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))

legend("bottomright", bty = "n", lwd = 1, y.intersp = 0.5, x.intersp = 0.5, legend = c(paste("All: ", formatC(all.mean.accu[3], format = "f", digits = 2), "", sep = ""), paste("L: ", formatC(landrace.mean.accu[3], format = "f", digits = 2), "", sep = ""), paste("DLY: ", formatC(dly.mean.accu[3], format = "f", digits = 2), "", sep = ""), paste("LO: ", formatC(lo.mean.accu[3], format = "f", digits = 2), "", sep = ""), paste("PHARP: ", formatC(pharp.mean.accu[3], format = "f", digits = 2), "", sep = "")), col = brewer.pal(8, "Dark2")[c(1:4, 6)], text.col = brewer.pal(8, "Dark2")[c(1:4, 6)], cex = 6/par("ps")/par("cex"))

text(0.3, 0.5, expression(paste("Mean ", italic(r^2))), cex = 6/par("ps")/par("cex"))

text(grconvertX(0.05 + file.width/25.4/3*2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("e")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

dev.off()


