# ==========================================================
# = figure to compare imptuation for software combinations =
# ==========================================================

args <- commandArgs(TRUE) # args = c("../report/figure_imputation_software.pdf", "Myriad Pro", "../reportData/impute5.lr550.RData", "../reportData/beagle.lr550.RData", "../reportData/minimac4.lr550.RData", "../reportData/impute5.dly550.RData", "../reportData/beagle.dly550.RData", "../reportData/minimac4.dly550.RData", "../reportData/impute5.nodly550.RData", "../reportData/beagle.nodly550.RData", "../reportData/minimac4.nodly550.RData")
library("RColorBrewer")

# prepare file to plot
# ============================================================

file.width = 150
cairo_pdf(file = args[1], width = file.width/25.4, height = file.width/25.4, family = args[2])
layout(matrix(1:12, byrow = T, ncol = 4, nrow = 3), widths = c(0.3, 0.4, 0.4, 0.4))

# first design, LR only
# ============================================================
par(las = 1, tcl = -0.2, mai = c(0.1, 0.05, 0.02, 0.02)*file.width/25.4, ps = 7, lwd = 0.5, xpd = T)

plot(c(-100, 550), c(0, 2), type = "n", xlab = "", ylab = "", axes = FALSE)
polygon(c(0, 550, 550, 0), c(1.1, 1.1, 1.4, 1.4),  col = brewer.pal(8, "Dark2")[7], lty = 1)
text(275, 1.25, "550", col = "white", cex = 6/par("ps")/par("cex"))
text(-200, 1.25, "Reference", cex = 6/par("ps")/par("cex"))
text(275, 1.5, "Landrace", col = brewer.pal(8, "Dark2")[7], cex = 6/par("ps")/par("cex"))

polygon(c(0, 100, 100, 0), c(0.6, 0.6, 0.9, 0.9),  col = brewer.pal(8, "Dark2")[7], lty = 1)
text(-140, 0.75, "Target", cex = 6/par("ps")/par("cex"))
text(300, 0.75, "Landrace", col = brewer.pal(8, "Dark2")[7], cex = 6/par("ps")/par("cex"))
text(50, 0.75, "100", col = "white", cex = 6/par("ps")/par("cex"))

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# 1. concordance rate
# ============================================================
par(las = 1, tcl = -0.2, mai = c(0.12, 0.05, 0.02, 0.01)*file.width/25.4, ps = 7, lwd = 0.5, xpd = T)

plot(c(0, 0.5), c(0, 100), axes = FALSE, type = "n", xlab = "", ylab = "")
load(args[3])
impute5.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 1] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[1])
load(args[4])
beagle.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 1] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[2])
load(args[5])
minimac4.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 1] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[3])

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), at = seq(0, 0.5, 0.1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(1, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "Minor allele frequency (MAF)", mgp = c(0.9, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Concordance rate (%)", mgp = c(1.1, 0, 0), cex.lab = 7/par("ps")/par("cex"))

legend("bottomright", bty = "n", lwd = 1, seg.len = 1, y.intersp = 0.5, x.intersp = 0.5, legend = c(paste("SHAPEIT4/IMPUTE5: ", formatC(impute5.mean.accu[1]*100, format = "f", digits = 2), "%", sep = ""), paste("Beagle5.2/Beagle5.2: ", formatC(beagle.mean.accu[1]*100, format = "f", digits = 2), "%", sep = ""), paste("Eagle2.4/Minimac4: ", formatC(minimac4.mean.accu[1]*100, format = "f", digits = 2), "%", sep = "")), col = brewer.pal(8, "Dark2")[c(1:3)], text.col = brewer.pal(8, "Dark2")[c(1:3)], cex = 5/par("ps")/par("cex"))

text(0.25, 30, "Mean concordance rate", cex = 6/par("ps")/par("cex"))



# 2. NRC
# ============================================================
plot(c(0, 0.5), c(0, 100), axes = FALSE, type = "n", xlab = "", ylab = "")
load(args[3])
impute5.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 2] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[1])
load(args[4])
beagle.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 2] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[2])
load(args[5])
minimac4.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 2] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[3])

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), at = seq(0, 0.5, 0.1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(1, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "Minor allele frequency (MAF)", mgp = c(0.9, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Non-ref concordance rate (%)", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))

legend("bottomright", bty = "n", lwd = 1, seg.len = 1, y.intersp = 0.5, x.intersp = 0.5, legend = c(paste("SHAPEIT4/IMPUTE5: ", formatC(impute5.mean.accu[2]*100, format = "f", digits = 2), "%", sep = ""), paste("Beagle5.2/Beagle5.2: ", formatC(beagle.mean.accu[2]*100, format = "f", digits = 2), "%", sep = ""), paste("Eagle2.4/Minimac4: ", formatC(minimac4.mean.accu[2]*100, format = "f", digits = 2), "%", sep = "")), col = brewer.pal(8, "Dark2")[c(1:3)], text.col = brewer.pal(8, "Dark2")[c(1:3)], cex = 5/par("ps")/par("cex"))

text(0.25, 35, "Mean non-reference\nconcordance rate", cex = 6/par("ps")/par("cex"))

# 3. r2
# ============================================================
plot(c(0, 0.5), c(0, 1), axes = FALSE, type = "n", xlab = "", ylab = "")
load(args[3])
impute5.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 3], type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[1])
load(args[4])
beagle.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 3], type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[2])
load(args[5])
minimac4.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 3], type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[3])

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), at = seq(0, 0.5, 0.1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(1, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "Minor allele frequency (MAF)", mgp = c(0.9, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic(r^2), " (imputed versus observed)")), mgp = c(1.1, 0, 0), cex.lab = 7/par("ps")/par("cex"))

legend("bottomright", bty = "n", lwd = 1, seg.len = 1, y.intersp = 0.5, x.intersp = 0.5, legend = c(paste("SHAPEIT4/IMPUTE5: ", formatC(impute5.mean.accu[3], format = "f", digits = 2), "", sep = ""), paste("Beagle5.2/Beagle5.2: ", formatC(beagle.mean.accu[3], format = "f", digits = 2), "", sep = ""), paste("Eagle2.4/Minimac4: ", formatC(minimac4.mean.accu[3], format = "f", digits = 2), "", sep = "")), col = brewer.pal(8, "Dark2")[c(1:3)], text.col = brewer.pal(8, "Dark2")[c(1:3)], cex = 5/par("ps")/par("cex"))

text(0.33, 0.30, expression(paste("Mean ", italic(r^2))), cex = 6/par("ps")/par("cex"))

# second design, DLY
# ============================================================
par(las = 1, tcl = -0.2, mai = c(0.1, 0.05, 0.02, 0.02)*file.width/25.4, ps = 7, lwd = 0.5, xpd = T)

plot(c(-100, 550), c(0, 2), type = "n", xlab = "", ylab = "", axes = FALSE)
polygon(c(0, 250, 250, 0), c(1.1, 1.1, 1.4, 1.4),  col = brewer.pal(8, "Dark2")[7], lty = 1)
polygon(c(250, 400, 400, 250), c(1.1, 1.1, 1.4, 1.4),  col = brewer.pal(8, "Dark2")[3], lty = 1)
polygon(c(400, 550, 550, 400), c(1.1, 1.1, 1.4, 1.4),  col = brewer.pal(8, "Dark2")[8], lty = 1)
text(125, 1.25, "250", col = "white", cex = 6/par("ps")/par("cex"))
text(325, 1.25, "150", col = "white", cex = 6/par("ps")/par("cex"))
text(475, 1.25, "150", col = "white", cex = 6/par("ps")/par("cex"))

text(80, 1.5, "Landrace", srt = 45, pos = 4, cex = 6/par("ps")/par("cex"), col = brewer.pal(8, "Dark2")[7])
text(260, 1.5, "Duroc", srt = 45, pos = 4, cex = 6/par("ps")/par("cex"), col = brewer.pal(8, "Dark2")[3])
text(400, 1.5, "Yorkshire", srt = 45, pos = 4, cex = 6/par("ps")/par("cex"), col = brewer.pal(8, "Dark2")[8])

text(-200, 1.25, "Reference", cex = 6/par("ps")/par("cex"))

polygon(c(0, 100, 100, 0), c(0.6, 0.6, 0.9, 0.9),  col = brewer.pal(8, "Dark2")[7], lty = 1)
text(-140, 0.75, "Target", cex = 6/par("ps")/par("cex"))
text(300, 0.75, "Landrace", col = brewer.pal(8, "Dark2")[7], cex = 6/par("ps")/par("cex"))
text(50, 0.75, "100", col = "white", cex = 6/par("ps")/par("cex"))

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# 1. concordance rate
# ============================================================
par(las = 1, tcl = -0.2, mai = c(0.12, 0.05, 0.02, 0.01)*file.width/25.4, ps = 7, lwd = 0.5, xpd = T)

plot(c(0, 0.5), c(0, 100), axes = FALSE, type = "n", xlab = "", ylab = "")
load(args[6])
impute5.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 1] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[1])
load(args[7])
beagle.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 1] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[2])
load(args[8])
minimac4.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 1] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[3])

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), at = seq(0, 0.5, 0.1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(1, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "Minor allele frequency (MAF)", mgp = c(0.9, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Concordance rate (%)", mgp = c(1.1, 0, 0), cex.lab = 7/par("ps")/par("cex"))

legend("bottomright", bty = "n", lwd = 1, seg.len = 1, y.intersp = 0.5, x.intersp = 0.5, legend = c(paste("SHAPEIT4/IMPUTE5: ", formatC(impute5.mean.accu[1]*100, format = "f", digits = 2), "%", sep = ""), paste("Beagle5.2/Beagle5.2: ", formatC(beagle.mean.accu[1]*100, format = "f", digits = 2), "%", sep = ""), paste("Eagle2.4/Minimac4: ", formatC(minimac4.mean.accu[1]*100, format = "f", digits = 2), "%", sep = "")), col = brewer.pal(8, "Dark2")[c(1:3)], text.col = brewer.pal(8, "Dark2")[c(1:3)], cex = 5/par("ps")/par("cex"))

text(0.25, 30, "Mean concordance rate", cex = 6/par("ps")/par("cex"))



# 2. NRC
# ============================================================
plot(c(0, 0.5), c(0, 100), axes = FALSE, type = "n", xlab = "", ylab = "")
load(args[6])
impute5.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 2] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[1])
load(args[7])
beagle.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 2] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[2])
load(args[8])
minimac4.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 2] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[3])

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), at = seq(0, 0.5, 0.1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(1, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "Minor allele frequency (MAF)", mgp = c(0.9, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Non-ref concordance rate (%)", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))

legend("bottomright", bty = "n", lwd = 1, seg.len = 1, y.intersp = 0.5, x.intersp = 0.5, legend = c(paste("SHAPEIT4/IMPUTE5: ", formatC(impute5.mean.accu[2]*100, format = "f", digits = 2), "%", sep = ""), paste("Beagle5.2/Beagle5.2: ", formatC(beagle.mean.accu[2]*100, format = "f", digits = 2), "%", sep = ""), paste("Eagle2.4/Minimac4: ", formatC(minimac4.mean.accu[2]*100, format = "f", digits = 2), "%", sep = "")), col = brewer.pal(8, "Dark2")[c(1:3)], text.col = brewer.pal(8, "Dark2")[c(1:3)], cex = 5/par("ps")/par("cex"))

text(0.25, 35, "Mean non-reference\nconcordance rate", cex = 6/par("ps")/par("cex"))

# 3. r2
# ============================================================
plot(c(0, 0.5), c(0, 1), axes = FALSE, type = "n", xlab = "", ylab = "")
load(args[6])
impute5.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 3], type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[1])
load(args[7])
beagle.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 3], type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[2])
load(args[8])
minimac4.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 3], type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[3])

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), at = seq(0, 0.5, 0.1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(1, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "Minor allele frequency (MAF)", mgp = c(0.9, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic(r^2), " (imputed versus observed)")), mgp = c(1.1, 0, 0), cex.lab = 7/par("ps")/par("cex"))

legend("bottomright", bty = "n", lwd = 1, seg.len = 1, y.intersp = 0.5, x.intersp = 0.5, legend = c(paste("SHAPEIT4/IMPUTE5: ", formatC(impute5.mean.accu[3], format = "f", digits = 2), "", sep = ""), paste("Beagle5.2/Beagle5.2: ", formatC(beagle.mean.accu[3], format = "f", digits = 2), "", sep = ""), paste("Eagle2.4/Minimac4: ", formatC(minimac4.mean.accu[3], format = "f", digits = 2), "", sep = "")), col = brewer.pal(8, "Dark2")[c(1:3)], text.col = brewer.pal(8, "Dark2")[c(1:3)], cex = 5/par("ps")/par("cex"))

text(0.33, 0.30, expression(paste("Mean ", italic(r^2))), cex = 6/par("ps")/par("cex"))

# third design, LO
# ============================================================
par(las = 1, tcl = -0.2, mai = c(0.1, 0.05, 0.02, 0.02)*file.width/25.4, ps = 7, lwd = 0.5, xpd = T)

plot(c(-100, 550), c(0, 2), type = "n", xlab = "", ylab = "", axes = FALSE)
polygon(c(0, 250, 250, 0), c(1.1, 1.1, 1.4, 1.4),  col = brewer.pal(8, "Dark2")[7], lty = 1)
polygon(c(250, 550, 5500, 250), c(1.1, 1.1, 1.4, 1.4),  col = brewer.pal(8, "Dark2")[1], lty = 1)
text(125, 1.25, "250", col = "white", cex = 6/par("ps")/par("cex"))
text(400, 1.25, "300", col = "white", cex = 6/par("ps")/par("cex"))

text(80, 1.5, "Landrace", srt = 45, pos = 4, cex = 6/par("ps")/par("cex"), col = brewer.pal(8, "Dark2")[7])
text(400, 1.5, "Other", srt = 45, pos = 4, cex = 6/par("ps")/par("cex"), col = brewer.pal(8, "Dark2")[1])

text(-200, 1.25, "Reference", cex = 6/par("ps")/par("cex"))

polygon(c(0, 100, 100, 0), c(0.6, 0.6, 0.9, 0.9),  col = brewer.pal(8, "Dark2")[7], lty = 1)
text(-140, 0.75, "Target", cex = 6/par("ps")/par("cex"))
text(300, 0.75, "Landrace", col = brewer.pal(8, "Dark2")[7], cex = 6/par("ps")/par("cex"))
text(50, 0.75, "100", col = "white", cex = 6/par("ps")/par("cex"))

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# 1. concordance rate
# ============================================================
par(las = 1, tcl = -0.2, mai = c(0.12, 0.05, 0.02, 0.01)*file.width/25.4, ps = 7, lwd = 0.5, xpd = T)

plot(c(0, 0.5), c(0, 100), axes = FALSE, type = "n", xlab = "", ylab = "")
load(args[9])
impute5.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 1] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[1])
load(args[10])
beagle.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 1] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[2])
load(args[11])
minimac4.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 1] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[3])

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), at = seq(0, 0.5, 0.1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(1, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "Minor allele frequency (MAF)", mgp = c(0.9, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Concordance rate (%)", mgp = c(1.1, 0, 0), cex.lab = 7/par("ps")/par("cex"))

legend("bottomright", bty = "n", lwd = 1, seg.len = 1, y.intersp = 0.5, x.intersp = 0.5, legend = c(paste("SHAPEIT4/IMPUTE5: ", formatC(impute5.mean.accu[1]*100, format = "f", digits = 2), "%", sep = ""), paste("Beagle5.2/Beagle5.2: ", formatC(beagle.mean.accu[1]*100, format = "f", digits = 2), "%", sep = ""), paste("Eagle2.4/Minimac4: ", formatC(minimac4.mean.accu[1]*100, format = "f", digits = 2), "%", sep = "")), col = brewer.pal(8, "Dark2")[c(1:3)], text.col = brewer.pal(8, "Dark2")[c(1:3)], cex = 5/par("ps")/par("cex"))

text(0.25, 30, "Mean concordance rate", cex = 6/par("ps")/par("cex"))



# 2. NRC
# ============================================================
plot(c(0, 0.5), c(0, 100), axes = FALSE, type = "n", xlab = "", ylab = "")
load(args[9])
impute5.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 2] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[1])
load(args[10])
beagle.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 2] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[2])
load(args[11])
minimac4.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 2] * 100, type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[3])

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), at = seq(0, 0.5, 0.1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(1, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "Minor allele frequency (MAF)", mgp = c(0.9, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Non-ref concordance rate (%)", mgp = c(1.2, 0, 0), cex.lab = 7/par("ps")/par("cex"))

legend("bottomright", bty = "n", lwd = 1, seg.len = 1, y.intersp = 0.5, x.intersp = 0.5, legend = c(paste("SHAPEIT4/IMPUTE5: ", formatC(impute5.mean.accu[2]*100, format = "f", digits = 2), "%", sep = ""), paste("Beagle5.2/Beagle5.2: ", formatC(beagle.mean.accu[2]*100, format = "f", digits = 2), "%", sep = ""), paste("Eagle2.4/Minimac4: ", formatC(minimac4.mean.accu[2]*100, format = "f", digits = 2), "%", sep = "")), col = brewer.pal(8, "Dark2")[c(1:3)], text.col = brewer.pal(8, "Dark2")[c(1:3)], cex = 5/par("ps")/par("cex"))

text(0.25, 35, "Mean non-reference\nconcordance rate", cex = 6/par("ps")/par("cex"))

# 3. r2
# ============================================================
plot(c(0, 0.5), c(0, 1), axes = FALSE, type = "n", xlab = "", ylab = "")
load(args[9])
impute5.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 3], type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[1])
load(args[10])
beagle.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 3], type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[2])
load(args[11])
minimac4.mean.accu <- mean.accu
points(seq(0.005, 0.495, 0.01), accu[, 3], type = "l", lwd = 1, col = brewer.pal(8, "Dark2")[3])

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), at = seq(0, 0.5, 0.1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(1, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "Minor allele frequency (MAF)", mgp = c(0.9, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic(r^2), " (imputed versus observed)")), mgp = c(1.1, 0, 0), cex.lab = 7/par("ps")/par("cex"))

legend("bottomright", bty = "n", lwd = 1, seg.len = 1, y.intersp = 0.5, x.intersp = 0.5, legend = c(paste("SHAPEIT4/IMPUTE5: ", formatC(impute5.mean.accu[3], format = "f", digits = 2), "", sep = ""), paste("Beagle5.2/Beagle5.2: ", formatC(beagle.mean.accu[3], format = "f", digits = 2), "", sep = ""), paste("Eagle2.4/Minimac4: ", formatC(minimac4.mean.accu[3], format = "f", digits = 2), "", sep = "")), col = brewer.pal(8, "Dark2")[c(1:3)], text.col = brewer.pal(8, "Dark2")[c(1:3)], cex = 5/par("ps")/par("cex"))

text(0.33, 0.30, expression(paste("Mean ", italic(r^2))), cex = 6/par("ps")/par("cex"))


dev.off()


