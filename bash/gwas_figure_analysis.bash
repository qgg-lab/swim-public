# =============================
# = data for figures for gwas =
# =============================


# 1. bf, p value, make manhattan plot png
# get chip p values
# ============================================================

module purge
module load GCC/8.3.0 OpenMPI/3.1.4
module load R/4.0.2
Rscript -e 'impute.pval <- read.table("/mnt/research/qgg/share/dingr/swim_data/GWAS_data/S21_imputation_EBF.mlma", header = T, as.is = T); chr.len <- as.numeric(read.table("/mnt/research/qgg/resource/sscrofa11.1/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.fai", header = FALSE, as.is = TRUE)[1:18, 2]); names(chr.len) <- as.character(1:18); chr.cum <- c(0, cumsum(chr.len)[1:17]) + (0:17)*1000000; names(chr.cum) <- as.character(1:18); chr.col <- rep(c("grey50", "grey80"), 9); png(file = "bf.man.png", type = "cairo", width = 90/25.4 - 0.20*90/25.4/2, height = 90/25.4/2 - 0.22*90/25.4/2, units = "in", res = 1200); par(las = 1, tcl = -0.2, mai = c(0, 0, 0, 0), ps = 7, lwd = 0.5, xpd = T); plot(chr.cum[as.character(impute.pval[, 1])] + impute.pval[, 3], -log10(impute.pval[, 6]), bg = chr.col[impute.pval[, 1]], xaxs = "i", yaxs = "i", xlim = c(-2e7, chr.cum[18] + chr.len[18] + 2e7), ylim = c(-0.5, 14), axes = F, cex = 0.2, pch = 21, col = NULL); dev.off();'

Rscript -e 'chip.pval <- read.table("/mnt/research/qgg/share/dingr/swim_data/GWAS_data/S21_EBF.mlma", header = T, as.is = T); save(chip.pval, file = "bf.s21.pval.RData");'

# 2. gene models
# ============================================================

awk '$3 == "gene" && $1 == 1 && $4 >= 158500000 && $5 <= 162500000' /mnt/research/qgg/resource/sscrofa11.1/annot/Sus_scrofa.Sscrofa11.1.105.gff3 | perl -ne 'chomp $_; @line = split /\t/, $_; print join("\t", @line[(0, 3, 4, 6)]), "\t"; if ($line[8] =~ m/Name=(.*?);/) { print $1, "\n"; } else { print ".\n"; }' > chr1.mc4r.gene.bound

# 3. p value and r2 in the region
# ============================================================

module purge
module load GCC/8.3.0 OpenMPI/3.1.4
module load R/4.0.2
Rscript -e 'impute.pval <- read.table("/mnt/research/qgg/share/dingr/swim_data/GWAS_data/S21_EBF_imputation.chr1_158.5_162.5.mlma", header = T, as.is = T); r2 <- read.table("/mnt/research/qgg/share/dingr/swim_data/GWAS_data/EBF_topSNP_r2.ld", header = T, as.is = T); save(impute.pval, r2, file = "bf.s21.impute.pval.r2.RData");'

# 4. bl
# ============================================================

module purge
module load GCC/8.3.0 OpenMPI/3.1.4
module load R/4.0.2
Rscript -e 'impute.pval <- read.table("/mnt/research/qgg/share/dingr/swim_data/GWAS_data/w64_imputation_BL.mlma", header = T, as.is = T); chr.len <- as.numeric(read.table("/mnt/research/qgg/resource/sscrofa11.1/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.fai", header = FALSE, as.is = TRUE)[1:18, 2]); names(chr.len) <- as.character(1:18); chr.cum <- c(0, cumsum(chr.len)[1:17]) + (0:17)*1000000; names(chr.cum) <- as.character(1:18); chr.col <- rep(c("grey50", "grey80"), 9); png(file = "bl.man.png", type = "cairo", width = 90/25.4 - 0.20*90/25.4/2, height = 90/25.4/2 - 0.22*90/25.4/2, units = "in", res = 1200); par(las = 1, tcl = -0.2, mai = c(0, 0, 0, 0), ps = 7, lwd = 0.5, xpd = T); plot(chr.cum[as.character(impute.pval[, 1])] + impute.pval[, 3], -log10(impute.pval[, 6]), bg = chr.col[impute.pval[, 1]], xaxs = "i", yaxs = "i", xlim = c(-2e7, chr.cum[18] + chr.len[18] + 2e7), ylim = c(-1.4, 40), axes = F, cex = 0.2, pch = 21, col = NULL); dev.off();'

Rscript -e 'chip.pval <- read.table("/mnt/research/qgg/share/dingr/swim_data/GWAS_data/W64_BL.mlma", header = T, as.is = T); save(chip.pval, file = "bl.w64.pval.RData");'

awk '$3 == "gene" && $1 == 17 && $4 >= 15300000 && $5 <= 16300000' /mnt/research/qgg/resource/sscrofa11.1/annot/Sus_scrofa.Sscrofa11.1.105.gff3 | perl -ne 'chomp $_; @line = split /\t/, $_; print join("\t", @line[(0, 3, 4, 6)]), "\t"; if ($line[8] =~ m/Name=(.*?);/) { print $1, "\n"; } else { print ".\n"; }' > chr17.bmp2.gene.bound

Rscript -e 'impute.pval <- read.table("/mnt/research/qgg/share/dingr/swim_data/GWAS_data/W64_BL_imputation.chr17_15.3_16.3.mlma", header = T, as.is = T); r2 <- read.table("/mnt/research/qgg/share/dingr/swim_data/GWAS_data/BL_topSNP_r2.ld", header = T, as.is = T); save(impute.pval, r2, file = "bl.w64.impute.pval.r2.RData");'

# 5. bw
# ============================================================

module purge
module load GCC/8.3.0 OpenMPI/3.1.4
module load R/4.0.2
Rscript -e 'impute.pval <- read.table("/mnt/research/qgg/share/dingr/swim_data/GWAS_data/S22_imputation_BW0.mlma", header = T, as.is = T); chr.len <- as.numeric(read.table("/mnt/research/qgg/resource/sscrofa11.1/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.fai", header = FALSE, as.is = TRUE)[1:18, 2]); names(chr.len) <- as.character(1:18); chr.cum <- c(0, cumsum(chr.len)[1:17]) + (0:17)*1000000; names(chr.cum) <- as.character(1:18); chr.col <- rep(c("grey50", "grey80"), 9); png(file = "bw.man.png", type = "cairo", width = 90/25.4 - 0.20*90/25.4/2, height = 90/25.4/2 - 0.22*90/25.4/2, units = "in", res = 1200); par(las = 1, tcl = -0.2, mai = c(0, 0, 0, 0), ps = 7, lwd = 0.5, xpd = T); plot(chr.cum[as.character(impute.pval[, 1])] + impute.pval[, 3], -log10(impute.pval[, 6]), bg = chr.col[impute.pval[, 1]], xaxs = "i", yaxs = "i", xlim = c(-2e7, chr.cum[18] + chr.len[18] + 2e7), ylim = c(-0.33, 8), axes = F, cex = 0.2, pch = 21, col = NULL); dev.off();'

Rscript -e 'chip.pval <- read.table("/mnt/research/qgg/share/dingr/swim_data/GWAS_data/S22_BW0.mlma", header = T, as.is = T); save(chip.pval, file = "bw.s22.pval.RData");'

awk '$3 == "gene" && $1 == 2 && $4 >= 7500000 && $5 <= 8650000' /mnt/research/qgg/resource/sscrofa11.1/annot/Sus_scrofa.Sscrofa11.1.105.gff3 | perl -ne 'chomp $_; @line = split /\t/, $_; print join("\t", @line[(0, 3, 4, 6)]), "\t"; if ($line[8] =~ m/Name=(.*?);/) { print $1, "\n"; } else { print ".\n"; }' > chr2.slc.gene.bound

Rscript -e 'impute.pval <- read.table("/mnt/research/qgg/share/dingr/swim_data/GWAS_data/S22_BW0_imputation.chr2_7.50_8.65.mlma", header = T, as.is = T); r2 <- read.table("/mnt/research/qgg/share/dingr/swim_data/GWAS_data/BW0_topSNP_r2.ld", header = T, as.is = T); save(impute.pval, r2, file = "bw.s22.impute.pval.r2.RData");'
