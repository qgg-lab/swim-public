# ===============================================
# = make file for generating tables and figures =
# ===============================================

all: ../report/figure_genetics.pdf ../report/figure_imputation.pdf ../report/figure_imputation_software.pdf ../report/figure_bf.pdf ../report/figure_bl.pdf ../report/figure_bw.pdf
	
FONT="Myriad Pro"

# figure for PCA
# ============================================================

../report/figure_genetics.pdf: figure_genetics.R ../reportData/all.ids.random.snp.count ../reportData/figurePCA.RData ../reportData/Fig_muti2.Landrace ../reportData/Fig_muti2.Yorkshire ../reportData/Fig_muti2.Duroc ../reportData/Fig_muti2.Meishan ../reportData/Fig_muti2.Tibetan ../reportData/Fig_muti2.WB_Asia ../reportData/Fig_muti2.WB_European
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figure_genetics.R ../reportData/all.ids.random.snp.count ../reportData/Fig_muti.Landrace ../reportData/Fig_muti.Yorkshire ../reportData/Fig_muti.Duroc ../reportData/Fig_muti.Meishan ../reportData/Fig_muti.Tibetan ../reportData/Fig_muti.WB_Asia ../reportData/Fig_muti.WB_European ../reportData/figurePCA.RData ../reportData/SWIM_BreedInfo.txt ../reportData/select.all.pruned.fam ../reportData/select.all.pruned.2.Q ../reportData/select.all.pruned.4.Q ../reportData/select.all.pruned.6.Q ../reportData/breed.code.csv ../report/figure_genetics.pdf $(FONT)

../report/figure_imputation.pdf: figure_imputation.R ../reportData/all.maf.hist.RData ../reportData/chip.maf.hist.RData ../reportData/impute5.all.RData ../reportData/impute5.lr550.RData ../reportData/impute5.dly550.RData ../reportData/impute5.nodly550.RData ../reportData/pharp.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figure_imputation.R ../report/figure_imputation.pdf $(FONT) ../reportData/all.maf.hist.RData ../reportData/chip.maf.hist.RData ../reportData/impute5.all.RData ../reportData/impute5.lr550.RData ../reportData/impute5.dly550.RData ../reportData/impute5.nodly550.RData ../reportData/pharp.RData
	
../report/figure_imputation_software.pdf: figure_imputation_software.R ../reportData/impute5.lr550.RData ../reportData/beagle.lr550.RData ../reportData/minimac4.lr550.RData ../reportData/impute5.dly550.RData ../reportData/beagle.dly550.RData ../reportData/minimac4.dly550.RData ../reportData/impute5.nodly550.RData ../reportData/beagle.nodly550.RData ../reportData/minimac4.nodly550.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figure_imputation_software.R ../report/figure_imputation_software.pdf $(FONT) ../reportData/impute5.lr550.RData ../reportData/beagle.lr550.RData ../reportData/minimac4.lr550.RData ../reportData/impute5.dly550.RData ../reportData/beagle.dly550.RData ../reportData/minimac4.dly550.RData ../reportData/impute5.nodly550.RData ../reportData/beagle.nodly550.RData ../reportData/minimac4.nodly550.RData
	
../report/figure_bf.pdf: figure_bf.R ../reportData/bf.s21.pval.RData ../reportData/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.fai ../reportData/bf.man.png ../reportData/bf.s21.impute.pval.r2.RData ../reportData/chr1.mc4r.gene.bound
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figure_bf.R ../report/figure_bf.pdf $(FONT) ../reportData/bf.s21.pval.RData ../reportData/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.fai ../reportData/bf.man.png ../reportData/bf.s21.impute.pval.r2.RData ../reportData/chr1.mc4r.gene.bound
		
../report/figure_bl.pdf: figure_bl.R ../reportData/bl.w64.pval.RData ../reportData/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.fai ../reportData/bl.man.png ../reportData/bl.w64.impute.pval.r2.RData ../reportData/BL_topSNP_phe.csv
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figure_bl.R ../report/figure_bl.pdf $(FONT) ../reportData/bl.w64.pval.RData ../reportData/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.fai ../reportData/bl.man.png ../reportData/bl.w64.impute.pval.r2.RData ../reportData/chr17.bmp2.gene.bound ../reportData/BL_topSNP_phe.csv
			
../report/figure_bw.pdf: figure_bw.R ../reportData/bw.s22.pval.RData ../reportData/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.fai ../reportData/bw.man.png ../reportData/bw.s22.impute.pval.r2.RData ../reportData/chr2.slc.gene.bound
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figure_bw.R ../report/figure_bw.pdf $(FONT) ../reportData/bw.s22.pval.RData ../reportData/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.fai ../reportData/bw.man.png ../reportData/bw.s22.impute.pval.r2.RData ../reportData/chr2.slc.gene.bound
		