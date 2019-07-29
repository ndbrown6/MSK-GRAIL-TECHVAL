#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/figureS9")) {
	dir.create("../res/figureS9")
}

#==================================================
# heatmap of ch-related variants in wbc
#==================================================
if (!file.exists("../res/tables/Table_S8_ch_sorted_maftools.maf")) {
	maf = read_csv(file="../res/tables/Table_S8.tsv", col_types = cols(.default = col_character())) %>%
		  type_convert() %>%
		  dplyr::arrange(Variant_Classification,
	  				 	 Hugo_Symbol,
	  				 	 Tumor_Sample_Barcode,
	  				 	 Chromosome,
	  				 	 Start_Position)

	write_tsv(x=maf, path="../res/tables/Table_S8_ch_sorted_maftools.maf", append=FALSE, col_names=TRUE)
}
variants = read.maf(maf="../res/tables/Table_S8_ch_sorted_maftools.maf", clinicalData="../res/tables/clinical.tsv")

pdf(file="../res/figureS9/onco_plot_wbc.pdf", width=12)
oncoplot(maf = variants, genes = chip_genes,
		 drawRowBar = TRUE, drawColBar = TRUE,
		 clinicalFeatures = "Tissue",
		 annotationColor = list("Tissue" = c(Breast = as.character(cohort_cols["Breast"]),
		 									 Lung = as.character(cohort_cols["Lung"]),
		 									 Prostate = as.character(cohort_cols["Prostate"]),
		 									 Healthy = as.character(cohort_cols["Control"]))),
		 removeNonMutated = FALSE,
		 sortByMutation = TRUE,
		 sortByAnnotation = TRUE,		 
		 fontSize = 10)
dev.off()

#==================================================
# oncodrive of ch-related variants in wbc
#==================================================
variants_oncodrive = m0_oncodrive(maf = variants, ignoreGenes=c("MST1", "JAK2", "STAT5B", "GNAS")) %>%
					 type_convert() %>%
					 dplyr::select(symbol = Hugo_Symbol,
					 			   frac = fract_muts_in_clusters,
					 			   pval = pval,
					 			   fdr = fdr,
					 			   total = total) %>%
					 mutate(signif = ifelse(pval<.05, "+", "-"))

plot.0 = ggplot(variants_oncodrive, aes(x=frac, y=-log10(fdr), size=total, fill=signif, color=signif, label=symbol)) +
		 geom_point(alpha=1, shape=21) +
		 scale_fill_manual(values = c("+" = "#BE1E2D", "-" = "blue")) +
		 scale_color_manual(values = c("+" = "#BE1E2D", "-" = "blue")) +
		 theme_classic(base_size=16) +
		 labs(x = "\nFraction of mutations in cluster", y = expression(-Log[10]~"(FDR)")) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=13)) +
		 theme(legend.justification = c(1, 0),
		 	   legend.position = c(1, 0),
		 	   legend.title = element_blank(),
		 	   legend.background = element_blank(),
		 	   legend.text=element_text(size=8)) +
		 guides(fill=guide_legend(title="Significance")) +
		 guides(size=guide_legend(title="Number of\nmutations")) +
		 guides(color="none") +
		 geom_text_repel() +
		 ylim(0,1.25)

pdf(file="../res/figureS9/onco_drive_wbc.pdf", width=7, height=7)
print(plot.0)
dev.off()

#==================================================
# lollipop of ppm1d and dnmt3a
#==================================================
pdf(file="../res/figureS9/lol_ppm1d.pdf", height=4, width=12)
lollipopPlot(maf=variants, gene="PPM1D", pointSize=2.5)
dev.off()

pdf(file="../res/figureS9/lol_dnmt3a.pdf", height=4, width=12)
lollipopPlot(maf=variants, gene="DNMT3A", pointSize=2.5)
dev.off()
