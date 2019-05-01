#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/rebuttal")) {
	dir.create("../res/rebuttal")
}

#==================================================
# Heatmap of CH related variants in WBC
#==================================================
maf = read.csv(file="../res/tables/Table_S8.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
	  arrange(Variant_Classification, Hugo_Symbol, Tumor_Sample_Barcode, Chromosome, Start_Position)
write.table(maf, file="../res/tables/Table_S8_ch_sorted_maftools.maf", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
wbc = read.maf(maf="../res/tables/Table_S8_ch_sorted_maftools.maf", clinicalData="../res/tables/clinical.tsv")

pdf(file="../res/rebuttal/OncoPlot_WBC.pdf", width=12)
oncoplot(maf = wbc, genes = chip_genes,
		 removeNonMutated = FALSE,
		 sortByMutation = TRUE,
		 clinicalFeatures = "Tissue",
		 sortByAnnotation = TRUE,
		 annotationColor = list("Tissue" = c(Breast = "salmon",
		 						Lung = "#FDAE61",
		 						Prostate = "#ABDDA4",
		 						Healthy = "cadetblue")),
		 fontSize = 10)
dev.off()
