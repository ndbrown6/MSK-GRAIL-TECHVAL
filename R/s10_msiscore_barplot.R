#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
source('config.R')

if (!dir.exists("../res/figureS10")) {
	dir.create("../res/figureS10")
}

all_msi_df = read_tsv(file=msi_processed_data_file)

alias_df = read_tsv(file=subject_alias_file)

type_frame_df = data_frame(type=c("VB","VL","VP"), expand=c("Metastatic Breast","Metastatic Lung","Metastatic Prostate"), subject_type=c("Breast","Lung","Prostate"))

useful_msi_df = all_msi_df %>%
				left_join(alias_df %>% mutate(id=patient_id)) %>%
				mutate(id=alias) %>%
				select(patient_id,original_msi,fixed_msi,subject_type,sample)
				
pdf(file="../res/figureS10/all_msi_score_original.pdf", height=10, width=18)
par(mar = c(6.1, 6, 4.1, 1), mfcol=c(1,3))
subject_types = c("Breast", "Lung", "Prostate")
for (i in 1:length(subject_types)) {
	msi_df_tumor = useful_msi_df %>%
			 	   filter(subject_type == subject_types[i]) %>%
			 	   filter(sample == "tumor") %>%
			 	   .[["original_msi"]]
	msi_df_cfdna = useful_msi_df %>%
			 	   filter(subject_type == subject_types[i]) %>%
			 	   filter(sample == "cfdna") %>%
			 	   .[["original_msi"]]
	msi_df = cbind(msi_df_cfdna, msi_df_tumor)
	barplot(t(msi_df)+.001, beside=TRUE, horiz=TRUE, col=c("forestgreen","#612B5C"), border=NA, axes=FALSE, xlim=c(0,.31), ylim=c(1,nrow(msi_df)*2), space=rep(.1, nrow(msi_df)*2))
	axis(1, at = NULL, cex.axis = 1.75, padj = 0.25, lwd=1.55, lwd.ticks=1.55)
	mtext(side = 1, text = "MSI score", line = 4, cex = 1.35)
	legend("bottomright", pch=15, col=c("forestgreen","#612B5C"), legend=c(" cfDNA"," Tumor"), box.lwd=-1, cex=1.75, pt.cex=2.5)
}
dev.off()

pdf(file="../res/figureS10/all_msi_score_fixed.pdf", height=10, width=18)
par(mar = c(6.1, 6, 4.1, 1), mfcol=c(1,3))
subject_types = c("Breast", "Lung", "Prostate")
for (i in 1:length(subject_types)) {
	msi_df_tumor = useful_msi_df %>%
			 	   filter(subject_type == subject_types[i]) %>%
			 	   filter(sample == "tumor") %>%
			 	   .[["fixed_msi"]]
	msi_df_cfdna = useful_msi_df %>%
			 	   filter(subject_type == subject_types[i]) %>%
			 	   filter(sample == "cfdna") %>%
			 	   .[["fixed_msi"]]
	msi_df = cbind(msi_df_cfdna, msi_df_tumor)
	barplot(t(msi_df)+.001, beside=TRUE, horiz=TRUE, col=c("forestgreen","#612B5C"), border=NA, axes=FALSE, xlim=c(0,.31), ylim=c(1,nrow(msi_df)*2), space=rep(.1, nrow(msi_df)*2))
	axis(1, at = NULL, cex.axis = 1.75, padj = 0.25, lwd=1.55, lwd.ticks=1.55)
	mtext(side = 1, text = "MSI score", line = 4, cex = 1.35)
	legend(x="bottomright", pch=15, col=c("forestgreen","#612B5C"), legend=c(" cfDNA"," Tumor"), box.lwd=-1, cex=1.75, pt.cex=2.5)
}
dev.off()
