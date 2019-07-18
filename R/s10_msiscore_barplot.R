#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
source('config.R')

if (!dir.exists("../res/figureS10")) {
	dir.create("../res/figureS10")
}

all_msi_df = read_tsv(file=url_msi_processed_data)

alias_df = read_tsv(file=url_subject_alias)

type_frame_df = data_frame(type=c("VB","VL","VP"), expand=c("Metastatic Breast","Metastatic Lung","Metastatic Prostate"), subject_type=c("Breast","Lung","Prostate"))

useful_msi_df = all_msi_df %>%
				left_join(alias_df %>% mutate(id=patient_id)) %>%
				mutate(id=alias) %>%
				dplyr::select(patient_id,original_msi,fixed_msi,subject_type,sample)
				
cols = c("cfDNA"="#d7191c",
		 "Tumor"="#fdae61")
		 
msi_tumor = useful_msi_df %>%
		    filter(sample=="tumor") %>%
		    arrange(patient_id)
msi_cfdna = useful_msi_df %>%
		    filter(sample=="cfdna") %>%
		    arrange(patient_id)

index = order(msi_tumor$original_msi)
msi_tumor = msi_tumor[index,,drop=FALSE]
msi_cfdna = msi_cfdna[index,,drop=FALSE]
msi_tumor = bind_cols(msi_tumor, "level" = factor(msi_tumor$patient_id, levels=msi_tumor$patient_id, ordered=TRUE))
msi_cfdna = bind_cols(msi_cfdna, "level" = factor(msi_cfdna$patient_id, levels=msi_cfdna$patient_id, ordered=TRUE))
msi_df = bind_rows(msi_tumor, msi_cfdna)
		 
plot.0 = msi_df %>%
		 mutate(sample = ifelse(sample=="cfdna", "cfDNA", "Tumor")) %>%
		 ggplot(aes(x=level, y=original_msi+.001, group=sample, fill=sample)) +
		 geom_bar(stat="identity", position=position_dodge()) +
		 facet_wrap(~subject_type, scales = "free_y") +
		 theme_bw(base_size=15) +
		 coord_flip() +
		 labs(y="\nMSI score\n", x="\n") +
		 scale_fill_manual(values = cols) +
		 guides(fill=guide_legend(title=c("Sample type"))) +
		 theme(legend.title=element_text(size=14))


pdf(file="../res/figureS10/all_msi_score_original.pdf", height=10, width=18)
print(plot.0)
dev.off()

plot.0 = msi_df %>%
		 mutate(sample = ifelse(sample=="cfdna", "cfDNA", "Tumor")) %>%
		 ggplot(aes(x=level, y=fixed_msi+.001, group=sample, fill=sample)) +
		 geom_bar(stat="identity", position=position_dodge()) +
		 facet_wrap(~subject_type, scales = "free_y") +
		 theme_bw(base_size=15) +
		 coord_flip() +
		 labs(y="\nMSI score\n", x="\n") +
		 scale_fill_manual(values = cols) +
		 guides(fill=guide_legend(title=c("Sample type"))) +
		 theme(legend.title=element_text(size=14))


pdf(file="../res/figureS10/all_msi_score_fixed.pdf", height=10, width=18)
print(plot.0)
dev.off()
