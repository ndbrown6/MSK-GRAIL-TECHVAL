#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/rebuttal")) {
	dir.create("../res/rebuttal")
}

tracker_grail = read_csv(file=patient_tracker, col_types = cols(.default = col_character()))  %>%
 				type_convert()
tracker_impact = read_csv(impact_tracker, col_types = cols(.default = col_character()))  %>%
				 type_convert()
valid_patient_ids = tracker_grail %>%
 	  				filter(patient_id %in% tracker_impact$patient_id | tissue == "Healthy") %>%
 	  				filter(!(tissue %in% c("Breast", "Lung", "Prostate") & study=="Merlin")) %>%
 	  				.[["patient_id"]]
clinical = read_tsv(clinical_file, col_types = cols(.default = col_character())) %>%
 		   type_convert() %>%
 		   rename(localisation = tissue)
 			   
valid_patient_ids = intersect(valid_patient_ids, clinical$patient_id)

qc_metrics_cfdna = read.csv(file="../modified_v11/QC_metrics/TechVal_Merlin_QC_metrics.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
			 	   select(sample_id, patient_id, sample_type, tissue, volume_of_blood_mL, volume_of_DNA_source_mL, DNA_extraction_yield_ng, DNA_input_concentration_ng_uL, Library_preparation_input_ng, raw.MEAN_BAIT_COVERAGE, collapsed.MEAN_BAIT_COVERAGE, collapsed_fragment.MEAN_BAIT_COVERAGE, readErrorRate, readSubstErrorRate, Study) %>%
			 	   filter(sample_type=="cfDNA")
tracker_grail_cfdna = read.csv(file=patient_tracker, header=TRUE, sep=",", stringsAsFactors=FALSE) %>%
					  select(patient_id, cfdna_sample_id) %>%
					  rename(msk_id = patient_id, sample_id = cfdna_sample_id)
qc_metrics_cfdna = left_join(qc_metrics_cfdna, tracker_grail_cfdna, by="sample_id")

qc_metrics_wbc = read.csv(file="../modified_v11/QC_metrics/TechVal_Merlin_QC_metrics.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
			 	   select(sample_id, patient_id, sample_type, tissue, volume_of_blood_mL, volume_of_DNA_source_mL, DNA_extraction_yield_ng, DNA_input_concentration_ng_uL, Library_preparation_input_ng, raw.MEAN_BAIT_COVERAGE, collapsed.MEAN_BAIT_COVERAGE, collapsed_fragment.MEAN_BAIT_COVERAGE, readErrorRate, readSubstErrorRate, Study) %>%
			 	   filter(sample_type=="gDNA")
tracker_grail_wbc = read.csv(file=patient_tracker, header=TRUE, sep=",", stringsAsFactors=FALSE) %>%
					select(patient_id, gdna_sample_id) %>%
					rename(msk_id = patient_id, sample_id = gdna_sample_id)
qc_metrics_wbc = left_join(qc_metrics_wbc, tracker_grail_wbc, by="sample_id")

qc_metrics = rbind(qc_metrics_cfdna, qc_metrics_wbc) %>%
			 filter(msk_id %in% valid_patient_ids) %>%
			 select(-sample_id, -msk_id) %>%
			 rename(Patient_ID = patient_id,
			 		Sample_Type = sample_type,
			 		Tissue = tissue,
			 		Volume_of_blood_mL = volume_of_blood_mL,
			 		Volume_of_DNA_source_mL = volume_of_DNA_source_mL,
			 		Uncollapsed_Mean_Coverage = raw.MEAN_BAIT_COVERAGE,
				   	Collapsed_Mean_Coverage = collapsed.MEAN_BAIT_COVERAGE,
				   	Collapsed_Fragment_Mean_Coverage = collapsed_fragment.MEAN_BAIT_COVERAGE,
				   	Indel_and_Substitution_Error_Rate = readErrorRate,
				   	Substitution_Error_Rate = readSubstErrorRate,
				   	Assay_Version = Study) %>%
			mutate(Assay_Version = ifelse(Assay_Version=="TechVal", "V1", "V2")) %>%
			arrange(Patient_ID, Sample_Type) %>%
			mutate(Library_preparation_input_ng = ifelse(Library_preparation_input_ng>75, 75, Library_preparation_input_ng))

#==================================================
# Comparison of sequencing depth of cases and
# control cfDNA samples
#==================================================
tmp = qc_metrics %>%
	  filter(Sample_Type=="cfDNA") %>%
	  mutate(Tissue = ifelse(Tissue=="Healthy", "Control", "Cancer"))
plot.0 = ggplot(tmp, aes(x = Tissue, y = Uncollapsed_Mean_Coverage, fill = Tissue)) + 
		 geom_boxplot(alpha=1, outlier.size=2.5, outlier.shape=21) + 
		 scale_fill_manual(values = c("salmon", "cadetblue")) + 
		 facet_wrap(~Sample_Type) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="", y="Uncollapsed mean coverage\n") +
		 guides(fill=FALSE)
pdf(file="../res/rebuttal/UNCOLLAPSED_MEAN_BAIT_COVERAGE.pdf", width=4, height=6)
print(plot.0)
dev.off()

tmp = qc_metrics %>%
	  filter(Sample_Type=="cfDNA") %>%
	  mutate(Tissue = ifelse(Tissue=="Healthy", "Control", "Cancer"))
plot.0 = ggplot(tmp, aes(x = Tissue, y = Collapsed_Mean_Coverage, fill = Tissue)) + 
		 geom_boxplot(alpha=1, outlier.size=2.5, outlier.shape=21) + 
		 scale_fill_manual(values = c("salmon", "cadetblue")) + 
		 facet_wrap(~Sample_Type) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="", y="Collapsed mean coverage\n") +
		 guides(fill=FALSE)
pdf(file="../res/rebuttal/COLLAPSED_MEAN_BAIT_COVERAGE.pdf", width=4, height=6)
print(plot.0)
dev.off()

tmp = qc_metrics %>%
	  filter(Sample_Type=="cfDNA") %>%
	  mutate(Tissue = ifelse(Tissue=="Healthy", "Control", "Cancer"))
plot.0 = ggplot(tmp, aes(x = Tissue, y = Collapsed_Fragment_Mean_Coverage, fill = Tissue)) + 
		 geom_boxplot(alpha=1, outlier.size=2.5, outlier.shape=21) + 
		 scale_fill_manual(values = c("salmon", "cadetblue")) + 
		 facet_wrap(~Sample_Type) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="", y="Unique fragments mean coverage\n") +
		 guides(fill=FALSE)
pdf(file="../res/rebuttal/COLLAPSED_FRAGMENT_MEAN_BAIT_COVERAGE.pdf", width=4, height=6)
print(plot.0)
dev.off()
 
#==================================================
# Comparison of sequencing depth cfDNA and gDNA by
# tissue type
#==================================================
tmp = qc_metrics %>%
	  mutate(Tissue = ifelse(Tissue=="Healthy", "Control", Tissue)) %>%
	  mutate(Tissue = factor(Tissue, levels=c("Breast","Lung","Prostate", "Control")))
plot.0 = ggplot(tmp, aes(x = Tissue, y = Uncollapsed_Mean_Coverage, fill = Tissue)) + 
		 geom_boxplot(alpha=1, outlier.size=2.5, outlier.shape=21) + 
		 scale_fill_manual(values = c("salmon", "#FDAE61", "#ABDDA4", "cadetblue")) + 
		 facet_wrap(~Sample_Type) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="", y="Uncollapsed mean coverage\n") +
		 guides(fill=FALSE) +
		 scale_y_continuous(breaks=c(50000, 75000, 100000, 125000), labels=c("50000", "75000", "100000", "125000"))
pdf(file="../res/rebuttal/UNCOLLAPSED_MEAN_BAIT_COVERAGE_cfDNA_vs_gDNA_by_tissue.pdf", width=8, height=6)
print(plot.0)
dev.off()

tmp = qc_metrics %>%
	  mutate(Tissue = ifelse(Tissue=="Healthy", "Control", Tissue)) %>%
	  mutate(Tissue = factor(Tissue, levels=c("Breast","Lung","Prostate", "Control")))
plot.0 = ggplot(tmp, aes(x = Tissue, y = Collapsed_Mean_Coverage, fill = Tissue)) + 
		 geom_boxplot(alpha=1, outlier.size=2.5, outlier.shape=21) + 
		 scale_fill_manual(values = c("salmon", "#FDAE61", "#ABDDA4", "cadetblue")) + 
		 facet_wrap(~Sample_Type) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="", y="Collapsed mean coverage\n") +
		 guides(fill=FALSE)
pdf(file="../res/rebuttal/COLLAPSED_MEAN_BAIT_COVERAGE_cfDNA_vs_gDNA_by_tissue.pdf", width=8, height=6)
print(plot.0)
dev.off()

#==================================================
# Correlation of input DNA to library preparation
# with mean target coverage
#==================================================
tmp = qc_metrics %>%
	  filter(Sample_Type=="cfDNA") %>%
	  mutate(Tissue = ifelse(Tissue=="Healthy", "Control", Tissue)) %>%
	  mutate(Tissue = factor(Tissue, levels=c("Breast","Lung","Prostate", "Control")))
fit = lm(Collapsed_Mean_Coverage ~ Library_preparation_input_ng, data=tmp)
plot.0 = ggplot(tmp, aes(x = Library_preparation_input_ng, y = Collapsed_Mean_Coverage, fill = Tissue)) + 
		 geom_point(alpha=1, shape=21, size=3.5) +
		 scale_fill_manual(values = c("salmon", "#FDAE61", "#ABDDA4", "cadetblue")) + 
		 facet_wrap(~Sample_Type) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="\nAmount of DNA used\nfor library preparation (ng)\n", y="Collapsed mean coverage\n") +
		 geom_abline(slope = fit$coefficients[2], intercept = fit$coefficients[1], col="orange3") +
		 theme(legend.justification = c(1, 0),
		 	   legend.position = c(1, 0),
		 	   legend.title = element_blank(),
		 	   legend.background = element_blank(),
		 	   legend.text=element_text(size=8)) +
		 coord_cartesian(xlim = c(10, 80))

pdf(file="../res/rebuttal/COLLAPSED_MEAN_BAIT_COVERAGE_cfDNA_vs_input_DNA.pdf", width=6, height=6)
print(plot.0)
dev.off()

#==================================================
# Correlation of input DNA to library preparation
# with sample type
#==================================================
tmp = qc_metrics %>%
	  filter(Sample_Type=="cfDNA") %>%
	  mutate(Tissue = ifelse(Tissue=="Healthy", "Control", Tissue)) %>%
	  mutate(Tissue = factor(Tissue, levels=c("Breast","Lung","Prostate", "Control")))  
plot.0 = ggplot(tmp, aes(x = Tissue, y = Library_preparation_input_ng, fill = Tissue)) + 
		 geom_boxplot(alpha=1, outlier.shape=21, outlier.size=3.5) +
		 scale_fill_manual(values = c("salmon", "#FDAE61", "#ABDDA4", "cadetblue")) +
		 facet_wrap(~Sample_Type) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="", y="Amount of DNA used\nfor library preparation (ng)\n") +
		 coord_cartesian(ylim = c(10, 80)) +
		 guides(fill=FALSE)
pdf(file="../res/rebuttal/Input_cfDNA_by_Tissue.pdf", width=6, height=6)
print(plot.0)
dev.off()

#==================================================
# Comparison of error rate
#==================================================
tmp = qc_metrics %>%
	  mutate(Tissue = ifelse(Tissue=="Healthy", "Control", Tissue)) %>%
	  mutate(Tissue = factor(Tissue, levels=c("Breast","Lung","Prostate", "Control")))
plot.0 = ggplot(tmp, aes(x = Tissue, y = 1e5*Indel_and_Substitution_Error_Rate, fill = Tissue)) + 
		 geom_boxplot(alpha=1, outlier.size=2.5, outlier.shape=21) + 
		 scale_fill_manual(values = c("salmon", "#FDAE61", "#ABDDA4", "cadetblue")) + 
		 facet_wrap(~Sample_Type) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="", y="% collapsed bases (1E-05)\n") +
		 guides(fill=FALSE) +
		 coord_cartesian(ylim = c(5,11))
pdf(file="../res/rebuttal/COMBINED_ERROR_RATE_cfDNA_vs_gDNA_by_tissue.pdf", width=8, height=6)
print(plot.0)
dev.off()

tmp = qc_metrics %>%
	  mutate(Tissue = ifelse(Tissue=="Healthy", "Control", Tissue)) %>%
	  mutate(Tissue = factor(Tissue, levels=c("Breast","Lung","Prostate", "Control")))
plot.0 = ggplot(tmp, aes(x = Tissue, y = 1e5*Substitution_Error_Rate, fill = Tissue)) + 
		 geom_boxplot(alpha=1, outlier.size=2.5, outlier.shape=21) + 
		 scale_fill_manual(values = c("salmon", "#FDAE61", "#ABDDA4", "cadetblue")) + 
		 facet_wrap(~Sample_Type) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="", y="% collapsed bases (1E-05)\n") +
		 guides(fill=FALSE) +
		 coord_cartesian(ylim = c(1,7))
pdf(file="../res/rebuttal/SUBSTITUTION_ERROR_RATE_cfDNA_vs_gDNA_by_tissue.pdf", width=8, height=6)
print(plot.0)
dev.off()
