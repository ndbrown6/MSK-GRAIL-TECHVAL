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
# Age by cancer status and treatment history
#==================================================
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

index_control = grepl("W", valid_patient_ids, fixed=TRUE)
index_cancer = grepl("V", valid_patient_ids, fixed=TRUE)

tmp = clinical %>%
	  select(patient_id, age) %>%
	  mutate(cat = "") %>%
	  mutate(cat = ifelse(patient_id %in% valid_patient_ids[index_control], "Control", cat)) %>%
	  mutate(cat = ifelse(patient_id %in% valid_patient_ids[index_cancer], "Cancer", cat)) %>%
	  filter(cat != "") %>%
	  mutate(facet = "Age by cancer status")
	  
	  
p = signif(wilcox.test(tmp$age[tmp$cat=="Cancer"], tmp$age[tmp$cat=="Control"], alternative="two.sided", correct=FALSE)$p.value, 3)
		
plot.0 = ggplot(tmp, aes(x=cat, y=age, fill = cat)) +
		 geom_boxplot(stat = "boxplot") +
		 scale_fill_manual(values = rev(c("#FDAE61", "salmon"))) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="", y="Age (years)\n") +
 		 coord_cartesian(ylim = c(20,100)) +
		 guides(fill=FALSE) +
		 facet_wrap(~facet) +
		 theme_bw(base_size=15) +
  		 annotate(geom="text", x=1.5, y=100, label=paste0("P = ", p))
	  	
pdf(file="../res/rebuttal/Age_by_Cancer.pdf", width=5, height=5)
print(plot.0)
dev.off()

tx = read.csv(file="../res/etc/prior_tx_techval_0818.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
	 rename(patient_id = ID)
clinical = 	left_join(clinical, tx, by="patient_id") %>%
			mutate(prior_rt = ifelse(subj_type=="Healthy", 0, prior_rt)) %>%
			mutate(prior_chemo = ifelse(subj_type=="Healthy", 0, prior_chemo)) %>%
			filter(patient_id %in% valid_patient_ids)

tmp = clinical %>%
	  select(patient_id, age) %>%
	  mutate(cat = "No RT/CT") %>%
	  mutate(cat = ifelse(clinical$prior_rt==1 | clinical$prior_chemo==1, "RT/CT", cat)) %>%
	  mutate(cat = factor(cat, levels=c("RT/CT","No RT/CT"))) %>%
	  mutate(facet = "Age by treatment history")
	  
p = signif(wilcox.test(tmp$age[tmp$cat=="No RT/CT"], tmp$age[tmp$cat=="RT/CT"], alternative="two.sided", correct=FALSE)$p.value, 3)
		
plot.0 = ggplot(tmp, aes(x=cat, y=age, fill = cat)) +
		 geom_boxplot(stat = "boxplot") +
		 scale_fill_manual(values = rev(c("#FDAE61", "salmon"))) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="", y="Age (years)\n") +
 		 coord_cartesian(ylim = c(20,100)) +
		 guides(fill=FALSE) +
		 facet_wrap(~facet) +
		 theme_bw(base_size=15) +
  		 annotate(geom="text", x=1.5, y=100, label=paste0("P = ", p))
	  	
pdf(file="../res/rebuttal/Age_by_Treatment.pdf", width=5, height=5)
print(plot.0)
dev.off()

#==================================================
# Comparison of sequencing depth of cases and
# control cfDNA samples
#==================================================
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

tmp.0 = qc_metrics %>%
	  filter(Sample_Type=="cfDNA") %>%
	  arrange(Patient_ID)
tmp.1 = qc_metrics %>%
	  filter(Sample_Type=="gDNA") %>%
	  arrange(Patient_ID)

tmp = tmp.0 %>%
	  select(Tissue, Uncollapsed_Mean_Coverage) %>%
	  mutate(Ratio = Uncollapsed_Mean_Coverage/tmp.1$Uncollapsed_Mean_Coverage) %>%
	  mutate(Tissue = ifelse(Tissue == "Healthy", "Control", "Cancer"))

p = signif(wilcox.test(tmp$Ratio[tmp$Tissue=="Cancer"], tmp$Ratio[tmp$Tissue=="Control"], alternative="two.sided", correct=FALSE)$p.value, 3)
		
plot.0 = ggplot(tmp, aes(x=Tissue, y=Ratio)) +
		 geom_boxplot(outlier.colour=NA, outlier.fill=NA, outlier.shape=NA, fill="steelblue") +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="", y="Ratio uncollapsed coverage (cfDNA/gDNA)\n") +
		 guides(fill=FALSE) +
  		 annotate(geom="text", x=1.5, y=3, label=paste0("P = ", p)) +
  		 coord_cartesian(ylim = c(0,3))
		 
pdf(file="../res/rebuttal/Ratio_Uncollapsed_Coverage_by_Cancer.pdf", width=4, height=6)
print(plot.0)
dev.off()

tmp = tmp.0 %>%
	  select(Tissue, Collapsed_Mean_Coverage) %>%
	  mutate(Ratio = Collapsed_Mean_Coverage/tmp.1$Collapsed_Mean_Coverage) %>%
	  mutate(Tissue = ifelse(Tissue == "Healthy", "Control", "Cancer"))
	  
p = signif(wilcox.test(tmp$Ratio[tmp$Tissue=="Cancer"], tmp$Ratio[tmp$Tissue=="Control"], alternative="two.sided", correct=FALSE)$p.value, 3)

plot.0 = ggplot(tmp, aes(x=Tissue, y=Ratio)) +
		 geom_boxplot(outlier.colour=NA, outlier.fill=NA, outlier.shape=NA, fill="steelblue") +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="", y="Ratio collapsed coverage (cfDNA/gDNA)\n") +
		 guides(fill=FALSE) +
  		 annotate(geom="text", x=1.5, y=3, label=paste0("P = ", p)) +
   		 coord_cartesian(ylim = c(0,3))
		 
pdf(file="../res/rebuttal/Ratio_Collapsed_Coverage_by_Cancer.pdf", width=4, height=6)
print(plot.0)
dev.off()

tx = read.csv(file="../res/etc/prior_tx_techval_0818.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
	 rename(Patient_ID = ID)
tmp = tmp.0 %>%
	  left_join(tx) %>%
	  mutate(prior_rt = ifelse(Tissue=="Healthy", 0, prior_rt)) %>%
	  mutate(prior_chemo = ifelse(Tissue=="Healthy", 0, prior_chemo)) %>%
	  mutate(cat = "No RT/CT") %>%
	  mutate(cat = ifelse(prior_rt==1 | prior_chemo==1, "RT/CT", cat)) %>%
	  select(cat, Uncollapsed_Mean_Coverage) %>%
	  mutate(Ratio = Uncollapsed_Mean_Coverage/tmp.1$Uncollapsed_Mean_Coverage)
	  
p = signif(wilcox.test(tmp$Ratio[tmp$cat=="RT/CT"], tmp$Ratio[tmp$cat=="No RT/CT"], alternative="two.sided", correct=FALSE)$p.value, 3)
		
plot.0 = ggplot(tmp, aes(x=cat, y=Ratio)) +
		 geom_boxplot(outlier.colour=NA, outlier.fill=NA, outlier.shape=NA, fill="salmon") +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="", y="Ratio uncollapsed coverage (cfDNA/gDNA)\n") +
		 guides(fill=FALSE) +
  		 annotate(geom="text", x=1.5, y=3, label=paste0("P = ", p)) +
  		 coord_cartesian(ylim = c(0,3))
		 
pdf(file="../res/rebuttal/Ratio_Uncollapsed_Coverage_by_Treatment.pdf", width=4, height=6)
print(plot.0)
dev.off()

tmp = tmp.0 %>%
	  left_join(tx) %>%
	  mutate(prior_rt = ifelse(Tissue=="Healthy", 0, prior_rt)) %>%
	  mutate(prior_chemo = ifelse(Tissue=="Healthy", 0, prior_chemo)) %>%
	  mutate(cat = "No RT/CT") %>%
	  mutate(cat = ifelse(prior_rt==1 | prior_chemo==1, "RT/CT", cat)) %>%
	  select(cat, Collapsed_Mean_Coverage) %>%
	  mutate(Ratio = Collapsed_Mean_Coverage/tmp.1$Collapsed_Mean_Coverage)
	  
p = signif(wilcox.test(tmp$Ratio[tmp$cat=="RT/CT"], tmp$Ratio[tmp$cat=="No RT/CT"], alternative="two.sided", correct=FALSE)$p.value, 3)
		
plot.0 = ggplot(tmp, aes(x=cat, y=Ratio)) +
		 geom_boxplot(outlier.colour=NA, outlier.fill=NA, outlier.shape=NA, fill="salmon") +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="", y="Ratio collapsed coverage (cfDNA/gDNA)\n") +
		 guides(fill=FALSE) +
  		 annotate(geom="text", x=1.5, y=3, label=paste0("P = ", p)) +
  		 coord_cartesian(ylim = c(0,3))
		 
pdf(file="../res/rebuttal/Ratio_Collapsed_Coverage_by_Treatment.pdf", width=4, height=6)
print(plot.0)
dev.off()
	  
#==================================================
# Boxplot of number CH-derived variants by age in
# cfDNA
#==================================================
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
		   
snv_vars = read_tsv(snv_file$scored, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(level_2a = as.character(level_2a)) %>%
		   mutate(level_r1 = as.character(level_r1))
		   
indel_vars = read_tsv(indel_file$scored, col_types = cols(.default = col_character())) %>%
			 type_convert() %>%
		     mutate(level_2a = as.character(level_2a)) %>%
		     mutate(level_r1 = as.character(level_r1))

wbc_stack = read_tsv(wbc_variants$scored, col_types = cols(.default = col_character())) %>%
			type_convert()
			
msk_anno = read_tsv(msk_anno_joined, col_types = cols(.default = col_character())) %>%
  		   type_convert()
			
tracker_grail = read_csv(file=patient_tracker)

tracker_impact = read_csv(file=impact_tracker)

bed_file = rtracklayer::import.bed(con=common_bed)
bed_ranges = GenomicRanges::ranges(bed_file)
total_bed_Mb = round(sum(GenomicRanges::width(bed_ranges)) / 1e6, 1)

valid_patient_ids = tracker_grail %>%
					filter(patient_id %in% tracker_impact$patient_id) %>%
					.[["patient_id"]]
  
indel_vars = indel_vars %>%
			 mutate(filter = replace(filter,
             		patient_id == "MSK-VB-0001" &
             		gene == "GATA3" &
             		filter == "PASS",
             		"CSR_MATCH_ELIMINATED"),
         	 		ccd = replace(ccd,
             			   		  patient_id == "MSK-VB-0001" &
                           		  gene == "GATA3" &
                           		  filter == "CSR_MATCH_ELIMINATED",
                           		  0))

snv_plasma = snv_vars %>%
  			 filter(ccd == 1,
         			(c_panel == 1 | panel == 1),
         			study == "TechVal",
         			grail == 1,
         			patient_id %in% valid_patient_ids) %>%
			 mutate(vtype = "SNV")

indel_plasma = indel_vars %>%
			   filter(ccd == 1,
         			  (c_panel == 1 | panel == 1),
         			  study == "TechVal",
         			  grail == 1,
         			  patient_id %in% valid_patient_ids) %>%
  			   mutate(vtype = "INDEL",
         			  altenddistmedian = as.integer(altenddistmedian))
         			  
healthy_snv = snv_vars %>%
  			  filter((c_panel == 1 | panel == 1),
         			  subj_type == "Healthy",
         			  grail == 1) %>%
  			  mutate(vtype = "SNV")

healthy_indel = indel_vars %>%
  				filter((c_panel == 1 | panel == 1),
         			    subj_type == "Healthy",
         				grail == 1) %>%
  				mutate(vtype = "INDEL",
         			   altenddistmedian = as.integer(altenddistmedian))

small_vars_plasma = full_join(snv_plasma, indel_plasma) %>%
					full_join(healthy_snv) %>%
					full_join(healthy_indel)
small_vars_plasma = small_vars_plasma %>%
  					mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
  					
small_vars_plasma = small_vars_plasma %>%
					mutate(loc = str_c(chrom, ":", position_orig, "_", ref_orig, ">", alt_orig))  					

variants = label_bio_source(small_vars_plasma)

variants = left_join(variants, msk_anno %>% select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
variants = variants %>%
		   mutate(bio_source = case_when(
		   					   MSK == 1 & grail == 1 ~ "biopsy_matched",
		   					   MSK == 1 & grail == 0 ~ "biopsy_only",
		   					   category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
		   					   category %in% c("germline", "germlineish") ~ "germline",
		   					   category %in% c("blood", "bloodier") ~ "WBC_matched",
		   					   category == "somatic" & `IM-T.alt_count` > bam_read_count_cutoff ~ "IMPACT-BAM_matched",
		   					   category == "somatic" ~ "VUSo",
		   					   TRUE ~ "other"),
		   		  af = ifelse(is.na(af), 0, af),
		   		  af_nobaq = round(adnobaq / dpnobaq * 100, 2),
		   		  af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq))
		   		  
tmp.0 = variants %>%
		filter(bio_source == "WBC_matched") %>%
		count(patient_id)
		   		  
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

tmp.1 = clinical %>%
		filter(patient_id %in% valid_patient_ids) %>%
		select(patient_id, age)
		
tmp.1 = tmp.1 %>%
		left_join(tmp.0, by="patient_id") %>%
		mutate(n = ifelse(is.na(n), 0, n)) %>%
		mutate(age_cat = case_when(
			age >= 20 & age <= 39 ~ "20-39",
			age >= 40 & age <= 49 ~ "40-49",
			age >= 50 & age <= 59 ~ "50-59",
			age >= 60 & age <= 69 ~ "60-69",
			age >= 70 ~ ">70"
		)) %>%
		mutate(case_cat = ifelse(grepl("V", patient_id), "Cancer", "Control")) %>%
		mutate(facet_1 = "CH variants in cfDNA by cancer status") %>%
		mutate(facet_2 = "CH variants in cfDNA by treatment history") %>%
		mutate(age_cat = factor(age_cat, levels = c("20-39", "40-49", "50-59", "60-69", ">70"))) %>%
		mutate(case_cat = factor(case_cat, levels = c("Control", "Cancer"))) %>%
		left_join(tx, by="patient_id") %>%
		mutate(treatment_cat = "No RT/CT") %>%
	  	mutate(treatment_cat = ifelse(prior_rt==1 | prior_chemo==1, "RT/CT", treatment_cat)) %>%
	  	mutate(treatment_cat = ifelse(is.na(treatment_cat), "No RT/CT", treatment_cat)) %>%
	  	mutate(treatment_cat = factor(treatment_cat, levels = c("No RT/CT", "RT/CT")))

plot.0 = ggplot(tmp.1, aes(x = age_cat, y = n, fill = case_cat)) + 
		 geom_boxplot(alpha=1, outlier.size=NA, outlier.shape=NA) + 
		 scale_fill_manual(values = c("#1F6784", "#EAA411")) + 
		 facet_wrap(~facet_1) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.165, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="", y="Number of variants\n") +
		 guides(fill=guide_legend(title=c("Status")))
		 
pdf(file="../res/rebuttal/CH_mutations_cfDNA_Cancer_Control.pdf", width=7, height=5)
print(plot.0)
dev.off()

plot.0 = ggplot(tmp.1, aes(x = age_cat, y = n, fill = treatment_cat)) + 
		 geom_boxplot(alpha=1, outlier.size=NA, outlier.shape=NA) + 
		 scale_fill_manual(values = c("#1F6784", "#EAA411")) + 
		 facet_wrap(~facet_2) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="", y="Number of variants\n") +
		 guides(fill=guide_legend(title=c("Treatment history")))
		 
pdf(file="../res/rebuttal/CH_mutations_cfDNA_RT_CT.pdf", width=7, height=5)
print(plot.0)
dev.off()

#==================================================
# Boxplot of number CH-derived variants by age in
# WBC
#==================================================
load(all_vars_and_clinical)

#==================================================
# Save variants
#==================================================
save_vars = all_vars
 
#==================================================
# Default filters
#==================================================
all_vars = all_vars %>%
 		   filter(is_patient_valid) %>%
 		   filter(c_panel) %>%
 		   filter(!is_hypermutator) %>%
 		   filter(!is_lowdepth) %>%
 		   filter(!is_lowqual) %>%
 		   filter(!is_tumor_matched) %>%
 		   filter(!is_cfdna_matched)
 
#==================================================
# < 5% recurrence | is_hotspot | frame-shifting
# in CH related gene
#==================================================
n_samples = all_vars %>%
 			distinct(patient_id) %>%
  		    count()
recurrence = all_vars %>%
   			 group_by(loc_lng) %>%
   			 count() %>%
   			 ungroup() %>%
   			 rename(n_recurrence=n) %>%
   			 mutate(f_recurrence=n_recurrence/n_samples$n)
all_vars = all_vars %>%
   		   left_join(recurrence)
all_vars = all_vars %>%
 		   filter(f_recurrence < 0.05 | is_hotspot | (Variant_Classification %in% c("Frame_Shift_Ins", "Frame_Shift_Del", "Nonsense_Mutation", "Nonstop_Mutation") & (SYMBOL %in% chip_genes)))
  		   
#==================================================
# Co-occurring indels filter
#==================================================
recurrence = all_vars %>%
   			 group_by(patient_id, loc_srt) %>%
   			 count() %>%
   			 ungroup() %>%
   			 rename(n_indel=n)
all_vars = all_vars %>%
		   left_join(recurrence)
all_vars = all_vars %>%
		   filter(!(n_indel > 1 & indel))
		   
#==================================================
# Variant classe filter
#==================================================
all_vars = all_vars %>%
 		   filter(Variant_Classification!="3'Flank") %>%
  		   filter(Variant_Classification!="3'UTR") %>%
  		   filter(Variant_Classification!="5'Flank") %>%
  		   filter(Variant_Classification!="5'UTR") %>%
  		   filter(Variant_Classification!="In_Frame_Del") %>%
  		   filter(Variant_Classification!="In_Frame_Ins") %>%
  		   filter(Variant_Classification!="Intron") %>%
  		   filter(Variant_Classification!="RNA") %>%
  		   filter(Variant_Classification!="Silent") %>%
  		   filter(Variant_Classification!="IGR") %>%
  		   filter(Variant_Classification!="Translation_Start_Site")
  		   
#==================================================
# HLA-A
#==================================================
all_vars = all_vars %>%
 		   filter(SYMBOL!="HLA-A")

#==================================================
# Germline filter
#==================================================
all_vars = all_vars %>%
 		   filter((adnobaq/dpnobaq)<=.3 | (Variant_Classification %in% c("Frame_Shift_Ins", "Frame_Shift_Del", "Nonsense_Mutation", "Nonstop_Mutation") & SYMBOL %in% chip_genes))
		   
#==================================================
# ExAC filter
#==================================================
all_vars = all_vars %>%
 		   filter(!in_exac)
 		   
#==================================================
# gnomAD filter
#==================================================
all_vars = all_vars %>%
 		   filter(!in_gnomad)
 		   

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
 			   
tmp.0 = all_vars %>%
		count(patient_id)
 			   
tmp.1 = clinical %>%
		filter(patient_id %in% valid_patient_ids) %>%
		select(patient_id, age)
		
tmp.1 = tmp.1 %>%
		left_join(tmp.0, by="patient_id") %>%
		mutate(n = ifelse(is.na(n), 0, n)) %>%
		mutate(age_cat = case_when(
			age >= 20 & age <= 39 ~ "20-39",
			age >= 40 & age <= 49 ~ "40-49",
			age >= 50 & age <= 59 ~ "50-59",
			age >= 60 & age <= 69 ~ "60-69",
			age >= 70 ~ ">70"
		)) %>%
		mutate(case_cat = ifelse(grepl("V", patient_id), "Cancer", "Control")) %>%
		mutate(facet_1 = "CH variants in WBC by cancer status") %>%
		mutate(facet_2 = "CH variants in WBC by treatment history") %>%
		mutate(age_cat = factor(age_cat, levels = c("20-39", "40-49", "50-59", "60-69", ">70"))) %>%
		mutate(case_cat = factor(case_cat, levels = c("Control", "Cancer"))) %>%
		left_join(tx, by="patient_id") %>%
		mutate(treatment_cat = "No RT/CT") %>%
	  	mutate(treatment_cat = ifelse(prior_rt==1 | prior_chemo==1, "RT/CT", treatment_cat)) %>%
	  	mutate(treatment_cat = ifelse(is.na(treatment_cat), "No RT/CT", treatment_cat)) %>%
	  	mutate(treatment_cat = factor(treatment_cat, levels = c("No RT/CT", "RT/CT")))

plot.0 = ggplot(tmp.1, aes(x = age_cat, y = n, fill = case_cat)) + 
		 geom_boxplot(alpha=1, outlier.size=NA, outlier.shape=NA) + 
		 scale_fill_manual(values = c("#1F6784", "#EAA411")) + 
		 facet_wrap(~facet_1) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.165, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="", y="Number of variants\n") +
		 guides(fill=guide_legend(title=c("Status"))) +
		 coord_cartesian(ylim = c(0,31))
		 
pdf(file="../res/rebuttal/CH_mutations_WBC_Cancer_Control.pdf", width=9, height=5)
print(plot.0)
dev.off()

plot.0 = ggplot(tmp.1, aes(x = age_cat, y = n, fill = treatment_cat)) + 
		 geom_boxplot(alpha=1, outlier.size=NA, outlier.shape=NA) + 
		 scale_fill_manual(values = c("#1F6784", "#EAA411")) + 
		 facet_wrap(~facet_2) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="", y="Number of variants\n") +
		 guides(fill=guide_legend(title=c("Treatment history"))) +
		 coord_cartesian(ylim = c(0,31))
		 
pdf(file="../res/rebuttal/CH_mutations_WBC_RT_CT.pdf", width=9, height=5)
print(plot.0)
dev.off()
