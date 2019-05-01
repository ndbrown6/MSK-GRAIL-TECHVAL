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
# QC metrics
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

#==================================================
# LOD for increasing DNA input to library
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
         			grail == 1 | MSK == 1,
         			patient_id %in% valid_patient_ids) %>%
			 mutate(vtype = "SNV")

indel_plasma = indel_vars %>%
			   filter(ccd == 1,
         			  (c_panel == 1 | panel == 1),
         			  study == "TechVal",
         			  grail == 1 | MSK == 1,
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

all_patient_table = small_vars_plasma %>%
					distinct(subj_type, patient_id)
					
all_patient_table = cbind.data.frame(subj_type = rep(all_patient_table$subj_type, 4),
                                     patient_id = rep(all_patient_table$patient_id, 4),
                                     bio_source = rep(c("WBC_matched",
                                                        "VUSo",
                                                        "biopsy_matched",
                                                        "biopsy_only"),
                                                       each = nrow(all_patient_table)))

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

#==================================================
# LOD for VAF
#==================================================		   		  
tmp = variants %>%
	  filter(bio_source=="biopsy_matched" | bio_source=="IMPACT-BAM_matched") %>%
	  select(af_nobaq, adnobaq, dpnobaq, subj_type, bio_source, patient_id) %>%
	  rename(afnobaq = af_nobaq) %>%
	  mutate(bio_source = ifelse(bio_source=="biopsy_matched", "Biopsy matched", "Biopsy sub-\nthrehold")) %>%
	  mutate(sample_type = "Tissue matched cfDNA variants") %>%
	  filter(bio_source=="Biopsy matched") %>%
	  group_by(patient_id) %>%
	  summarize(LOD_af = min(afnobaq))
tmp2 = qc_metrics %>%
	   filter(Sample_Type=="cfDNA") %>%
	   filter(Tissue!="Healthy") %>%
	   rename(patient_id = Patient_ID)
tmp = full_join(tmp, tmp2, by="patient_id") %>%
	  arrange(Library_preparation_input_ng) %>%
	  filter(!is.na(LOD_af)) %>%
	  mutate(LOD_af = ifelse(LOD_af==0, 0.01, LOD_af)) %>%
	  mutate(Sample_Type = "VAF of tissue matched cfDNA variants")

plot.0 = ggplot(tmp, aes(x=Library_preparation_input_ng, y=LOD_af, fill=Tissue)) +
  		 geom_point(color="black", size=3, shape=21) +
  		 scale_fill_manual(values = c("salmon", "#FDAE61", "#ABDDA4")) +
  		 facet_wrap(~Sample_Type) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="\nAmount of DNA used\nfor library preparation (ng)\n", y="VAF (%)\n") +
		 coord_cartesian(xlim=c(10,80), ylim = c(0.01,100)) +
		 scale_y_log10(
 			 	breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
 			 	labels = function(x) { c("0", "0.1", "1", "10", "100") }
 		 ) +
 		 theme(legend.justification = c(1, 0),
		 	   legend.position = c(1, 0),
		 	   legend.title = element_blank(),
		 	   legend.background = element_blank(),
		 	   legend.text=element_text(size=8))
        	   
pdf(file="../res/rebuttal/LOD_af.pdf", width=6, height=6)
print(plot.0)
dev.off()

#==================================================
# LOD for AD
#==================================================		   		  
tmp = variants %>%
	  filter(bio_source=="biopsy_matched" | bio_source=="IMPACT-BAM_matched") %>%
	  select(af_nobaq, adnobaq, dpnobaq, subj_type, bio_source, patient_id) %>%
	  rename(afnobaq = af_nobaq) %>%
	  mutate(bio_source = ifelse(bio_source=="biopsy_matched", "Biopsy matched", "Biopsy sub-\nthrehold")) %>%
	  mutate(sample_type = "Tissue matched cfDNA variants") %>%
	  filter(bio_source=="Biopsy matched") %>%
	  group_by(patient_id) %>%
	  summarize(LOD_af = min(afnobaq))
LOD_ad = LOD_dp = vector(mode="numeric", length=nrow(tmp))
for (i in 1:nrow(tmp)) {
	index = which(variants$patient_id==tmp$patient_id[i] & variants$af_nobaq==tmp$LOD_af[i])
	index = index[which.min(variants$adnobaq[index])]
	LOD_ad[i] = variants$adnobaq[index]
	LOD_dp[i] = variants$dpnobaq[index]
}
tmp = cbind(tmp, LOD_ad, LOD_dp)
tmp2 = qc_metrics %>%
	   filter(Sample_Type=="cfDNA") %>%
	   filter(Tissue!="Healthy") %>%
	   rename(patient_id = Patient_ID)
tmp = full_join(tmp, tmp2, by="patient_id") %>%
	  arrange(Library_preparation_input_ng) %>%
	  filter(!is.na(LOD_ad)) %>%
	  mutate(LOD_ad = ifelse(LOD_ad==0, 0.1, LOD_ad)) %>%
	  mutate(Sample_Type = "Allelic Depth of tissue matched cfDNA variants")

plot.0 = ggplot(tmp, aes(x=Library_preparation_input_ng, y=LOD_ad, fill=Tissue)) +
  		 geom_point(color="black", size=3, shape=21) +
  		 scale_fill_manual(values = c("salmon", "#FDAE61", "#ABDDA4")) +
  		 facet_wrap(~Sample_Type) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="\nAmount of DNA used\nfor library preparation (ng)\n", y="Base level collapsed depth\n") +
		 coord_cartesian(xlim=c(10,80), ylim = c(0.1,3000)) +
		 scale_y_log10(
 			 	breaks = function(x) { c(0.1, 1, 10, 100, 1000, 3000) },
 			 	labels = function(x) { c("0", "1", "10", "100", "1000", "3000") }
 		 ) +
 		 theme(legend.justification = c(1, 0),
		 	   legend.position = c(1, 0),
		 	   legend.title = element_blank(),
		 	   legend.background = element_blank(),
		 	   legend.text=element_text(size=8))
		 
pdf(file="../res/rebuttal/LOD_ad.pdf", width=6, height=6)
print(plot.0)
dev.off()

#==================================================
# LOD for DP
#==================================================		   		  
tmp = variants %>%
	  filter(bio_source=="biopsy_matched" | bio_source=="IMPACT-BAM_matched") %>%
	  select(af_nobaq, adnobaq, dpnobaq, subj_type, bio_source, patient_id) %>%
	  rename(afnobaq = af_nobaq) %>%
	  mutate(bio_source = ifelse(bio_source=="biopsy_matched", "Biopsy matched", "Biopsy sub-\nthrehold")) %>%
	  mutate(sample_type = "Tissue matched cfDNA variants") %>%
	  filter(bio_source=="Biopsy matched") %>%
	  group_by(patient_id) %>%
	  summarize(LOD_af = min(afnobaq))
LOD_ad = LOD_dp = vector(mode="numeric", length=nrow(tmp))
for (i in 1:nrow(tmp)) {
	index = which(variants$patient_id==tmp$patient_id[i] & variants$af_nobaq==tmp$LOD_af[i])
	index = index[which.min(variants$adnobaq[index])]
	LOD_ad[i] = variants$adnobaq[index]
	LOD_dp[i] = variants$dpnobaq[index]
}
tmp = cbind(tmp, LOD_ad, LOD_dp)
tmp2 = qc_metrics %>%
	   filter(Sample_Type=="cfDNA") %>%
	   filter(Tissue!="Healthy") %>%
	   rename(patient_id = Patient_ID)
tmp = full_join(tmp, tmp2, by="patient_id") %>%
	  arrange(Library_preparation_input_ng) %>%
	  filter(!is.na(LOD_dp)) %>%
	  mutate(LOD_dp = ifelse(LOD_dp==0, 1, LOD_dp)) %>%
	  mutate(Sample_Type = "Total depth of tissue matched cfDNA variants")

plot.0 = ggplot(tmp, aes(x=Library_preparation_input_ng, y=LOD_dp, fill=Tissue)) +
  		 geom_point(color="black", size=3, shape=21) +
  		 scale_fill_manual(values = c("salmon", "#FDAE61", "#ABDDA4")) +
  		 facet_wrap(~Sample_Type) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="\nAmount of DNA used\nfor library preparation (ng)\n", y="Base level collapsed depth\n") +
		 coord_cartesian(xlim=c(10,80), ylim = c(250,10000)) +
		 scale_y_log10(
 			 	breaks = function(x) { c(250, 500, 1000, 2500, 5000, 10000) },
 			 	labels = function(x) { c("250", "500", "1000", "2500", "5000", "10000") }
 		 ) +
 		 theme(legend.justification = c(1, 0),
		 	   legend.position = c(1, 0),
		 	   legend.title = element_blank(),
		 	   legend.background = element_blank(),
		 	   legend.text=element_text(size=8))
		 
pdf(file="../res/rebuttal/LOD_dp.pdf", width=6, height=6)
print(plot.0)
dev.off()

#==================================================
# LOD for VAF and DP
#==================================================		   		  
tmp = variants %>%
	  filter(bio_source=="biopsy_matched" | bio_source=="IMPACT-BAM_matched") %>%
	  select(af_nobaq, adnobaq, dpnobaq, subj_type, bio_source, patient_id) %>%
	  rename(afnobaq = af_nobaq) %>%
	  mutate(bio_source = ifelse(bio_source=="biopsy_matched", "Biopsy matched", "Biopsy sub-\nthrehold")) %>%
	  mutate(sample_type = "Tissue matched cfDNA variants") %>%
	  filter(bio_source=="Biopsy matched") %>%
	  group_by(patient_id) %>%
	  summarize(LOD_af = min(afnobaq))
LOD_ad = LOD_dp = vector(mode="numeric", length=nrow(tmp))
for (i in 1:nrow(tmp)) {
	index = which(variants$patient_id==tmp$patient_id[i] & variants$af_nobaq==tmp$LOD_af[i])
	index = index[which.min(variants$adnobaq[index])]
	LOD_ad[i] = variants$adnobaq[index]
	LOD_dp[i] = variants$dpnobaq[index]
}
tmp = cbind(tmp, LOD_ad, LOD_dp)
tmp2 = qc_metrics %>%
	   filter(Sample_Type=="cfDNA") %>%
	   filter(Tissue!="Healthy") %>%
	   rename(patient_id = Patient_ID)
tmp = full_join(tmp, tmp2, by="patient_id") %>%
	  arrange(Library_preparation_input_ng) %>%
	  filter(!is.na(LOD_af)) %>%
	  filter(!is.na(LOD_ad)) %>%
	  filter(!is.na(LOD_dp)) %>%
	  mutate(Sample_Type = "Total depth of tissue matched cfDNA variants")

plot.0 = ggplot(tmp, aes(x = LOD_af, y = LOD_dp, fill = Tissue)) + 
		 geom_point(alpha=1, aes(size = Library_preparation_input_ng), shape=21) +
		 scale_fill_manual(values = c("salmon", "#FDAE61", "#ABDDA4")) +
		 facet_wrap(~Sample_Type) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="\nVAF (%)\n\n", y="Base level collapsed depth\n") +
		 coord_cartesian(xlim=c(0.01,100), ylim = c(250,10000)) +
		 scale_x_log10(
 			 	breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
 			 	labels = function(x) { c("0", "0.1", "1", "10", "100") }
 		 ) +
		 scale_y_log10(
 			 	breaks = function(x) { c(250, 500, 1000, 2500, 5000, 10000) },
 			 	labels = function(x) { c("250", "500", "1000", "2500", "5000", "10000") }
 		 ) +
		 theme(legend.justification = c(1, 0),
		 	   legend.position = c(1, 0),
		 	   legend.title = element_blank(),
		 	   legend.background = element_blank(),
		 	   legend.text=element_text(size=8))

pdf(file="../res/rebuttal/LOD_DP_versus_LOD_VAF.pdf", width=6, height=6)
print(plot.0)
dev.off()