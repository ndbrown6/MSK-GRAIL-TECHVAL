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
# Scatter plot of VAF in WBC and cfDNA
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


variants_nohyper = variants %>%
		   		   filter(!patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id))
cfdna_vs_gdna_vaf_plot = variants_nohyper %>%
						 filter(is_nonsyn) %>%
						 filter(!is.na(afmeancfdna)) %>%
						 filter(bio_source %in% c("WBC_matched", "biopsy_matched","VUSo", "IMPACT-BAM_matched")) %>%
						 mutate(afmeancfdna = afmeancfdna * 100) %>%
						 mutate(afmeangdna = afmeangdna * 100) %>%
						 mutate(afcfdna_nobaq = 100 * (adnobaq+2)/(dpnobaq+4)) %>%
						 mutate(afgdna_nobaq = 100 * (adgdna+2)/(dpgdna+4)) %>%
						 mutate(afcfdna_nobaq_nos = 100 * (adnobaq)/(dpnobaq)) %>%
						 mutate(afgdna_nobaq_nos = 100 * (adgdna)/(dpgdna)) %>%
						 mutate(afcfdna_nobaq_nos = ifelse(is.na(afcfdna_nobaq_nos) | afcfdna_nobaq_nos==0, 0.01, afcfdna_nobaq_nos)) %>%
						 mutate(afgdna_nobaq_nos = ifelse(is.na(afgdna_nobaq_nos) | afgdna_nobaq_nos==0, 0.01, afgdna_nobaq_nos)) %>%
						 mutate(facets_1 = "Mean posterior VAF") %>%
						 mutate(facets_2 = "VAF with pseudocounts without BAQ") %>%
						 mutate(facets_3 = "VAF without pseudocounts without BAQ") %>%
						 mutate(bio_source = case_when(
						 	bio_source == "biopsy_matched" ~ "Biopsy-matched",
						 	bio_source == "IMPACT-BAM_matched" ~ "Biopsy-subthreshold",
						 	bio_source == "VUSo" ~ "VUSo",
						 	bio_source == "WBC_matched" ~ "WBC-matched",
						 ))

cols = c("Biopsy-matched" = "#297AA3",
         "Biopsy-subthreshold" = "#F4AC33",
         "VUSo" = "#3FBC45",
         "WBC-matched"="#D68EAF")

plot.0 = ggplot(cfdna_vs_gdna_vaf_plot, aes(x = afmeancfdna, y = afmeangdna, fill = bio_source)) +
			 geom_abline(linetype = 1, color = "goldenrod3") +
			 geom_point(alpha=1, size=2.5, shape=21, color = "black") +
 			 scale_fill_manual(values = cols) +
 			 facet_wrap(~facets_1) +
 			 theme_bw(base_size=15) +
 			 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.3, 0.8), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
 			 labs(x="\nVAF in cfDNA (%)\n", y="VAF in WBC (%)\n") +
 			 scale_x_log10(
  			 	breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
  			 	labels = function(x) { c("0.01", "0.1", "1", "10", "100") }
  			 ) + 
  			 scale_y_log10(
  			 	breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
  			 	labels = function(x) { c("0.01", "0.1", "1", "10", "100") }
  			 ) +
 			 coord_cartesian(xlim = c(0.01, 100), ylim = c(0.01,100)) +
 			 annotation_logticks() +
			 guides(fill=guide_legend(title=c("Variant category")))
			 

pdf(file="../res/rebuttal/VAF_VAF_Posterior.pdf", width=5, height=6)
print(plot.0)
dev.off()

plot.0 = ggplot(cfdna_vs_gdna_vaf_plot, aes(x = afcfdna_nobaq, y = afgdna_nobaq, fill = bio_source)) +
			 geom_abline(linetype = 1, color = "goldenrod3") +
			 geom_point(alpha=1, size=2.5, shape=21, color = "black") +
 			 scale_fill_manual(values = cols) +
 			 facet_wrap(~facets_2) +
 			 theme_bw(base_size=15) +
 			 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.3, 0.8), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
 			 labs(x="\nVAF in cfDNA (%)\n", y="VAF in WBC (%)\n") +
 			 scale_x_log10(
  			 	breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
  			 	labels = function(x) { c("0.01", "0.1", "1", "10", "100") }
  			 ) + 
  			 scale_y_log10(
  			 	breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
  			 	labels = function(x) { c("0.01", "0.1", "1", "10", "100") }
  			 ) +
 			 coord_cartesian(xlim = c(0.01, 100), ylim = c(0.01,100)) +
 			 annotation_logticks() +
			 guides(fill=guide_legend(title=c("Variant category")))
		 
pdf(file="../res/rebuttal/VAF_VAF_pseudo_no_BAQ.pdf", width=5, height=6)
print(plot.0)
dev.off()

plot.0 = ggplot(cfdna_vs_gdna_vaf_plot, aes(x = afcfdna_nobaq_nos, y = afgdna_nobaq_nos, fill = bio_source)) +
			 geom_abline(linetype = 1, color = "goldenrod3") +
			 geom_point(alpha=1, size=2.5, shape=21, color = "black") +
 			 scale_fill_manual(values = cols) +
 			 facet_wrap(~facets_3) +
 			 theme_bw(base_size=15) +
 			 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.3, 0.8), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
 			 labs(x="\nVAF in cfDNA (%)\n", y="VAF in WBC (%)\n") +
 			 scale_x_log10(
  			 	breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
  			 	labels = function(x) { c("0.01", "0.1", "1", "10", "100") }
  			 ) + 
  			 scale_y_log10(
  			 	breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
  			 	labels = function(x) { c("0.01", "0.1", "1", "10", "100") }
  			 ) +
 			 coord_cartesian(xlim = c(0.01, 100), ylim = c(0.01,100)) +
 			 annotation_logticks() +
			 guides(fill=guide_legend(title=c("Variant category")))
		 
pdf(file="../res/rebuttal/VAF_VAF_nopsedo__noBAQ.pdf", width=5, height=6)
print(plot.0)
dev.off()
