#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
source('config.R')

if (!dir.exists("../res/figureS4")) {
	dir.create("../res/figureS4")
}

#==================================================
# Barplot of recurrent genes
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
					
variants = label_bio_source(small_vars_plasma) %>%
		   left_join(msk_anno %>% dplyr::select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
                             
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
 		   
# sum the number of times a gene occurs in a bio_source in each tissue type
subj_num_smry = variants %>%
				group_by(subj_type) %>%
    			summarise(subj_num = length(unique(patient_id))) %>%
    			ungroup()
  
gene_recurrences = variants %>%
				   filter(!is.na(gene), bio_source %in% c("biopsy_matched", "biopsy_only") | (bio_source %in% c("IMPACT-BAM_matched", "VUSo") & is_nonsyn)) %>%
 				   group_by(bio_source, subj_type, gene) %>%
 				   summarize(number = n(), num_patient = length(unique(patient_id))) %>%
 				   ungroup %>%
 				   left_join(subj_num_smry) %>%
 				   mutate(percent_patient = round(num_patient / subj_num * 100, 0))
				   
  
gene_list = c()
for (subj in subj_num_smry$subj_type[subj_num_smry$subj_type != "Control"]) {
	top_cancer_gene_ids = gene_recurrences %>%
						  filter(subj_type == subj,
						  bio_source %in% c("biopsy_matched", "biopsy_only")) %>%
 						  group_by(gene) %>%
 						  summarise(all = sum(percent_patient)) %>%
 						  ungroup() %>%
 						  arrange(-all) %>%
 						  slice(1:15) %>%
 						  .[["gene"]]
 	gene_list = c(gene_list, top_cancer_gene_ids)
}
gene_list = unique(gene_list)

top_cancer_genes_ordered = gene_recurrences %>%
 						   filter(bio_source %in% c("biopsy_matched", "biopsy_only"), gene %in% gene_list) %>%
 						   group_by(gene) %>%
 						   summarise(perc = sum(percent_patient)) %>%
 						   ungroup() %>%
 						   arrange(-perc) %>%
 						   .[["gene"]]
   
top_cancer_genes_table = gene_recurrences %>%
 						 filter(bio_source %in% c("biopsy_matched", "IMPACT-BAM_matched", "VUSo"), gene %in% gene_list) %>%
 						 dplyr::select(subj_type, gene, percent_patient, bio_source)
   
top_cancer_genes_table = top_cancer_genes_table %>%
 						 mutate(gene = factor(gene, levels = top_cancer_genes_ordered),
 						 subj_type = factor(subj_type, levels = c("Control", "Breast", "Lung", "Prostate")),
 						 bio_source = factor(bio_source, levels = c("VUSo", "biopsy_matched", "IMPACT-BAM_matched")))
  
bar_fill_cols = c("biopsy_matched" = "#297AA3", "IMPACT-BAM_matched" = "#F4AC33", "VUSo" = "#3FBC45")

plot.0 = top_cancer_genes_table %>%
		 mutate(bio_source = factor(bio_source, levels=c("VUSo", "IMPACT-BAM_matched", "biopsy_matched"), ordered=TRUE)) %>%
		 ggplot(aes(x=gene, y=percent_patient+1, group=bio_source, fill=bio_source)) +
		 geom_bar(stat="identity", color="black") +
		 scale_fill_manual(values = bar_fill_cols) +
		 theme_bw(base_size=15) +
 		 facet_wrap(~subj_type, nrow=4, ncol=1) +
		 labs(x="\n", y="\n% of patients\n") +
		 coord_cartesian(ylim=c(0,60)) +
		 scale_y_continuous(
		 	breaks = function(x) { c(0, 20, 40, 60) },
		 	labels = function(x) { c(0, 20, 40, 60) }
		 ) +
		 theme(legend.position="none",
		 	   axis.text.x = element_text(angle = 90, hjust = 1, size = 10))

pdf(file="../res/figureS4/recurrent_genes_combined_no_hyper.pdf", width=10, height=8.7)
print(plot.0)
dev.off()
