#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
source('config.R')

if (!dir.exists("../res/figureS8")) {
	dir.create("../res/figureS8")
}

#==================================================
# Heat map of genes mutated as VUSo for control and
# non-hypermutated cancer patients
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
variants = left_join(variants, msk_anno %>% dplyr::select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
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

subj_num_smry = variants %>%
				filter(!patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id)) %>%
				group_by(subj_type) %>%
				summarise(subj_num = length(unique(patient_id))) %>%
				ungroup() %>%
				mutate(subj_type_num = paste0(subj_type, "(", subj_num, ")"))

gene_recurrences = variants %>%
 				   filter(!is.na(gene)) %>%
 				   filter(!patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id)) %>%
				   group_by(bio_source, subj_type, gene, is_nonsyn) %>%
 				   summarize(number = n(), num_patient = length(unique(patient_id))) %>%
 				   ungroup %>%
 				   left_join(subj_num_smry) %>%
 				   mutate(percent_patient = round(num_patient / subj_num * 100, 0))
   
gene_list = c()
gene_stats_table_list = list()

for (subj in subj_num_smry$subj_type) {
   top_vuso_gene_ids = gene_recurrences %>%
  					   filter(subj_type == subj, bio_source == "VUSo", is_nonsyn) %>%
   					   group_by(gene) %>%
    				   summarise(all = sum(percent_patient)) %>%
   					   ungroup() %>%
   					   arrange(-all) %>%
   					   dplyr::slice(1:20) %>%
   					   .[["gene"]]
   
   gene_list = c(gene_list, top_vuso_gene_ids)
}
 
for (subj in c("Control", "Breast", "Lung", "Prostate")) {
  subj_smry = gene_recurrences %>%
   			  filter(subj_type == subj, bio_source == "VUSo", is_nonsyn, gene %in% unique(gene_list)) %>%
  			  dplyr::select(gene, num_patient)
   
  gene_stats_table_list[[subj]] = subj_smry
}

gene_stats_smry = data.frame(gene = unique(gene_list))
gene_stats_smry = left_join(gene_stats_smry, gene_stats_table_list[["Control"]], by = "gene") %>%
				  left_join(gene_stats_table_list[["Breast"]], by = "gene") %>%
				  left_join(gene_stats_table_list[["Lung"]], by = "gene") %>%
				  left_join(gene_stats_table_list[["Prostate"]], by = "gene")
colnames(gene_stats_smry) = c("gene", "Control", "Breast", "Lung", "Prostate")
gene_stats_smry[is.na(gene_stats_smry)] = 0
gene_stats_mean = gene_stats_smry %>%
				  dplyr::select(-gene) %>%
 				  rowMeans()
gene_stats_smry = bind_cols(gene_stats_smry, Mean = gene_stats_mean) %>%
				  arrange(-Mean)
				  
gene_stats_smry = gene_stats_smry[1:41,,drop=FALSE]

gene_stats_matrix = as.matrix(gene_stats_smry %>%
					dplyr::select(-gene, -Mean))
rownames(gene_stats_matrix) = gene_stats_smry$gene
idx_names = match(colnames(gene_stats_matrix), subj_num_smry$subj_type)
colnames(gene_stats_matrix) = subj_num_smry$subj_type_num[idx_names]

max_num = max(gene_stats_matrix) + 3
my_palette = colorRampPalette(c("white","#1F6784"))(n = 29)
col_breaks = c(seq(0, 5, length = 10), seq(5.1, max_num, length = 20))

colnames(gene_stats_matrix) = gsub(pattern="(", replacement="\n(n = ", x=colnames(gene_stats_matrix), fixed=TRUE)

pdf(file="../res/figureS8/heatmap_vuso_matched_variants.pdf", width=21, height=14)
corrplot(corr=t(gene_stats_matrix), method="color", is.corr=FALSE, addgrid.col="white",
		 col=colorRampPalette(c(rep("#ffffff",2), "#19547b"))(100), cl.pos = "n",
		 addCoef.col = "black", number.font = 1)

dev.off()


vuso_controls = variants %>%
				filter(subj_type == "Control" & bio_source == "VUSo" & is_nonsyn)

gene_recurrences = vuso_controls %>%
				   count(gene) %>%
				   arrange(desc(n)) %>%
				   dplyr::slice(1:41)
				   
gene_stats_matrix = gene_recurrences[,"n",drop=FALSE]
rownames(gene_stats_matrix) = gene_recurrences$gene

pdf(file="../res/figureS8/heatmap_vuso_variants_controls.pdf", width=21, height=5)
corrplot(corr=t(gene_stats_matrix), method="color", is.corr=FALSE, addgrid.col="white",
		 col=colorRampPalette(c(rep("#ffffff",2), "#19547b"))(100), cl.pos = "n",
		 addCoef.col = "black", number.font = 1)

dev.off()

#==================================================
# Heat map of genes mutated as VUSo for
# hypermutated cancer patients
#==================================================
subj_num_smry = variants %>%
				filter((patient_id %in% hypermutators$patient_id) | (patient_id %in% msi_hypermutators$patient_id)) %>%
				group_by(subj_type) %>%
				summarise(subj_num = length(unique(patient_id))) %>%
				ungroup() %>%
				mutate(subj_type_num = paste0(subj_type, "(", subj_num, ")"))

gene_recurrences = variants %>%
 				   filter(!is.na(gene)) %>%
 				   filter((patient_id %in% hypermutators$patient_id) | (patient_id %in% msi_hypermutators$patient_id)) %>%
				   group_by(bio_source, subj_type, gene, is_nonsyn) %>%
 				   summarize(number = n(), num_patient = length(unique(patient_id))) %>%
 				   ungroup %>%
 				   left_join(subj_num_smry) %>%
 				   mutate(percent_patient = round(num_patient / subj_num * 100, 0))
   
gene_list = c()
gene_stats_table_list = list()

for (subj in subj_num_smry$subj_type) {
   top_vuso_gene_ids = gene_recurrences %>%
  					   filter(subj_type == subj, bio_source == "VUSo", is_nonsyn) %>%
   					   group_by(gene) %>%
    				   summarise(all = sum(percent_patient)) %>%
   					   ungroup() %>%
   					   arrange(-all) %>%
   					   .[["gene"]]
   
   gene_list = c(gene_list, top_vuso_gene_ids)
}
 
for (subj in c("Breast", "Lung", "Prostate")) {
  subj_smry = gene_recurrences %>%
   			  filter(subj_type == subj, bio_source == "VUSo", is_nonsyn, gene %in% unique(gene_list)) %>%
  			  dplyr::select(gene, num_patient)
   
  gene_stats_table_list[[subj]] = subj_smry
}

gene_stats_smry = data.frame(gene = unique(gene_list))
gene_stats_smry = left_join(gene_stats_smry, gene_stats_table_list[["Breast"]], by = "gene") %>%
				  left_join(gene_stats_table_list[["Lung"]], by = "gene") %>%
				  left_join(gene_stats_table_list[["Prostate"]], by = "gene")
colnames(gene_stats_smry) = c("gene", "Breast", "Lung", "Prostate")
gene_stats_smry[is.na(gene_stats_smry)] = 0
gene_stats_mean = gene_stats_smry %>%
				  dplyr::select(-gene) %>%
 				  rowMeans()
gene_stats_smry = bind_cols(gene_stats_smry, Mean = gene_stats_mean) %>%
				  arrange(-Mean)
				  
gene_stats_smry = gene_stats_smry[1:41,,drop=FALSE]

gene_stats_matrix = as.matrix(gene_stats_smry %>%
					dplyr::select(-gene, -Mean))
rownames(gene_stats_matrix) = gene_stats_smry$gene
idx_names = match(colnames(gene_stats_matrix), subj_num_smry$subj_type)
colnames(gene_stats_matrix) = subj_num_smry$subj_type_num[idx_names]

max_num = max(gene_stats_matrix) + 3
my_palette = colorRampPalette(c("white","#1F6784"))(n = 29)
col_breaks = c(seq(0, 5, length = 10), seq(5.1, max_num, length = 20))

colnames(gene_stats_matrix) = gsub(pattern="(", replacement="\n(n = ", x=colnames(gene_stats_matrix), fixed=TRUE)

pdf(file="../res/figureS8/heatmap_vuso_matched_variants_hypermutators.pdf", width=21, height=14)
corrplot(corr=t(gene_stats_matrix), method="color", is.corr=FALSE, addgrid.col="white",
		 col=colorRampPalette(c(rep("#ffffff",2), "#19547b"))(100), cl.pos = "n",
		 addCoef.col = "black", number.font = 1)

dev.off()
