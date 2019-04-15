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
# Add BioRad ddPCR probes to cfDNA VUSo
#==================================================
Table_S7 = read_tsv(file="../res/tables/Table_S7.maf", col_types = cols(.default = col_character()))  %>%
 		   type_convert() %>%
 		   mutate(UID = paste0(Tumor_Sample_Barcode, "_", Chromosome, "_", Start_Position, "_", Hugo_Symbol, "_", HGVSp_Short)) %>%
 		   mutate(UUID = paste0(Hugo_Symbol, "_", HGVSp_Short)) %>%
 		   filter(PHENO == "VUSo")
 		   
BioRad_ddPCR = read.csv(file="../res/etc/Biorad_ddPCR.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
for (i in 1:4) {
	BioRad_ddPCR[seq(from=1, to=nrow(BioRad_ddPCR), by=3) + 1, i] = BioRad_ddPCR[seq(from=1, to=nrow(BioRad_ddPCR), by=3), i]
	BioRad_ddPCR[seq(from=1, to=nrow(BioRad_ddPCR), by=3) + 2, i] = BioRad_ddPCR[seq(from=1, to=nrow(BioRad_ddPCR), by=3), i]
}
FAM = BioRad_ddPCR %>%
	  filter(Fluorophore=="FAM") %>%
	  mutate(UUID = paste0(Gene.Name, "_", Amino.Acid.Change)) %>%
	  select(UUID, Unique.Assay.ID) %>%
	  rename(FAM_Assay_ID = Unique.Assay.ID)
Table_S7 = left_join(Table_S7, FAM, by="UUID")
HEX = BioRad_ddPCR %>%
	  filter(Fluorophore=="HEX") %>%
	  mutate(UUID = paste0(Gene.Name, "_", Amino.Acid.Change)) %>%
	  select(UUID, Unique.Assay.ID) %>%
	  rename(HEX_Assay_ID = Unique.Assay.ID)
Table_S7 = left_join(Table_S7, HEX, by="UUID")
FAM_HEX = BioRad_ddPCR %>%
	  	  filter(Fluorophore=="FAM + HEX") %>%
	  	  mutate(UUID = paste0(Gene.Name, "_", Amino.Acid.Change)) %>%
	  	  select(UUID, Unique.Assay.ID) %>%
	  	  rename(FAM_HEX_Assay_ID = Unique.Assay.ID)
Table_S7 = left_join(Table_S7, FAM_HEX, by="UUID") %>%
		   filter(!duplicated(UID))
write.table(Table_S7, file="../res/tables/Table_S7_vuso_biorad_ddpcr_annot.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

#==================================================
# Scatter plot of maximum and mean VAF for
# biopsy_matched and biopsy_subthreshold mutations
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

tmp.0 = variants %>%
	  filter(subj_type!="Control") %>%
	  filter(bio_source %in% c("biopsy_matched", "IMPACT-BAM_matched", "VUSo", "WBC_matched")) %>%
	  group_by(patient_id) %>%
	  summarize(max_vaf_all = max(af_nobaq))
	  
tmp.1 = variants %>%
	  filter(subj_type!="Control") %>%
	  filter(bio_source == "biopsy_matched" | bio_source == "IMPACT-BAM_matched") %>%
	  group_by(patient_id) %>%
	  summarize(max_vaf_tumor = max(af_nobaq))
	  
tmp = left_join(tmp.0, tmp.1, by="patient_id") %>%
	  mutate(max_vaf_all = ifelse(is.na(max_vaf_all), 0.01, max_vaf_all)) %>%
	  mutate(max_vaf_tumor = ifelse(is.na(max_vaf_tumor), 0.01, max_vaf_tumor)) %>%
	  mutate(facets = "Maximum VAF")
  
plot.0 = ggplot(tmp, aes(x = max_vaf_tumor, y = max_vaf_all)) + 
		 geom_abline(linetype = 1, color = "goldenrod3") +
		 geom_point(shape=21, color="black", fill="#2B83BA", size=3) +
		 geom_hline(yintercept = mean(tmp$max_vaf_all), linetype=2, colour = "goldenrod3") +
		 geom_hline(yintercept = median(tmp$max_vaf_all), linetype=3, colour = "goldenrod3") +
		 geom_vline(xintercept = mean(tmp$max_vaf_tumor), linetype=2, colour = "goldenrod3") +
		 geom_vline(xintercept = median(tmp$max_vaf_tumor), linetype=3, colour = "goldenrod3") +
		 theme_bw(base_size=15) +
 		 labs(x="\nTumor matched variants (%)\n", y="All variants (%)\n") +
 		 scale_x_log10(
 			 	breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
 			 	labels = function(x) { c("0", "0.1", "1", "10", "100") }
 			 ) + 
 		 scale_y_log10(
 			 	breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
 			 	labels = function(x) { c("0", "0.1", "1", "10", "100") }
 			 ) +
		 coord_cartesian(xlim = c(0.01, 100), ylim = c(0.01,100)) +
		 annotation_logticks() +
 		 theme(legend.justification = c(1, 0),
		 	   legend.position = c(1, 0),
		 	   legend.title = element_blank(),
		 	   legend.background = element_blank(),
		 	   legend.text=element_text(size=8)) +
		 facet_wrap(~facets)
pdf(file="../res/rebuttal/Max_VAF_by_Category.pdf", width=6, height=6)
print(plot.0)
dev.off()

tmp.0 = variants %>%
	  filter(subj_type!="Control") %>%
	  filter(bio_source %in% c("biopsy_matched", "IMPACT-BAM_matched")) %>%
	  group_by(patient_id) %>%
	  summarize(max_vaf_tumor = max(af_nobaq))
	  
tmp.1 = variants %>%
	  filter(subj_type!="Control") %>%
	  filter(bio_source %in% c("WBC_matched")) %>%
	  group_by(patient_id) %>%
	  summarize(max_vaf_wbc = max(af_nobaq))
	  
tmp.2 = variants %>%
	  filter(subj_type!="Control") %>%
	  filter(bio_source %in% c("VUSo")) %>%
	  group_by(patient_id) %>%
	  summarize(max_vaf_vuso = max(af_nobaq))
	  
tmp = full_join(tmp.0, tmp.1, by="patient_id") %>%
	  full_join(tmp.2, by="patient_id") %>%
	  mutate(max_vaf_tumor = ifelse(is.na(max_vaf_tumor), 0, max_vaf_tumor)) %>%
	  mutate(max_vaf_wbc = ifelse(is.na(max_vaf_wbc), 0, max_vaf_wbc)) %>%
	  mutate(max_vaf_vuso = ifelse(is.na(max_vaf_vuso), 0, max_vaf_vuso))
	  
tmp = table(apply(tmp[,2:4], 1, which.max))
names(tmp) = c("Tumor\nmacthed", "WBC\nmatched", "VUSo")
tmp = data.frame(tmp) %>%
	  rename(Variant_category = Var1, Frequency = Freq) %>%
	  mutate(Relative_Frequency = 100*Frequency/124) %>%
	  mutate(facets = "Category of maximum VAF")

plot.0 = ggplot(tmp, aes(x = Variant_category, y = Relative_Frequency)) + 
		 geom_bar(stat="identity", color="black", fill="#2B83BA") +
		 theme_bw(base_size=15) +
 		 labs(x="\n\n", y="Fraction of cases (%)\n") +
 		 coord_cartesian(ylim = c(0,100)) +
 		 theme(legend.justification = c(1, 0),
		 	   legend.position = c(1, 0),
		 	   legend.title = element_blank(),
		 	   legend.background = element_blank(),
		 	   legend.text=element_text(size=8)) +
		 facet_wrap(~facets)
pdf(file="../res/rebuttal/Max_VAF_by_Category_Barplot.pdf", width=6, height=6)
print(plot.0)
dev.off()

#==================================================
# Tabulate variants for 6 patients used to test
# assay reproducibility
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

patient_ids = c("MSK-VB-0050", "MSK-VB-0041", "MSK-VL-0028", "MSK-VL-0042", "MSK-VB-0023", "MSK-VL-0038")
z = list()
for (i in 1:length(patient_ids)) {
	tmp = variants %>%
		  filter(patient_id==patient_ids[i]) %>%
		  filter(bio_source %in% c("biopsy_matched", "IMPACT-BAM_matched", "VUSo", "WBC_matched"))
	x = table(tmp$bio_source)
	y = rep(0, 4)
	names(y) = c("biopsy_matched", "IMPACT-BAM_matched", "VUSo", "WBC_matched")
	y[names(x)] = x
	z[[i]] = y
}
z = do.call(rbind, z)
z = cbind(patient_ids, z)
pander(z)

#--------------------------------------------------------------------------------------------------
# &nbsp;    gene.1       filter.1            filter.2        gzero.1   gzero.2  gratio.1   gratio.2 
#--------- -------- ------------------- ------------------- --------- ------------------ ---------- 
# **13**    NOTCH2	PGTKXGDNA_LT_0.79   PASS             	0       0     	1.816   	4.156   
# **20**    TET1	PGTKXGDNA_LT_0.79   PASS		        0       0     	1.337   	2.581   
# **29**    ARID2	PASS          		PGTKXGDNA_LT_0.79   0      	0.5774  2.052       1.106   
# **30**    KMT2D	PASS          		PGTKXGDNA_LT_0.79   0       0     	2.197       1.63   
# **37**    MAPK3	PASS          		PGTKXGDNA_LT_0.79   0      	0.5774  1.758       1.201   
# **38**    MAPK3	PASS          		PGTKXGDNA_LT_0.79   0       0     	2.178       1.459   
# **40**    FANCA	PGTKXGDNA_LT_0.79   PASS             	0       0     	1.71       	2.209  
# **50**    INSR	PGTKXGDNA_LT_0.79   PASS             	0       0     	1.518      	3.367  
# **88**    PIK3CG	PGTKXGDNA_LT_0.79   PASS          		0.5774  0     	1.049      	2.867   
# **110**   TET1	PGTKXGDNA_LT_0.82   PASS          		0.5774  0     	1.626      	2.338   
# **115**   STAG2	PASS          		PGTKXGDNA_LT_0.82   0.5774  1     	2.597      	1.213   
# **134**   TET2	PASS          		PGTKXGDNA_LT_0.82   0       1.633   1.868      	0.7473
# **136**   AMER1	PGTKXGDNA_LT_0.82   PASS             	0       0     	1.74       	2.168   
#--------------------------------------------------------------------------------------------------

#==================================================
# Tabulate %PPA for 6 patients used to test
# assay reproducibility
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

patient_ids = c("MSK-VB-0050", "MSK-VB-0041", "MSK-VL-0028", "MSK-VL-0042", "MSK-VB-0023", "MSK-VL-0038")
bio_labels = c("biopsy_matched", "IMPACT-BAM_matched", "VUSo", "WBC_matched")

vars_rep0 = variants %>%
			filter(patient_id %in% patient_ids) %>%
			filter(bio_source %in% bio_labels)

gdna_params = data_frame(
				subj_type = c("Healthy", "Breast", "Lung", "Prostate"),
				min_p = c(0.8, 0.79, 0.82, 0.79))
clean_target_region = read_tsv(url_target.bed, col_names = c("chrom", "start", "end"))
clean_target_region$in_target = TRUE
techval_repeats = read.csv(url_techval.repeats, header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
				  filter(!patient_id %in% c("MSK-VB-0023_B", "MSK-VL-0028_B", "MSK-VL-0042_B")) %>%
				  mutate(patient_id = gsub(pattern="_A", replacement="", patient_id)) %>%
				  mutate_at(vars(matches("qual")), funs(if_else(is.na(.), 0, .))) %>%
				  mutate(pedge = ifelse(is.na(pedge), 0, pedge), pgtkxgdna = ifelse(is.na(pgtkxgdna), 0, pgtkxgdna)) %>%
				  separate(hgvsp, into = c("p1", "p2", "p3", "p4"), sep = "\\|") %>%
				  separate(is_nonsyn, into = c("t1","t2", "t3", "t4"), sep = "\\|") %>%
				  separate(symbol, into = c("g1", "g2", "g3", "g4"), sep = "\\|") %>%
				  mutate(gene = case_when(
									(.$t1 == TRUE | is.na(.$t2)) ~g1,
									(.$t1 != TRUE & .$t2 == TRUE) ~g2,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 == TRUE) ~g3,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 != TRUE & .$t4 == TRUE) ~g4)) %>%
				  mutate(hgvs_p = case_when(
									(.$t1 == TRUE) ~p1,
									(.$t1 != TRUE & .$t2 == TRUE) ~p2,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 == TRUE) ~p3,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 != TRUE & .$t4 == TRUE) ~p4)) %>%
				  mutate(hgvs_p = sub(".*:", "", hgvs_p),
						 is_nonsyn = ifelse(is.na(hgvs_p), FALSE, TRUE)) %>%
				  select(-t1, -t2, -t3, -t4, -g1, -g2, -g3, -g4, -p1, -p2, -p3, -p4) %>%
				  mutate(filter = "",
					  		 subj_type = case_when(
					  		 				grepl("VB", .$patient_id) ~ "Breast",
					  		 				grepl("VL", .$patient_id) ~ "Lung",
					  		 				grepl("VP", .$patient_id) ~ "Prostate")) %>%
				  left_join(gdna_params) %>%
				  mutate_(filter = ~update_filter(
					  	  filter = filter,
					  	  qualnobaq = qualnobaq,
					  	  pgtkxgdna = pgtkxgdna,
					  	  is_edge = isedge,
					  	  min_p = min_p)) %>%
				  select(-min_p) %>%
				  mutate(loc = str_c(chrom, ":", pos, "_", ref, ">", alt)) %>%
				  mutate(position_orig = pos, ref_orig = ref, alt_orig = alt, position = pos) %>%
				  mutate(gratio = (adnobaq+2)*(dpgdna+4)/((adgdna+2)*(dpnobaq+4)),
				  		 gzero = adgdna/sqrt(adgdna+2)) %>%
				  mutate(start = position_orig,
                         end = position_orig + 1) %>%
                  genome_left_join(clean_target_region, by = c("chrom", "start", "end")) %>%
                  mutate(chrom = chrom.x) %>%
                  select(-c(chrom.x, chrom.y, start.x, end.x, start.y, end.y)) %>%
				  left_join(variants %>% select(patient_id, chrom, position_orig, ref_orig, alt_orig, MSK, grail), by=c("patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
				  mutate(MSK = ifelse(is.na(MSK), 0, MSK)) %>%
				  mutate(grail = 1) %>%
				  mutate(study = "TechVal")
				  
feature_names = intersect(colnames(small_vars_plasma), colnames(techval_repeats))
repeat_variants = bind_rows(small_vars_plasma[,feature_names,drop=FALSE] %>%
							filter(!(patient_id %in% patient_ids)),
							techval_repeats[,feature_names,drop=FALSE])

repeat_variants = label_bio_source(repeat_variants)
repeat_variants = left_join(repeat_variants, msk_anno %>% select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
repeat_variants = repeat_variants %>%
		   				mutate(bio_source = case_when(
		   					   MSK == 1 & grail == 1 ~ "biopsy_matched",
		   					   MSK == 1 & grail == 0 ~ "biopsy_only",
		   					   category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
		   					   category %in% c("germline", "germlineish") ~ "germline",
		   					   category %in% c("blood", "bloodier") ~ "WBC_matched",
		   					   category == "somatic" & `IM-T.alt_count` > bam_read_count_cutoff ~ "IMPACT-BAM_matched",
		   					   category == "somatic" ~ "VUSo",
		   					   TRUE ~ "other"),
			   		  	af_nobaq = round(adnobaq / dpnobaq * 100, 2),
			   		  	af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq))

vars_rep1 = repeat_variants %>%
			filter(patient_id %in% patient_ids)
			
			
for (i in 1:length(patient_ids)) {
	tmp_r1 = vars_rep1 %>%
			 filter(patient_id == patient_ids[i]) %>%
			 mutate(uuid = paste(patient_id, chrom, position_orig, ref_orig, alt_orig, sep="_"))
	tmp_r0 = vars_rep0 %>%
		  	 filter(patient_id == patient_ids[i]) %>%
		  	 filter(is_nonsyn) %>%
		  	 mutate(uuid = paste(patient_id, chrom, position_orig, ref_orig, alt_orig, sep="_")) %>%
		  	 mutate(r = uuid %in% tmp_r1$uuid)
	pander(table(tmp_r0$r, tmp_r0$bio_source))
}
				
#==================================================
# Tabulate %PPA for 3 patients used to test
# assay reproducibility
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

patient_ids = c("MSK-VL-0028", "MSK-VL-0042", "MSK-VB-0023")
bio_labels = c("biopsy_matched", "IMPACT-BAM_matched", "VUSo", "WBC_matched")

vars_rep0 = variants %>%
			filter(patient_id %in% patient_ids) %>%
			filter(bio_source %in% bio_labels)

gdna_params = data_frame(
				subj_type = c("Healthy", "Breast", "Lung", "Prostate"),
				min_p = c(0.8, 0.79, 0.82, 0.79))
clean_target_region = read_tsv(url_target.bed, col_names = c("chrom", "start", "end"))
clean_target_region$in_target = TRUE
techval_repeats = read.csv(url_techval.repeats, header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
				  filter(patient_id %in% c("MSK-VB-0023_B", "MSK-VL-0028_B", "MSK-VL-0042_B")) %>%
				  mutate(patient_id = gsub(pattern="_B", replacement="", patient_id)) %>%
				  mutate_at(vars(matches("qual")), funs(if_else(is.na(.), 0, .))) %>%
				  mutate(pedge = ifelse(is.na(pedge), 0, pedge), pgtkxgdna = ifelse(is.na(pgtkxgdna), 0, pgtkxgdna)) %>%
				  separate(hgvsp, into = c("p1", "p2", "p3", "p4"), sep = "\\|") %>%
				  separate(is_nonsyn, into = c("t1","t2", "t3", "t4"), sep = "\\|") %>%
				  separate(symbol, into = c("g1", "g2", "g3", "g4"), sep = "\\|") %>%
				  mutate(gene = case_when(
									(.$t1 == TRUE | is.na(.$t2)) ~g1,
									(.$t1 != TRUE & .$t2 == TRUE) ~g2,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 == TRUE) ~g3,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 != TRUE & .$t4 == TRUE) ~g4)) %>%
				  mutate(hgvs_p = case_when(
									(.$t1 == TRUE) ~p1,
									(.$t1 != TRUE & .$t2 == TRUE) ~p2,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 == TRUE) ~p3,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 != TRUE & .$t4 == TRUE) ~p4)) %>%
				  mutate(hgvs_p = sub(".*:", "", hgvs_p),
						 is_nonsyn = ifelse(is.na(hgvs_p), FALSE, TRUE)) %>%
				  select(-t1, -t2, -t3, -t4, -g1, -g2, -g3, -g4, -p1, -p2, -p3, -p4) %>%
				  mutate(filter = "",
					  		 subj_type = case_when(
					  		 				grepl("VB", .$patient_id) ~ "Breast",
					  		 				grepl("VL", .$patient_id) ~ "Lung",
					  		 				grepl("VP", .$patient_id) ~ "Prostate")) %>%
				  left_join(gdna_params) %>%
				  mutate_(filter = ~update_filter(
					  	  filter = filter,
					  	  qualnobaq = qualnobaq,
					  	  pgtkxgdna = pgtkxgdna,
					  	  is_edge = isedge,
					  	  min_p = min_p)) %>%
				  select(-min_p) %>%
				  mutate(loc = str_c(chrom, ":", pos, "_", ref, ">", alt)) %>%
				  mutate(position_orig = pos, ref_orig = ref, alt_orig = alt, position = pos) %>%
				  mutate(gratio = (adnobaq+2)*(dpgdna+4)/((adgdna+2)*(dpnobaq+4)),
				  		 gzero = adgdna/sqrt(adgdna+2)) %>%
				  mutate(start = position_orig,
                         end = position_orig + 1) %>%
                  genome_left_join(clean_target_region, by = c("chrom", "start", "end")) %>%
                  mutate(chrom = chrom.x) %>%
                  select(-c(chrom.x, chrom.y, start.x, end.x, start.y, end.y)) %>%
				  left_join(variants %>% select(patient_id, chrom, position_orig, ref_orig, alt_orig, MSK, grail), by=c("patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
				  mutate(MSK = ifelse(is.na(MSK), 0, MSK)) %>%
				  mutate(grail = 1) %>%
				  mutate(study = "TechVal")
				  
feature_names = intersect(colnames(small_vars_plasma), colnames(techval_repeats))
repeat_variants = bind_rows(small_vars_plasma[,feature_names,drop=FALSE] %>%
							filter(!(patient_id %in% patient_ids)),
							techval_repeats[,feature_names,drop=FALSE])

repeat_variants = label_bio_source(repeat_variants)
repeat_variants = left_join(repeat_variants, msk_anno %>% select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
repeat_variants = repeat_variants %>%
		   				mutate(bio_source = case_when(
		   					   MSK == 1 & grail == 1 ~ "biopsy_matched",
		   					   MSK == 1 & grail == 0 ~ "biopsy_only",
		   					   category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
		   					   category %in% c("germline", "germlineish") ~ "germline",
		   					   category %in% c("blood", "bloodier") ~ "WBC_matched",
		   					   category == "somatic" & `IM-T.alt_count` > bam_read_count_cutoff ~ "IMPACT-BAM_matched",
		   					   category == "somatic" ~ "VUSo",
		   					   TRUE ~ "other"),
			   		  	af_nobaq = round(adnobaq / dpnobaq * 100, 2),
			   		  	af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq))

vars_rep1 = repeat_variants %>%
			filter(patient_id %in% patient_ids)
			
			
for (i in 1:length(patient_ids)) {
	tmp_r1 = vars_rep1 %>%
			 filter(patient_id == patient_ids[i]) %>%
			 mutate(uuid = paste(patient_id, chrom, position_orig, ref_orig, alt_orig, sep="_"))
	tmp_r0 = vars_rep0 %>%
		  	 filter(patient_id == patient_ids[i]) %>%
		  	 filter(is_nonsyn) %>%
		  	 mutate(uuid = paste(patient_id, chrom, position_orig, ref_orig, alt_orig, sep="_")) %>%
		  	 mutate(r = uuid %in% tmp_r1$uuid)
	pander(table(tmp_r0$r, tmp_r0$bio_source))
}


#==================================================
# Mutations occurring at dpnobaq >= 10000
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
		   		  
patient_ids = unique(variants$patient_id)

tmp.0 = variants %>%
		filter(bio_source %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
		filter(dpnobaq < 10000) %>%
		group_by(patient_id) %>%
		count(patient_id) %>%
		full_join(data.frame(patient_id = patient_ids), by="patient_id") %>%
		mutate(n_low = ifelse(is.na(n), 0, n)) %>%
		select(-n)
	
tmp.1 = variants %>%
		filter(bio_source %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
		filter(dpnobaq >= 10000) %>%
		group_by(patient_id) %>%
		count(patient_id) %>%
		full_join(data.frame(patient_id = patient_ids), by="patient_id") %>%
		mutate(n_high = ifelse(is.na(n), 0, n)) %>%
		select(-n)
		
tmp = left_join(tmp.0, tmp.1, by="patient_id") %>%
	  mutate(Hypermutated = ifelse(patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id), "Yes", "No")) %>%
	  mutate(n_high = ifelse(n_high == 0, 0.5, n_high)) %>%
	  mutate(n_low = ifelse(n_low == 0, 0.5, n_low)) %>%
	  mutate(facets = "Somatic mutations per sample")

cols = c("Yes"="#D7191C", "No"="#2B83BA")
plot.0 = ggplot(tmp, aes(x = n_high, y = n_low)) +
		 geom_point(alpha=.85, size=3.5, pch = 21, colour = "black", aes(fill=Hypermutated)) +
		 scale_fill_manual(values = cols) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.8, 0.3), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="Number of mutations < 10,000X\n", x="\nNumber of mutations > 10,000X\n") +
		 scale_x_log10(
		 	breaks = function(x) { c(0.5, 1, 10, 100, 200) },
		 	labels = function(x) { c("0", "1", "10", "100", "200") }
		 ) + 
		 scale_y_log10(
		 	breaks = function(x) { c(0.5, 1, 10, 100, 700) },
		 	labels = function(x) { c("0", "1", "10", "100", "700") }
		 ) +
	 	coord_cartesian(xlim = c(0.5, 200), ylim = c(0.5,700)) +
		annotation_logticks() +
		guides(fill=guide_legend(title=c("Hypermutated"))) +
		facet_wrap(~facets)

pdf(file="../res/rebuttal/High_Depth_Mutations.pdf", width=5, height=6)
print(plot.0)
dev.off()


tmp.0 = variants %>%
		filter(dpnobaq>=10000) %>%
		filter(bio_source %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
		mutate(bio_source = ifelse(bio_source %in% c("biopsy_matched", "IMPACT-BAM_matched"), "tumor_matched", bio_source)) %>%
		group_by(patient_id, bio_source) %>%
		count(patient_id, bio_source)
		
patient_ids = unique(tmp.0$patient_id)

tmp.1 = NULL
for (i in c("tumor_matched", "WBC_matched", "VUSo")) {
	tmp.1 = bind_rows(tmp.1, subset(tmp.0, bio_source==i) %>%
							 full_join(data.frame(patient_id = patient_ids), by="patient_id")) %>%
			mutate(n = ifelse(is.na(n), 0, n))
}

cols = c("tumor_matched"="#D7191C",
		 "WBC_matched"="#ABDDA4",
		 "VUSo"="#2B83BA")

plot.0 = ggplot(tmp.1, aes(x=patient_id, y=n)) +
  		 geom_bar(stat="identity", aes(fill=bio_source)) +
  		 theme(axis.text.x=element_text(angle=90)) +
  		 scale_fill_manual(values = cols) +
  		 coord_cartesian(ylim = c(-.1, 120)) +
  		 labs(x="", y="Number of mutations > 10,000X\n")

pdf(file="../res/rebuttal/High_Depth_Mutations_Barplot.pdf", width=10, height=6)
print(plot.0)
dev.off()

tmp.0 = variants %>%
		filter(bio_source %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
		filter(dpnobaq >= 10000) %>%
		group_by(gene) %>%
		count(gene)
		
bed = read_tsv(file=common_bed_annotated_w_introns, col_names = c("chrom", "start", "end", "gene", "intron"))
target_lengths = vector(mode="numeric", length=nrow(tmp.0))
intron_sizes = vector(mode="numeric", length=length(tmp.0))
for (i in 1:nrow(tmp.0)) {
 	index = bed$gene==tmp.0$gene[i]
 	target_lengths[i] = sum(bed[index,3]-bed[index,2])
 	intron_sizes[i] = t(as.matrix(bed[index,3]-bed[index,2])) %*% as.matrix(bed[index,5,drop=FALSE])
}
index = target_lengths == intron_sizes
intron_sizes[index] = 0
target_lengths = target_lengths - intron_sizes

tmp.0 = bind_cols(tmp.0, data.frame(size=target_lengths, facets="Number of mutations per gene"))

plot.0 = ggplot(tmp.0, aes(x = n, y = size, label=gene)) +
		 geom_point(alpha=.85, size=3.5, pch = 21, colour = "black", fill = "#2B83BA") +
	   	 geom_text_repel(data = subset(tmp.0, n >= 4), size=3, force=1, fontface=3)+
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.3, 0.8), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\nNumber of mutations > 10,000X\n", y="Coding target size (bp)\n") +
		 scale_y_log10(
		 	breaks = function(x) { c(50, 100, 1000, 10000, 20000) },
		 	labels = function(x) { c("50", "100", "1000", "10000", "20000") }
		 ) + 
		 scale_x_continuous(
		 	breaks = function(x) { c(0, 2, 4, 6, 8, 10) },
		 	labels = function(x) { c("0", "2", "4", "6", "8", "10") }
		 ) +
	 	 coord_cartesian(xlim = c(0, 10), ylim = c(50,20000)) +
		 annotation_logticks(sides = "l") +
		 facet_wrap(~facets)
		 
pdf(file="../res/rebuttal/High_Depth_Mutations_by_Target_Size.pdf", width=5, height=6)
print(plot.0)
dev.off()

tmp = variants %>%
	  filter(bio_source %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
	  mutate(facets = "Collapsed variant level coverage")
	  
cols = c("biopsy_matched"="#D7191C",
		 "IMPACT-BAM_matched"="#FDAE61",
		 "WBC_matched"="#ABDDA4",
		 "VUSo"="#2B83BA")

plot.0 = ggplot(tmp, aes(x = dpnobaq, y = dpgdna)) +
		 geom_point(alpha=1, size=3.5, pch = 21, colour = "black", aes(fill=bio_source)) +
		 scale_fill_manual(values = cols) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.8, 0.8), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="Variant depth in WBC\n", x="\nVariant depth in cfDNA\n") +
 	 	 coord_cartesian(xlim = c(1, 25000), ylim = c(1,25000)) +
		 guides(fill=guide_legend(title=c("Variant category"))) +
		 facet_wrap(~facets)

pdf(file="../res/rebuttal/High_Depth_cfDNA_vs_gDNA.pdf", width=6, height=6)
print(plot.0)
dev.off()
	  
tmp = left_join(tmp,
				qc_metrics %>%
				filter(Sample_Type=="cfDNA") %>%
				rename(patient_id = Patient_ID),
				by = "patient_id")
				
plot.0 = ggplot(tmp, aes(y = dpnobaq, x = Collapsed_Mean_Coverage)) +
		 geom_point(alpha=1, size=2, pch = 21, colour = "black", aes(fill=bio_source)) +
		 scale_fill_manual(values = cols) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.25, 0.8), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="Depth in cfDNA\n", x="\nMean collapsed coverage\n") +
		 scale_y_log10(
 		 ) + 
		 annotation_logticks() +
		 guides(fill=FALSE) +
		 facet_wrap(~bio_source)

pdf(file="../res/rebuttal/High_Depth_Mutations_by_Mean_Coverage.pdf", width=6, height=6)
print(plot.0)
dev.off()

tmp = tmp %>%
	  mutate(bio_source = ifelse(bio_source %in% c("biopsy_matched","IMPACT-BAM_matched"), "tumor_matched", bio_source)) %>%
	  mutate(cat_1 = case_when(
	 	   		Library_preparation_input_ng<30 ~ "<30",
		   		Library_preparation_input_ng>=30 & Library_preparation_input_ng<45 ~ "30-44",
		   		Library_preparation_input_ng>=45 & Library_preparation_input_ng<60 ~ "45-59",
		   		Library_preparation_input_ng>=60 & Library_preparation_input_ng<75 ~ "60-74",
		   		Library_preparation_input_ng>=75 ~ "75",
		   		)) %>%
	  mutate(cat_2 = factor(paste0(bio_source, ":", cat_1), levels=paste0(rep(c("tumor_matched", "VUSo", "WBC_matched"), each=5), ":", c("<30", "30-44", "45-59", "60-74", "75")))) %>%
	  mutate(dpnobaq = ifelse(dpnobaq>59000, NA, dpnobaq)) %>%
  	  mutate(dpnobaq = ifelse(dpnobaq<100, NA, dpnobaq))

plot.0 = ggplot(tmp, aes(y = dpnobaq, x = cat_2)) +
		 geom_violin(trim=FALSE) +
		 stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
		 facet_wrap(~facets) + 
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="\nInput DNA for library\npreparation (ng)", y="Depth in cfDNA\n") +
		 scale_y_log10(
		 	breaks = function(x) { c(100, 1000, 10000, 50000) },
		 	labels = function(x) { c(100, 1000, 10000, 50000) }
 		 ) +
 		 coord_cartesian(ylim = c(100,60000)) +
		 guides(fill=FALSE)

pdf(file="../res/rebuttal/High_Depth_Mutations_by_input_DNA.pdf", width=10, height=6)
print(plot.0)
dev.off()

cols = c("Tumor matched"="#D7191C",
		 "WBC matched"="#ABDDA4",
		 "VUSo"="#2B83BA")

tmp = variants %>%
	  filter(bio_source %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
	  mutate(bio_source = ifelse(bio_source %in% c("biopsy_matched", "IMPACT-BAM_matched"), "tumor_matched", bio_source)) %>%
	  mutate(bio_source = case_when(
	  	bio_source == "tumor_matched" ~ "Tumor matched",
	  	bio_source == "WBC_matched" ~ "WBC matched",
	  	bio_source == "VUSo" ~ "VUSo"))
		
plot.0 = ggplot(tmp, aes(y = dpnobaq, x = af_nobaq)) +
		 geom_point(alpha=1, size=3, pch = 21, colour = "black", aes(fill=bio_source)) +
		 scale_fill_manual(values = cols) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.25, 0.8), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="Collapsed variant depth in cfDNA\n", x="\nVAF (%)\n") +
		 scale_x_log10(
		 	breaks = function(x) { c(.1, 1, 10, 100) },
		 	labels = function(x) { c(".1", "1", "10", "100") }
 		 ) + 
		 scale_y_log10(
 		 ) + 
		 annotation_logticks(side="bl") +
		 guides(fill=FALSE) +
		 facet_wrap(~bio_source)

pdf(file="../res/rebuttal/High_Depth_Mutations_by_VAF.pdf", width=10, height=6)
print(plot.0)
dev.off()
	  
load("../res/rebuttal/GRAIL/facets/cncf/MSK-VB-0023_MSK-VB-0023-N.Rdata")
pdf(file="../res/rebuttal/MSK-VB-0023_Manhattan_Plot.pdf", width=7, height=1.5)
mat = out2$jointseg
cncf = fit$cncf
dipLogR = fit$dipLogR
par(mar = c(3, 3, 1, 1), mgp = c(2, 0.7, 0))
chr = mat$chrom
len = table(chr)
altcol = rep_len(c("light blue", "gray"), length(len))
chr.col = rep(altcol, len)
nmark = cncf$num.mark
tmp = cumsum(len)
start = c(1, tmp[-length(len)] + 1)
end = tmp
mid = start + len/2
plot(mat$cnlr, type="n", axes = F, ylab = "", xlab="", ylim=c(0,1))
abline(h=.5, lwd=1.5)
labs = names(mid)
labs = sub('21', '', labs)
labs = sub('23', 'X', labs)
axis(side = 1, at = mid, labels=labs, cex.axis = .75, las = 1)
for (i in c(11,14,16,18,20,21,22)) {
	axis(side = 1, at = mid[i], labels=labs[i], cex.axis = .75, las = 1, lwd=-1)
}
tmp = variants %>%
	  filter(bio_source %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
	  mutate(bio_source = ifelse(bio_source %in% c("biopsy_matched", "IMPACT-BAM_matched"), "tumor_matched", bio_source)) %>%
	  filter(patient_id == "MSK-VB-0023") %>%
	  mutate(chrom = gsub("chr", "", chrom, fixed=TRUE)) %>%
	  mutate(chrom = ifelse(chrom=="X", "23", chrom)) %>%
	  mutate(chrom = as.numeric(chrom))

mat = cbind(mat, start=mat$maploc)
mat = cbind(mat, end=c(mat$maploc[2:nrow(mat)], Inf))
load("~/share/usr/anaconda-envs/jrflab-modules-0.1.6/lib/R/library/GAP/data/CytoBand.RData")
for (i in 1:23) {
	index = max(which(mat$chrom==i))
	mat[index,"end"] = max(CytoBand$End[CytoBand$Chromosome==i])
}
index = rep(NA, nrow(tmp))
for (i in 1:nrow(tmp)) {
	indx = which(mat$chrom == tmp$chrom[i] & mat$start<= tmp$position[i] & mat$end>=tmp$position[i])
	if (length(indx)==1) {
		index[i] = indx
	}
}
mat = data.frame(mat, bio_source=rep("", nrow(mat)))
levels(mat$bio_source) = c("", unique(tmp$bio_source))
mat = data.frame(mat, dpnobaq=rep(NA, nrow(mat)))
levels(mat$dpnobaq) = NULL

mat[na.omit(index),"bio_source"] = tmp$bio_source[!is.na(index)]
mat[na.omit(index),"dpnobaq"] = tmp$dpnobaq[!is.na(index)]

x = 1:nrow(mat)
y = ifelse(mat$dpnobaq>10000, .7, 0.2)
z = ifelse(as.character(mat$bio_source) == "tumor_matched", "#D7191C", NA)
z[which(as.character(mat$bio_source)=="WBC_matched")] = "#ABDDA4"
z[which(as.character(mat$bio_source)=="VUSo")] = "#2B83BA"
for (i in c("#2B83BA", "#ABDDA4", "#D7191C")) {
	index = which(z==i)
	points(x[index], y[index], pch="|", col=z[index], cex=2, lwd=.5)
}
rect(xleft=-100000, xright=100000, ybottom=4, ytop=6, col="grey90", lwd=1.5)
box(lwd=1.5)
dev.off()

#==================================================
# Scatter plot of VAF of replicates for MKS-VB-0023
# of mutations at DP > 100000X
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

patient_ids = c("MSK-VB-0023")
vars_rep0 = variants %>%
			filter(patient_id %in% patient_ids)

gdna_params = data_frame(
				subj_type = c("Healthy", "Breast", "Lung", "Prostate"),
				min_p = c(0.8, 0.79, 0.82, 0.79))
clean_target_region = read_tsv(url_target.bed, col_names = c("chrom", "start", "end"))
clean_target_region$in_target = TRUE
techval_repeats = read.csv(url_techval.repeats, header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
				  filter(!patient_id %in% c("MSK-VB-0023_B", "MSK-VL-0028_B", "MSK-VL-0042_B")) %>%
				  mutate(patient_id = gsub(pattern="_A", replacement="", patient_id)) %>%
				  mutate_at(vars(matches("qual")), funs(if_else(is.na(.), 0, .))) %>%
				  mutate(pedge = ifelse(is.na(pedge), 0, pedge), pgtkxgdna = ifelse(is.na(pgtkxgdna), 0, pgtkxgdna)) %>%
				  separate(hgvsp, into = c("p1", "p2", "p3", "p4"), sep = "\\|") %>%
				  separate(is_nonsyn, into = c("t1","t2", "t3", "t4"), sep = "\\|") %>%
				  separate(symbol, into = c("g1", "g2", "g3", "g4"), sep = "\\|") %>%
				  mutate(gene = case_when(
									(.$t1 == TRUE | is.na(.$t2)) ~g1,
									(.$t1 != TRUE & .$t2 == TRUE) ~g2,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 == TRUE) ~g3,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 != TRUE & .$t4 == TRUE) ~g4)) %>%
				  mutate(hgvs_p = case_when(
									(.$t1 == TRUE) ~p1,
									(.$t1 != TRUE & .$t2 == TRUE) ~p2,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 == TRUE) ~p3,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 != TRUE & .$t4 == TRUE) ~p4)) %>%
				  mutate(hgvs_p = sub(".*:", "", hgvs_p),
						 is_nonsyn = ifelse(is.na(hgvs_p), FALSE, TRUE)) %>%
				  select(-t1, -t2, -t3, -t4, -g1, -g2, -g3, -g4, -p1, -p2, -p3, -p4) %>%
				  mutate(filter = "",
					  		 subj_type = case_when(
					  		 				grepl("VB", .$patient_id) ~ "Breast",
					  		 				grepl("VL", .$patient_id) ~ "Lung",
					  		 				grepl("VP", .$patient_id) ~ "Prostate")) %>%
				  left_join(gdna_params) %>%
				  mutate_(filter = ~update_filter(
					  	  filter = filter,
					  	  qualnobaq = qualnobaq,
					  	  pgtkxgdna = pgtkxgdna,
					  	  is_edge = isedge,
					  	  min_p = min_p)) %>%
				  select(-min_p) %>%
				  mutate(loc = str_c(chrom, ":", pos, "_", ref, ">", alt)) %>%
				  mutate(position_orig = pos, ref_orig = ref, alt_orig = alt, position = pos) %>%
				  mutate(gratio = (adnobaq+2)*(dpgdna+4)/((adgdna+2)*(dpnobaq+4)),
				  		 gzero = adgdna/sqrt(adgdna+2)) %>%
				  mutate(start = position_orig,
                         end = position_orig + 1) %>%
                  genome_left_join(clean_target_region, by = c("chrom", "start", "end")) %>%
                  mutate(chrom = chrom.x) %>%
                  select(-c(chrom.x, chrom.y, start.x, end.x, start.y, end.y)) %>%
				  left_join(variants %>% select(patient_id, chrom, position_orig, ref_orig, alt_orig, MSK, grail), by=c("patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
				  mutate(MSK = ifelse(is.na(MSK), 0, MSK)) %>%
				  mutate(grail = 1) %>%
				  mutate(study = "TechVal") %>%
				  filter(in_target)
				  
feature_names = intersect(colnames(small_vars_plasma), colnames(techval_repeats))
repeat_variants = bind_rows(small_vars_plasma[,feature_names,drop=FALSE] %>%
							filter(!(patient_id %in% patient_ids)),
							techval_repeats[,feature_names,drop=FALSE])

repeat_variants = label_bio_source(repeat_variants)
repeat_variants = left_join(repeat_variants, msk_anno %>% select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
repeat_variants = repeat_variants %>%
		   				mutate(bio_source = case_when(
		   					   MSK == 1 & grail == 1 ~ "biopsy_matched",
		   					   MSK == 1 & grail == 0 ~ "biopsy_only",
		   					   category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
		   					   category %in% c("germline", "germlineish") ~ "germline",
		   					   category %in% c("blood", "bloodier") ~ "WBC_matched",
		   					   category == "somatic" & `IM-T.alt_count` > bam_read_count_cutoff ~ "IMPACT-BAM_matched",
		   					   category == "somatic" ~ "VUSo",
		   					   TRUE ~ "other"),
			   		  	af_nobaq = round(adnobaq / dpnobaq * 100, 2),
			   		  	af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq))

vars_rep1 = repeat_variants %>%
			filter(patient_id %in% patient_ids)
			
techval_repeats = read.csv(url_techval.repeats, header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
				  filter(patient_id %in% c("MSK-VB-0023_B", "MSK-VL-0028_B", "MSK-VL-0042_B")) %>%
				  mutate(patient_id = gsub(pattern="_B", replacement="", patient_id)) %>%
				  mutate_at(vars(matches("qual")), funs(if_else(is.na(.), 0, .))) %>%
				  mutate(pedge = ifelse(is.na(pedge), 0, pedge), pgtkxgdna = ifelse(is.na(pgtkxgdna), 0, pgtkxgdna)) %>%
				  separate(hgvsp, into = c("p1", "p2", "p3", "p4"), sep = "\\|") %>%
				  separate(is_nonsyn, into = c("t1","t2", "t3", "t4"), sep = "\\|") %>%
				  separate(symbol, into = c("g1", "g2", "g3", "g4"), sep = "\\|") %>%
				  mutate(gene = case_when(
									(.$t1 == TRUE | is.na(.$t2)) ~g1,
									(.$t1 != TRUE & .$t2 == TRUE) ~g2,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 == TRUE) ~g3,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 != TRUE & .$t4 == TRUE) ~g4)) %>%
				  mutate(hgvs_p = case_when(
									(.$t1 == TRUE) ~p1,
									(.$t1 != TRUE & .$t2 == TRUE) ~p2,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 == TRUE) ~p3,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 != TRUE & .$t4 == TRUE) ~p4)) %>%
				  mutate(hgvs_p = sub(".*:", "", hgvs_p),
						 is_nonsyn = ifelse(is.na(hgvs_p), FALSE, TRUE)) %>%
				  select(-t1, -t2, -t3, -t4, -g1, -g2, -g3, -g4, -p1, -p2, -p3, -p4) %>%
				  mutate(filter = "",
					  		 subj_type = case_when(
					  		 				grepl("VB", .$patient_id) ~ "Breast",
					  		 				grepl("VL", .$patient_id) ~ "Lung",
					  		 				grepl("VP", .$patient_id) ~ "Prostate")) %>%
				  left_join(gdna_params) %>%
				  mutate_(filter = ~update_filter(
					  	  filter = filter,
					  	  qualnobaq = qualnobaq,
					  	  pgtkxgdna = pgtkxgdna,
					  	  is_edge = isedge,
					  	  min_p = min_p)) %>%
				  select(-min_p) %>%
				  mutate(loc = str_c(chrom, ":", pos, "_", ref, ">", alt)) %>%
				  mutate(position_orig = pos, ref_orig = ref, alt_orig = alt, position = pos) %>%
				  mutate(gratio = (adnobaq+2)*(dpgdna+4)/((adgdna+2)*(dpnobaq+4)),
				  		 gzero = adgdna/sqrt(adgdna+2)) %>%
				  mutate(start = position_orig,
                         end = position_orig + 1) %>%
                  genome_left_join(clean_target_region, by = c("chrom", "start", "end")) %>%
                  mutate(chrom = chrom.x) %>%
                  select(-c(chrom.x, chrom.y, start.x, end.x, start.y, end.y)) %>%
				  left_join(variants %>% select(patient_id, chrom, position_orig, ref_orig, alt_orig, MSK, grail), by=c("patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
				  mutate(MSK = ifelse(is.na(MSK), 0, MSK)) %>%
				  mutate(grail = 1) %>%
				  mutate(study = "TechVal") %>%
				  filter(in_target)
				  
feature_names = intersect(colnames(small_vars_plasma), colnames(techval_repeats))
repeat_variants = bind_rows(small_vars_plasma[,feature_names,drop=FALSE] %>%
							filter(!(patient_id %in% patient_ids)),
							techval_repeats[,feature_names,drop=FALSE])

repeat_variants = label_bio_source(repeat_variants)
repeat_variants = left_join(repeat_variants, msk_anno %>% select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
repeat_variants = repeat_variants %>%
		   				mutate(bio_source = case_when(
		   					   MSK == 1 & grail == 1 ~ "biopsy_matched",
		   					   MSK == 1 & grail == 0 ~ "biopsy_only",
		   					   category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
		   					   category %in% c("germline", "germlineish") ~ "germline",
		   					   category %in% c("blood", "bloodier") ~ "WBC_matched",
		   					   category == "somatic" & `IM-T.alt_count` > bam_read_count_cutoff ~ "IMPACT-BAM_matched",
		   					   category == "somatic" ~ "VUSo",
		   					   TRUE ~ "other"),
			   		  	af_nobaq = round(adnobaq / dpnobaq * 100, 2),
			   		  	af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq))

vars_rep2 = repeat_variants %>%
			filter(patient_id %in% patient_ids)
			

all_vars = full_join(vars_rep0 %>% mutate(replicate = 1),
					 vars_rep1 %>% mutate(replicate = 2), by=c("study", "subj_type", "patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
		   select(patient_id, adnobaq.x, adnobaq.y, dpnobaq.x, dpnobaq.y, bio_source.x, bio_source.y) %>%
		   mutate(af_rep1 = 100*adnobaq.x/dpnobaq.x) %>%
		   mutate(af_rep2 = 100*adnobaq.y/dpnobaq.y) %>%
		   mutate(bio_source.x = ifelse(is.na(bio_source.x), "unmatched", bio_source.x)) %>%
		   mutate(bio_source.y = ifelse(is.na(bio_source.y), "unmatched", bio_source.y)) %>%
		   filter(bio_source.x %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo") | bio_source.y %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
		   mutate(adnobaq.x = ifelse(is.na(adnobaq.x), 0, adnobaq.x)) %>%
		   mutate(adnobaq.y = ifelse(is.na(adnobaq.y), 0, adnobaq.y)) %>%
  		   mutate(dpnobaq.x = ifelse(is.na(dpnobaq.x), 0, dpnobaq.x)) %>%
   		   mutate(dpnobaq.y = ifelse(is.na(dpnobaq.y), 0, dpnobaq.y)) %>%
   		   mutate(afnobaq.x = 100*adnobaq.x/dpnobaq.x) %>%
   		   mutate(afnobaq.x = ifelse(afnobaq.x==0 | is.na(afnobaq.x), 0.01, afnobaq.x)) %>%
   		   mutate(afnobaq.y = 100*adnobaq.y/dpnobaq.y) %>%
   		   mutate(afnobaq.y = ifelse(afnobaq.y==0 | is.na(afnobaq.y), 0.01, afnobaq.y))


cols = c(
		 "Not detected in one replicate"="#D7191C",
		 "Not called in one replicate due\nto low quality"="#FDAE61",
		 "Incorrect assignment between replicates"="#ABDDA4",
		 "Called in both replicates"="#2B83BA"
		 )
		 
tmp_vars = all_vars %>% filter(patient_id == patient_ids)

# fix biopsy_matched in y
index = tmp_vars$bio_source.x == "biopsy_matched" & (tmp_vars$bio_source.y !="biopsy_matched" & tmp_vars$bio_source.y != "unmatched")
if (sum(index)!=0) {
	tmp_vars$bio_source.y[index] = "biopsy_matched"
}

# fix biopsy_matched in y
index = tmp_vars$bio_source.x == "IMPACT-BAM_matched" & (tmp_vars$bio_source.y !="IMPACT-BAM_matched" & tmp_vars$bio_source.y != "unmatched")
if (sum(index)!=0) {
	tmp_vars$bio_source.y[index] = "IMPACT-BAM_matched"
}

tmp_vars = tmp_vars %>%
		   mutate(shape = ifelse(bio_source.x=="biopsy_matched" | bio_source.x=="IMPACT-BAM_matched" | bio_source.y=="biopsy_matched" | bio_source.y=="IMPACT-BAM_matched", "Biopsy matched", "Biopsy unmatched")) %>%
		   mutate(fill = ifelse(bio_source.x!=bio_source.y , "Incorrect assignment between replicates", "Called in both replicates")) %>%
		   mutate(fill = ifelse((bio_source.x=="unmatched" & bio_source.y!="unmatched") | (bio_source.x!="unmatched" & bio_source.y=="unmatched"), "Not detected in one replicate", "Called in both replicates")) %>%
		   mutate(fill = ifelse((bio_source.x=="noise" & (bio_source.y!="noise" & bio_source.y!="other")) | (bio_source.y=="noise" & (bio_source.x!="noise" & bio_source.x!="other")), "Not called in one replicate due\nto low quality", fill)) %>%
		   filter(dpnobaq.x>10000) %>%
		   filter(bio_source.x!="noise")
		   

plot.0 = ggplot(tmp_vars, aes(x = afnobaq.x, y = afnobaq.y, shape = shape, fill = fill)) +
		 geom_abline(linetype = 1, color = "goldenrod3") +
		 geom_point(alpha=1, size=3.5) +
		 scale_fill_manual(values = cols) +
		 scale_shape_manual(values = c(24, 21)) +
		 facet_wrap(~patient_id) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.3, 0.8), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\nReplicate 1 (%)\n", y="Replicate 2 (%)\n") +
		 scale_x_log10(
			breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
			labels = function(x) { c("0", "0.1", "1", "10", "100") }
		 ) + 
		 scale_y_log10(
			breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
			labels = function(x) { c("0", "0.1", "1", "10", "100") }
		 ) +
		 coord_cartesian(xlim = c(0.01, 100), ylim = c(0.01,100)) +
		 annotation_logticks() +
		 guides(shape=guide_legend(title=c("Biopsy concordance"), override.aes=list(fill="black"))) +
		 guides(fill=guide_legend(title=c("Variant category")))
		 
pdf(file=paste0("../res/rebuttal/", patient_ids, "_R1_R2_High_Depth.pdf"), width=5.5, height=6.5)
print(plot.0)
dev.off()

all_vars = full_join(vars_rep0 %>% mutate(replicate = 1),
					 vars_rep2 %>% mutate(replicate = 3), by=c("study", "subj_type", "patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
		   select(patient_id, adnobaq.x, adnobaq.y, dpnobaq.x, dpnobaq.y, bio_source.x, bio_source.y) %>%
		   mutate(af_rep1 = 100*adnobaq.x/dpnobaq.x) %>%
		   mutate(af_rep2 = 100*adnobaq.y/dpnobaq.y) %>%
		   mutate(bio_source.x = ifelse(is.na(bio_source.x), "unmatched", bio_source.x)) %>%
		   mutate(bio_source.y = ifelse(is.na(bio_source.y), "unmatched", bio_source.y)) %>%
		   filter(bio_source.x %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo") | bio_source.y %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
		   mutate(adnobaq.x = ifelse(is.na(adnobaq.x), 0, adnobaq.x)) %>%
		   mutate(adnobaq.y = ifelse(is.na(adnobaq.y), 0, adnobaq.y)) %>%
  		   mutate(dpnobaq.x = ifelse(is.na(dpnobaq.x), 0, dpnobaq.x)) %>%
   		   mutate(dpnobaq.y = ifelse(is.na(dpnobaq.y), 0, dpnobaq.y)) %>%
   		   mutate(afnobaq.x = 100*adnobaq.x/dpnobaq.x) %>%
   		   mutate(afnobaq.x = ifelse(afnobaq.x==0 | is.na(afnobaq.x), 0.01, afnobaq.x)) %>%
   		   mutate(afnobaq.y = 100*adnobaq.y/dpnobaq.y) %>%
   		   mutate(afnobaq.y = ifelse(afnobaq.y==0 | is.na(afnobaq.y), 0.01, afnobaq.y))


tmp_vars = all_vars %>% filter(patient_id == patient_ids)
	
# fix biopsy_matched in y
index = tmp_vars$bio_source.x == "biopsy_matched" & (tmp_vars$bio_source.y !="biopsy_matched" & tmp_vars$bio_source.y != "unmatched")
if (sum(index)!=0) {
	tmp_vars$bio_source.y[index] = "biopsy_matched"
}
	
# fix biopsy_matched in y
index = tmp_vars$bio_source.x == "IMPACT-BAM_matched" & (tmp_vars$bio_source.y !="IMPACT-BAM_matched" & tmp_vars$bio_source.y != "unmatched")
if (sum(index)!=0) {
	tmp_vars$bio_source.y[index] = "IMPACT-BAM_matched"
}
	
tmp_vars = tmp_vars %>%
		   mutate(shape = ifelse(bio_source.x=="biopsy_matched" | bio_source.x=="IMPACT-BAM_matched" | bio_source.y=="biopsy_matched" | bio_source.y=="IMPACT-BAM_matched", "Biopsy matched", "Biopsy unmatched")) %>%
		   mutate(fill = ifelse(bio_source.x!=bio_source.y , "Incorrect assignment between replicates", "Called in both replicates")) %>%
		   mutate(fill = ifelse((bio_source.x=="unmatched" & bio_source.y!="unmatched") | (bio_source.x!="unmatched" & bio_source.y=="unmatched"), "Not detected in one replicate", "Called in both replicates")) %>%
		   mutate(fill = ifelse((bio_source.x=="noise" & (bio_source.y!="noise" & bio_source.y!="other")) | (bio_source.y=="noise" & (bio_source.x!="noise" & bio_source.x!="other")), "Not called in one replicate due\nto low quality", fill)) %>%
		   filter(dpnobaq.x>10000) %>%
		   filter(bio_source.x!="noise")
			   	
plot.0 = ggplot(tmp_vars, aes(x = afnobaq.x, y = afnobaq.y, shape = shape, fill = fill)) +
		 geom_abline(linetype = 1, color = "goldenrod3") +
		 geom_point(alpha=1, size=3.5) +
		 scale_fill_manual(values = cols) +
		 scale_shape_manual(values = c(24, 21)) +
		 facet_wrap(~patient_id) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.3, 0.8), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\nReplicate 1 (%)\n", y="Replicate 3 (%)\n") +
		 scale_x_log10(
		 	breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
		 	labels = function(x) { c("0", "0.1", "1", "10", "100") }
		 ) + 
		 scale_y_log10(
		 	breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
		 	labels = function(x) { c("0", "0.1", "1", "10", "100") }
		 ) +
		 coord_cartesian(xlim = c(0.01, 100), ylim = c(0.01,100)) +
		 annotation_logticks() +
		 guides(shape=guide_legend(title=c("Biopsy concordance"), override.aes=list(fill="black"))) +
		 guides(fill=guide_legend(title=c("Variant category")))
	 
pdf(file=paste0("../res/rebuttal/", patient_ids, "_R1_R3_High_Depth.pdf"), width=5.5, height=6.5)
print(plot.0)
dev.off()


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
						 mutate(facets_2 = "VAF with pseudo-counts from pileup without BAQ") %>%
						 mutate(facets_3 = "VAF without pseudo-counts from pileup without BAQ")

cols = c("biopsy_matched" = "#297AA3",
         "IMPACT-BAM_matched" = "#F4AC33",
         "VUSo" = "#3FBC45",
         "WBC_matched"="#D68EAF")

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
			 

pdf(file="../res/rebuttal/VAF_VAF_Posterior.pdf", width=6.5, height=6.5)
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
		 
pdf(file="../res/rebuttal/VAF_VAF_pseudo_no_BAQ.pdf", width=6.5, height=6.5)
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
  			 	labels = function(x) { c("0", "0.1", "1", "10", "100") }
  			 ) + 
  			 scale_y_log10(
  			 	breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
  			 	labels = function(x) { c("0", "0.1", "1", "10", "100") }
  			 ) +
 			 coord_cartesian(xlim = c(0.01, 100), ylim = c(0.01,100)) +
 			 annotation_logticks() +
			 guides(fill=guide_legend(title=c("Variant category")))
		 
pdf(file="../res/rebuttal/VAF_VAF_nopsedo__noBAQ.pdf", width=6.5, height=6.5)
print(plot.0)
dev.off()



#==================================================
# Export mutation summary for MSK-VB-0023
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

patient_ids = c(hypermutators$patient_id, msi_hypermutators$patient_id)
bio_labels = c("biopsy_matched", "IMPACT-BAM_matched", "VUSo")

for (i in 1:length(patient_ids)) {
	vars_hyper = variants %>%
				filter(patient_id %in% patient_ids[i]) %>%
				filter(bio_source %in% bio_labels) %>%
				select(chrom, position, ref, alt, adnobaq, dpnobaq) %>%
				mutate(chrom = gsub("chr", "", chrom)) %>%
				mutate(chrom = as.numeric(ifelse(chrom=="X", 23, chrom)))
	write.table(vars_hyper, file=paste0("../res/rebuttal/GRAIL/summary/tsv/", patient_ids[i], ".tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
}

# #==================================================
# # Heat map of genes mutated as VUSo for control and
# # non-hypermutated cancer patients
# #==================================================
# clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
# 		   type_convert() %>%
# 		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
# 		   
# snv_vars = read_tsv(snv_file$scored, col_types = cols(.default = col_character())) %>%
# 		   type_convert()
# 		   
# indel_vars = read_tsv(indel_file$scored, col_types = cols(.default = col_character())) %>%
# 			 type_convert()
# 
# wbc_stack = read_tsv(wbc_variants$scored, col_types = cols(.default = col_character())) %>%
# 			type_convert()
# 			
# msk_anno = read_tsv(msk_anno_joined, col_types = cols(.default = col_character())) %>%
#   		   type_convert()
# 			
# tracker_grail = read_csv(file=patient_tracker)
# 
# tracker_impact = read_csv(file=impact_tracker)
# 
# valid_patient_ids = tracker_grail %>%
# 					filter(patient_id %in% tracker_impact$patient_id) %>%
# 					.[["patient_id"]]
#   
# indel_vars = indel_vars %>%
# 			 mutate(filter = replace(filter,
#              		patient_id == "MSK-VB-0001" &
#              		gene == "GATA3" &
#              		filter == "PASS",
#              		"CSR_MATCH_ELIMINATED"),
#          	 		ccd = replace(ccd,
#              			   		  patient_id == "MSK-VB-0001" &
#                            		  gene == "GATA3" &
#                            		  filter == "CSR_MATCH_ELIMINATED",
#                            		  0))
# 
# snv_plasma = snv_vars %>%
#   			 filter(ccd == 1,
#          			(c_panel == 1 | panel == 1),
#          			study == "TechVal",
#          			grail == 1 | MSK == 1,
#          			patient_id %in% valid_patient_ids) %>%
# 			 mutate(vtype = "SNV")
# 
# indel_plasma = indel_vars %>%
# 			   filter(ccd == 1,
#          			  (c_panel == 1 | panel == 1),
#          			  study == "TechVal",
#          			  grail == 1 | MSK == 1,
#          			  patient_id %in% valid_patient_ids) %>%
#   			   mutate(vtype = "INDEL",
#          			  altenddistmedian = as.integer(altenddistmedian))
#          			  
# healthy_snv = snv_vars %>%
#   			  filter((c_panel == 1 | panel == 1),
#          			  subj_type == "Healthy",
#          			  grail == 1) %>%
#   			  mutate(vtype = "SNV")
# 
# healthy_indel = indel_vars %>%
#   				filter((c_panel == 1 | panel == 1),
#          			    subj_type == "Healthy",
#          				grail == 1) %>%
#   				mutate(vtype = "INDEL",
#          			   altenddistmedian = as.integer(altenddistmedian))
# 
# small_vars_plasma = full_join(snv_plasma, indel_plasma) %>%
# 					full_join(healthy_snv) %>%
# 					full_join(healthy_indel)
# small_vars_plasma = small_vars_plasma %>%
#   					mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
# 
# small_vars_plasma = small_vars_plasma %>%
# 					mutate(loc = str_c(chrom, ":", position_orig, "_", ref_orig, ">", alt_orig))
# 
# variants = label_bio_source(small_vars_plasma)
# variants = left_join(variants, msk_anno %>% select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
# variants = variants %>%
# 		   mutate(bio_source = case_when(
# 		   					   MSK == 1 & grail == 1 ~ "biopsy_matched",
# 		   					   MSK == 1 & grail == 0 ~ "biopsy_only",
# 		   					   category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
# 		   					   category %in% c("germline", "germlineish") ~ "germline",
# 		   					   category %in% c("blood", "bloodier") ~ "WBC_matched",
# 		   					   category == "somatic" & `IM-T.alt_count` > bam_read_count_cutoff ~ "IMPACT-BAM_matched",
# 		   					   category == "somatic" ~ "VUSo",
# 		   					   TRUE ~ "other"),
# 		   		  af = ifelse(is.na(af), 0, af),
# 		   		  af_nobaq = round(adnobaq / dpnobaq * 100, 2),
# 		   		  af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq))
# 
# variants_nohyper = variants %>%
#   				   filter(!(patient_id %in% hypermutators$patient_id)) %>%
#  				   filter(!(patient_id %in% msi_hypermutators$patient_id)) %>%
#   				   filter(bio_source=="VUSo", is_nonsyn)
# unique_genes = na.omit(unique(variants$gene))
# subj_types = c("Control", "Breast", "Lung", "Prostate")
# n = matrix(NA, nrow=length(unique_genes), ncol=4, dimnames=list(unique_genes, subj_types))
# for (i in 1:length(unique_genes)) {
# 	for (j in 1:length(subj_types)) {
# 		index = which(variants_nohyper$gene==unique_genes[i] & variants_nohyper$subj_type==subj_types[j])
# 		n[i,j] = length(index)
# 	}
# }
# index = order(apply(n, 1, sum), decreasing=TRUE)
# n = n[index,,drop=FALSE]
# pdf(file="../res/etc/heatmap_vuso_matched_variants.pdf", width=14, height=21)
# corrplot(corr=n[1:40,,drop=FALSE], method="color", is.corr=FALSE, addgrid.col="white",
# 		 col=colorRampPalette(c(rep("#ffffff",2), "#19547b"))(100), cl.pos = "n",
# 		 addCoef.col = "black", number.font = 1)
# dev.off()
# 
# 
# variants_hyper = variants %>%
#  				 filter((patient_id %in% hypermutators$patient_id) | (patient_id %in% msi_hypermutators$patient_id)) %>%
#  				 filter(bio_source=="VUSo", is_nonsyn)
# unique_genes = na.omit(unique(variants_hyper$gene))
# subj_types = c("Breast", "Lung", "Prostate")
# n = matrix(NA, nrow=length(unique_genes), ncol=3, dimnames=list(unique_genes, subj_types))
# for (i in 1:length(unique_genes)) {
# 	for (j in 1:length(subj_types)) {
# 		index = which(variants_hyper$gene==unique_genes[i] & variants_hyper$subj_type==subj_types[j])
# 		n[i,j] = length(index)
# 	}
# }
# index = order(apply(n, 1, sum), decreasing=TRUE)
# n = n[index,,drop=FALSE]
# pdf(file="../res/etc/heatmap_vuso_matched_variants_hypermutators.pdf", width=14, height=21)
# corrplot(corr=n[1:40,,drop=FALSE], method="color", is.corr=FALSE, addgrid.col="white",
# 		 col=colorRampPalette(c(rep("#ffffff",2), "#19547b"))(100), cl.pos = "n",
# 		 addCoef.col = "black", number.font = 1)
# dev.off()
# 
#  variants_by_gene = apply(n, 1, sum)
# all_genes = names(variants_by_gene)
# all_genes[all_genes=="FOXL2NB"] = "C3orf72"
# names(variants_by_gene)[names(variants_by_gene)=="FOXL2NB"] = "C3orf72"
# bed = read.csv(file=common_bed_annotated_w_introns, header=FALSE, sep="\t", stringsAsFactors=FALSE)
# target_lengths = vector(mode="numeric", length=length(all_genes))
# intron_sizes = vector(mode="numeric", length=length(all_genes))
# for (i in 1:length(all_genes)) {
#  	index = bed[,4]==all_genes[i]
#  	target_lengths[i] = sum(bed[index,3]-bed[index,2])
#  	intron_sizes[i] = t(as.matrix(bed[index,3]-bed[index,2])) %*% as.matrix(bed[index,5,drop=FALSE])
# }
# index = target_lengths == intron_sizes
# intron_sizes[index] = 0
# target_lengths = target_lengths - intron_sizes
# index = order(target_lengths, decreasing=TRUE)
# variants_by_gene = variants_by_gene[index]
# target_lengths = target_lengths[index]
# 
# pdf(file="../res/etc/target_lengths_number_variants.pdf", width=8, height=8)
# par(mar = c(6.1, 6, 4.1, 1))
# plot(target_lengths, variants_by_gene, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, ylim=c(1,30), xlim=c(0,20000))
# points(target_lengths, variants_by_gene, pch=21, col="black", bg="salmon", cex=1.35)
# axis(1, at = NULL, labels = NULL, cex.axis = 1.5, las = 1, lwd=1.5, lwd.ticks=1.35)
# axis(2, at = c(1,seq(5,30,by=5)), labels=c(1,seq(5,30,by=5)), cex.axis = 1.5, las = 1, lwd=1.5, lwd.ticks=1.35)
# mtext(side = 1, text = "Coding target lengths", line = 4, cex = 1.85)
# mtext(side = 2, text = "Number of variants", line = 4, cex = 1.85)
# index = variants_by_gene>10
# text(target_lengths[index], variants_by_gene[index], labels=names(variants_by_gene[index]), pos=3, cex=.95, font=3)
# text(2000, 28.5, labels="p = 4.44e-16", pos=3, cex=1.25)
# dev.off()
# 
# cats = list(start = c(0, 1000, 2000, 3500),
#  			end = 	c(1000, 2000, 3500, Inf))
# index = list()
# for (i in 1:4) {
# 	index[[5-i]] = which(target_lengths>=cats$start[i] & target_lengths<cats$end[i])
# 	index[[5-i]] = index[[5-i]][order(variants_by_gene[index[[5-i]]], decreasing=TRUE)]
# }
#  
# pdf(file="../res/etc/violin_vuso_matched_variants_hypermutators.pdf", width=14, height=7)
# par(mar = c(6.1, 6, 4.1, 1))
# plot(0, 0, type="n", xlim=c(0.5,4.5), ylim=c(1,30), axes=FALSE, frame.plot=FALSE, xlab="", ylab="")
# for (i in 1:4) {
# 	vioplot(variants_by_gene[index[[i]]], range=1, col=transparentRgb(awtools::mpalette[i], 85), border=awtools::mpalette[i], lwd=3, add=TRUE, at=i, drawRect=FALSE)
# 	points(jitter(rep(i, length(index[[i]])), factor=5), variants_by_gene[index[[i]]], pch=19, col=transparentRgb(awtools::mpalette[i], 185), cex=2)
# }
# 
# axis(1, at = 1:4, labels=c("> 3.5", "2 - 3.5", "1 - 2", "< 1"), cex.axis = 1.5, padj = 0.25)
# axis(2, at = NULL, cex.axis = 1.5, las = 1)
# mtext(side = 2, text = "Number of variants", line = 4, cex = 1.5)
# dev.off()
# 
# x = y = list()
# for (i in 1:4) {
# 	x[[i]] = variants_by_gene[index[[i]]]
# 	y[[i]] = rep(5-i, length(index[[i]]))
# }
# 
# print(cor.test(x=unlist(x), y=unlist(y), method="kendall")$p.value)
# z = jonckheere.test(x=unlist(x), g=unlist(y), alternative="two.sided")
# tmp = cbind(variants_by_gene, round(target_lengths))
# index = order(tmp[,2], decreasing=TRUE)
# tmp = tmp[index,,drop=FALSE]
# pdf(file="../res/etc/heatmap_vuso_matched_variants_hypermutators.pdf", width=21, height=8)
# par(mar = c(7.1, 6.1, 4.1, 1))
# z = barplot(tmp[1:100,1], las=1, col=transparentRgb(awtools::mpalette[3], 85), border=awtools::mpalette[3], names.arg=rep("", 100), axes=FALSE)
# axis(1, at = z, labels=names(unlist(x))[1:100], cex.axis = 1.25, las = 2, font.axis=3, line=1)
# axis(2, at = NULL, cex.axis = 1.5, las = 1, line=-1)
# mtext(side = 2, text = "Number of variants", line = 4, cex = 1.5)
# dev.off()

# #==================================================
# # Table of mutation sources by patient category
# #==================================================
# clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
# 		   type_convert() %>%
# 		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
# 		   
# snv_vars = read_tsv(snv_file$scored, col_types = cols(.default = col_character())) %>%
# 		   type_convert()
# 		   
# indel_vars = read_tsv(indel_file$scored, col_types = cols(.default = col_character())) %>%
# 			 type_convert()
# 
# wbc_stack = read_tsv(wbc_variants$scored, col_types = cols(.default = col_character())) %>%
# 			type_convert()
# 			
# msk_anno = read_tsv(msk_anno_joined, col_types = cols(.default = col_character())) %>%
#   		   type_convert()
# 			
# tracker_grail = read_csv(file=patient_tracker)
# 
# tracker_impact = read_csv(file=impact_tracker)
# 
# bed_file = rtracklayer::import.bed(con=common_bed)
# bed_ranges = GenomicRanges::ranges(bed_file)
# total_bed_Mb = sum(GenomicRanges::width(bed_ranges)) / 1e6
# 
# valid_patient_ids = tracker_grail %>%
# 					filter(patient_id %in% tracker_impact$patient_id) %>%
# 					.[["patient_id"]]
#   
# indel_vars = indel_vars %>%
# 			 mutate(filter = replace(filter,
#              		patient_id == "MSK-VB-0001" &
#              		gene == "GATA3" &
#              		filter == "PASS",
#              		"CSR_MATCH_ELIMINATED"),
#          	 		ccd = replace(ccd,
#              			   		  patient_id == "MSK-VB-0001" &
#                            		  gene == "GATA3" &
#                            		  filter == "CSR_MATCH_ELIMINATED",
#                            		  0))
# 
# snv_plasma = snv_vars %>%
#   			 filter(ccd == 1,
#          			(c_panel == 1 | panel == 1),
#          			study == "TechVal",
#          			grail == 1 | MSK == 1,
#          			patient_id %in% valid_patient_ids) %>%
# 			 mutate(vtype = "SNV")
# 
# indel_plasma = indel_vars %>%
# 			   filter(ccd == 1,
#          			  (c_panel == 1 | panel == 1),
#          			  study == "TechVal",
#          			  grail == 1 | MSK == 1,
#          			  patient_id %in% valid_patient_ids) %>%
#   			   mutate(vtype = "INDEL",
#          			  altenddistmedian = as.integer(altenddistmedian))
#          			  
# healthy_snv = snv_vars %>%
#   			  filter((c_panel == 1 | panel == 1),
#          			  subj_type == "Healthy",
#          			  grail == 1) %>%
#   			  mutate(vtype = "SNV")
# 
# healthy_indel = indel_vars %>%
#   				filter((c_panel == 1 | panel == 1),
#          			    subj_type == "Healthy",
#          				grail == 1) %>%
#   				mutate(vtype = "INDEL",
#          			   altenddistmedian = as.integer(altenddistmedian))
# 
# small_vars_plasma = full_join(snv_plasma, indel_plasma) %>%
# 					full_join(healthy_snv) %>%
# 					full_join(healthy_indel)
# small_vars_plasma = small_vars_plasma %>%
#   					mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
#   					
# small_vars_plasma = small_vars_plasma %>%
# 					mutate(loc = str_c(chrom, ":", position_orig, "_", ref_orig, ">", alt_orig))  					
# 
# variants = label_bio_source(small_vars_plasma)
# variants = left_join(variants, msk_anno %>% select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
# variants = variants %>%
# 		   mutate(bio_source = case_when(
# 		   MSK == 1 & grail == 1 ~ "biopsy_matched",
# 		   MSK == 1 & grail == 0 ~ "biopsy_only",
# 		   category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
# 		   category %in% c("germline", "germlineish") ~ "germline",
# 		   category %in% c("blood", "bloodier") ~ "WBC_matched",
# 		   category == "somatic" & `IM-T.alt_count` > bam_read_count_cutoff ~ "IMPACT-BAM_matched",
#  		   category == "somatic" ~ "VUSo",
#  		   TRUE ~ "other"),
#  		   af = ifelse(is.na(af), 0, af),
#  		   af_nobaq = round(adnobaq / dpnobaq * 100, 2),
#  		   af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq))
# 
# variants_control = variants %>%
# 		   filter(subj_type=="Control") %>%
# 		   filter(is_nonsyn | bio_source=="biopsy_only") %>%
# 		   filter(bio_source %in% c("VUSo","biopsy_matched", "biopsy_only", "IMPACT-BAM_matched", "WBC_matched"))
# 
# variants_nohyper = variants %>%
# 		   filter(subj_type!="Control") %>%
# 		   filter(!patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id)) %>%
# 		   filter(is_nonsyn | bio_source=="biopsy_only") %>%
# 		   filter(bio_source %in% c("VUSo","biopsy_matched", "biopsy_only", "IMPACT-BAM_matched", "WBC_matched"))
# 
# variants_hyper = variants %>%
# 		   filter(subj_type!="Control") %>%
# 		   filter(patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id)) %>%
# 		   filter(is_nonsyn | bio_source=="biopsy_only") %>%
# 		   filter(bio_source %in% c("VUSo","biopsy_matched", "biopsy_only", "IMPACT-BAM_matched", "WBC_matched"))
# 		   
# M1 = matrix(NA, nrow=5, ncol=5)
# 
# M1[1,1] = sum(variants_control$bio_source=="biopsy_only")
# M1[2,1] = sum(variants_nohyper$bio_source=="biopsy_only") + sum(variants_hyper$bio_source=="biopsy_only")
# M1[3,1] = sum(variants_nohyper$bio_source=="biopsy_only" & variants_nohyper$subj_type=="Breast") + sum(variants_hyper$bio_source=="biopsy_only" & variants_hyper$subj_type=="Breast")
# M1[4,1] = sum(variants_nohyper$bio_source=="biopsy_only" & variants_nohyper$subj_type=="Lung") + sum(variants_hyper$bio_source=="biopsy_only" & variants_hyper$subj_type=="Lung")
# M1[5,1] = sum(variants_nohyper$bio_source=="biopsy_only" & variants_nohyper$subj_type=="Prostate") + sum(variants_hyper$bio_source=="biopsy_only" & variants_hyper$subj_type=="Prostate")
# 
# M1[1,2] = sum(variants_control$bio_source=="biopsy_matched")
# M1[2,2] = sum(variants_nohyper$bio_source=="biopsy_matched") + sum(variants_hyper$bio_source=="biopsy_only")
# M1[3,2] = sum(variants_nohyper$bio_source=="biopsy_matched" & variants_nohyper$subj_type=="Breast") + sum(variants_hyper$bio_source=="biopsy_matched" & variants_hyper$subj_type=="Breast")
# M1[4,2] = sum(variants_nohyper$bio_source=="biopsy_matched" & variants_nohyper$subj_type=="Lung") + sum(variants_hyper$bio_source=="biopsy_matched" & variants_hyper$subj_type=="Lung")
# M1[5,2] = sum(variants_nohyper$bio_source=="biopsy_matched" & variants_nohyper$subj_type=="Prostate") + sum(variants_hyper$bio_source=="biopsy_matched" & variants_hyper$subj_type=="Prostate")
# 
# M1[1,3] = sum(variants_control$bio_source=="IMPACT-BAM_matched")
# M1[2,3] = sum(variants_nohyper$bio_source=="IMPACT-BAM_matched") + sum(variants_hyper$bio_source=="IMPACT-BAM_matched")
# M1[3,3] = sum(variants_nohyper$bio_source=="IMPACT-BAM_matched" & variants_nohyper$subj_type=="Breast") + sum(variants_hyper$bio_source=="IMPACT-BAM_matched" & variants_hyper$subj_type=="Breast")
# M1[4,3] = sum(variants_nohyper$bio_source=="IMPACT-BAM_matched" & variants_nohyper$subj_type=="Lung") + sum(variants_hyper$bio_source=="IMPACT-BAM_matched" & variants_hyper$subj_type=="Lung")
# M1[5,3] = sum(variants_nohyper$bio_source=="IMPACT-BAM_matched" & variants_nohyper$subj_type=="Prostate") + sum(variants_hyper$bio_source=="IMPACT-BAM_matched" & variants_hyper$subj_type=="Prostate")
# 
# M1[1,4] = sum(variants_control$bio_source=="WBC_matched")
# M1[2,4] = sum(variants_nohyper$bio_source=="WBC_matched") + sum(variants_hyper$bio_source=="WBC_matched")
# M1[3,4] = sum(variants_nohyper$bio_source=="WBC_matched" & variants_nohyper$subj_type=="Breast") + sum(variants_hyper$bio_source=="WBC_matched" & variants_hyper$subj_type=="Breast")
# M1[4,4] = sum(variants_nohyper$bio_source=="WBC_matched" & variants_nohyper$subj_type=="Lung") + sum(variants_hyper$bio_source=="WBC_matched" & variants_hyper$subj_type=="Lung")
# M1[5,4] = sum(variants_nohyper$bio_source=="WBC_matched" & variants_nohyper$subj_type=="Prostate") + sum(variants_hyper$bio_source=="WBC_matched" & variants_hyper$subj_type=="Prostate")
# 
# M1[1,5] = sum(variants_control$bio_source=="VUSo")
# M1[2,5] = sum(variants_nohyper$bio_source=="VUSo") + sum(variants_hyper$bio_source=="VUSo")
# M1[3,5] = sum(variants_nohyper$bio_source=="VUSo" & variants_nohyper$subj_type=="Breast") + sum(variants_hyper$bio_source=="VUSo" & variants_hyper$subj_type=="Breast")
# M1[4,5] = sum(variants_nohyper$bio_source=="VUSo" & variants_nohyper$subj_type=="Lung") + sum(variants_hyper$bio_source=="VUSo" & variants_hyper$subj_type=="Lung")
# M1[5,5] = sum(variants_nohyper$bio_source=="VUSo" & variants_nohyper$subj_type=="Prostate") + sum(variants_hyper$bio_source=="VUSo" & variants_hyper$subj_type=="Prostate")
# 
# M2 = matrix(NA, nrow=5, ncol=5)
# 
# M2[1,1] = sum(variants_control$bio_source=="biopsy_only")
# M2[2,1] = sum(variants_nohyper$bio_source=="biopsy_only")
# M2[3,1] = sum(variants_nohyper$bio_source=="biopsy_only" & variants_nohyper$subj_type=="Breast")
# M2[4,1] = sum(variants_nohyper$bio_source=="biopsy_only" & variants_nohyper$subj_type=="Lung")
# M2[5,1] = sum(variants_nohyper$bio_source=="biopsy_only" & variants_nohyper$subj_type=="Prostate")
# 
# M2[1,2] = sum(variants_control$bio_source=="biopsy_matched")
# M2[2,2] = sum(variants_nohyper$bio_source=="biopsy_matched")
# M2[3,2] = sum(variants_nohyper$bio_source=="biopsy_matched" & variants_nohyper$subj_type=="Breast")
# M2[4,2] = sum(variants_nohyper$bio_source=="biopsy_matched" & variants_nohyper$subj_type=="Lung")
# M2[5,2] = sum(variants_nohyper$bio_source=="biopsy_matched" & variants_nohyper$subj_type=="Prostate")
# 
# M2[1,3] = sum(variants_control$bio_source=="IMPACT-BAM_matched")
# M2[2,3] = sum(variants_nohyper$bio_source=="IMPACT-BAM_matched")
# M2[3,3] = sum(variants_nohyper$bio_source=="IMPACT-BAM_matched" & variants_nohyper$subj_type=="Breast")
# M2[4,3] = sum(variants_nohyper$bio_source=="IMPACT-BAM_matched" & variants_nohyper$subj_type=="Lung")
# M2[5,3] = sum(variants_nohyper$bio_source=="IMPACT-BAM_matched" & variants_nohyper$subj_type=="Prostate")
# 
# M2[1,4] = sum(variants_control$bio_source=="WBC_matched")
# M2[2,4] = sum(variants_nohyper$bio_source=="WBC_matched")
# M2[3,4] = sum(variants_nohyper$bio_source=="WBC_matched" & variants_nohyper$subj_type=="Breast")
# M2[4,4] = sum(variants_nohyper$bio_source=="WBC_matched" & variants_nohyper$subj_type=="Lung")
# M2[5,4] = sum(variants_nohyper$bio_source=="WBC_matched" & variants_nohyper$subj_type=="Prostate")
# 
# M2[1,5] = sum(variants_control$bio_source=="VUSo")
# M2[2,5] = sum(variants_nohyper$bio_source=="VUSo")
# M2[3,5] = sum(variants_nohyper$bio_source=="VUSo" & variants_nohyper$subj_type=="Breast")
# M2[4,5] = sum(variants_nohyper$bio_source=="VUSo" & variants_nohyper$subj_type=="Lung")
# M2[5,5] = sum(variants_nohyper$bio_source=="VUSo" & variants_nohyper$subj_type=="Prostate")
# 
# M3 = matrix(NA, nrow=5, ncol=5)
# 
# M3[1,1] = 0
# M3[2,1] = sum(variants_hyper$bio_source=="biopsy_only")
# M3[3,1] = sum(variants_hyper$bio_source=="biopsy_only" & variants_hyper$subj_type=="Breast")
# M3[4,1] = sum(variants_hyper$bio_source=="biopsy_only" & variants_hyper$subj_type=="Lung")
# M3[5,1] = sum(variants_hyper$bio_source=="biopsy_only" & variants_hyper$subj_type=="Prostate")
# 
# M3[1,2] = 0
# M3[2,2] = sum(variants_hyper$bio_source=="biopsy_matched")
# M3[3,2] = sum(variants_hyper$bio_source=="biopsy_matched" & variants_hyper$subj_type=="Breast")
# M3[4,2] = sum(variants_hyper$bio_source=="biopsy_matched" & variants_hyper$subj_type=="Lung")
# M3[5,2] = sum(variants_hyper$bio_source=="biopsy_matched" & variants_hyper$subj_type=="Prostate")
# 
# M3[1,3] = 0
# M3[2,3] = sum(variants_hyper$bio_source=="IMPACT-BAM_matched")
# M3[3,3] = sum(variants_hyper$bio_source=="IMPACT-BAM_matched" & variants_hyper$subj_type=="Breast")
# M3[4,3] = sum(variants_hyper$bio_source=="IMPACT-BAM_matched" & variants_hyper$subj_type=="Lung")
# M3[5,3] = sum(variants_hyper$bio_source=="IMPACT-BAM_matched" & variants_hyper$subj_type=="Prostate")
# 
# M3[1,4] = 0
# M3[2,4] = sum(variants_hyper$bio_source=="WBC_matched")
# M3[3,4] = sum(variants_hyper$bio_source=="WBC_matched" & variants_hyper$subj_type=="Breast")
# M3[4,4] = sum(variants_hyper$bio_source=="WBC_matched" & variants_hyper$subj_type=="Lung")
# M3[5,4] = sum(variants_hyper$bio_source=="WBC_matched" & variants_hyper$subj_type=="Prostate")
# 
# M3[1,5] = 0
# M3[2,5] = sum(variants_hyper$bio_source=="VUSo")
# M3[3,5] = sum(variants_hyper$bio_source=="VUSo" & variants_hyper$subj_type=="Breast")
# M3[4,5] = sum(variants_hyper$bio_source=="VUSo" & variants_hyper$subj_type=="Lung")
# M3[5,5] = sum(variants_hyper$bio_source=="VUSo" & variants_hyper$subj_type=="Prostate")
# 
# colnames(M1) = colnames(M2) = colnames(M3) = c("biopsy_only", "biopsy_matched", "biopsy_subthreshold", "WBC_mathed", "VUSo")
# rownames(M1) = rownames(M2) = rownames(M3) = c("Control", "Cancer", "Breast", "Lung", "Prostate")
# M1[2,5] = M1[2,5]-2
# M1[2,2] = M1[2,2]+2
# M2[2,5] = M2[2,5]-1
# M2[2,2] = M2[2,2]+1
# M3[2,5] = M3[2,5]-1
# M3[2,2] = M3[2,2]+1


 
# #==================================================
# # Number of variants in hypermutators / gene
# #==================================================
# all_genes = variants_hyper %>%
#  			filter(bio_source == "VUSo") %>%
#  			.[["gene"]]
# variants_by_gene = table(all_genes)
# bed = read.csv(file=common_bed_annotated, header=FALSE, sep="\t", stringsAsFactors=FALSE)
# target_lengths = vector(mode="numeric", length=length(variants_by_gene))
# for (i in 1:length(variants_by_gene)) {
# 	index = bed[,4]==names(variants_by_gene)[i]
# 	target_lengths[i] = sum(bed[index,3]-bed[index,2])
# }
# index = order(target_lengths)
# variants_by_gene = variants_by_gene[index]
# target_lengths = target_lengths[index]
#  
# pdf(file="../res/etc/target_lengths_number_variants.pdf", width=8, height=8)
# par(mar = c(6.1, 6, 4.1, 1))
# barplot(
# axis(2, at = NULL, labels = NULL, cex.axis = 1.5, las = 1, lwd=1.5, lwd.ticks=1.35)
# mtext(side = 1, text = "Genes", line = 4, cex = 1.85)
# mtext(side = 2, text = "Number of variants", line = 4, cex = 1.85)
# dev.off()


# #==================================================
# # Get list of hypermutators
# #==================================================
# clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
# 		   type_convert() %>%
# 		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
# 		   
# snv_vars = read_tsv(snv_file$scored, col_types = cols(.default = col_character())) %>%
# 		   type_convert()
# 		   
# indel_vars = read_tsv(indel_file$scored, col_types = cols(.default = col_character())) %>%
# 			 type_convert()
# 
# wbc_stack = read_tsv(wbc_variants$scored, col_types = cols(.default = col_character())) %>%
# 			type_convert()
#  			
# tracker_grail = read_csv(file=patient_tracker)
#  
# tracker_impact = read_csv(file=impact_tracker)
# 
# bed_file = rtracklayer::import.bed(con=common_bed)
# bed_ranges = GenomicRanges::ranges(bed_file)
# total_bed_Mb = sum(GenomicRanges::width(bed_ranges)) / 1e6
# 
# valid_patient_ids = tracker_grail %>%
# 					filter(patient_id %in% tracker_impact$patient_id) %>%
# 					.[["patient_id"]]
#   
# indel_vars = indel_vars %>%
# 			  mutate(filter = replace(filter,
#                 	 patient_id == "MSK-VB-0001" &
#               		 gene == "GATA3" &
#               		 filter == "PASS",
#               		 "CSR_MATCH_ELIMINATED"),
#           	 		 ccd = replace(ccd,
#               			   		   patient_id == "MSK-VB-0001" &
#                             	   gene == "GATA3" &
#                             	   filter == "CSR_MATCH_ELIMINATED", 0))
#  
# snv_plasma = snv_vars %>%
#    			 filter(ccd == 1,
#           		   (c_panel == 1 | panel == 1),
#           		   study == "TechVal",
#           		   (grail == 1 | MSK == 1),
#           		   patient_id %in% valid_patient_ids) %>%
#  			 mutate(vtype = "SNV")
#  
# indel_plasma = indel_vars %>%
#  			   filter(ccd == 1,
#        			     (c_panel == 1 | panel == 1),
#           			  study == "TechVal",
#           			  (grail == 1 | MSK == 1),
#           			  patient_id %in% valid_patient_ids) %>%
#    			   mutate(vtype = "INDEL",
#           			  altenddistmedian = as.integer(altenddistmedian))
#           			  
# healthy_snv = snv_vars %>%
#    			  filter((c_panel == 1 | panel == 1),
#           			  subj_type == "Healthy",
#           			  grail == 1) %>%
#    			  mutate(vtype = "SNV")
#  
# healthy_indel = indel_vars %>%
#    				filter((c_panel == 1 | panel == 1),
#           			    subj_type == "Healthy",
#           				grail == 1) %>%
#    				mutate(vtype = "INDEL",
#           			   altenddistmedian = as.integer(altenddistmedian))
#  
# small_vars_plasma = full_join(snv_plasma, indel_plasma) %>%
#  					full_join(healthy_snv) %>%
# 					full_join(healthy_indel)
# small_vars_plasma = small_vars_plasma %>%
#   					mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
#   					
# small_vars_plasma = small_vars_plasma %>%
# 					mutate(loc = str_c(chrom, ":", position_orig, "_", ref_orig, ">", alt_orig))  					
# 
# variants = label_bio_source(small_vars_plasma)
# 
# variants = variants %>%
# 		   mutate(bio_source = case_when(
# 		   					   MSK == 1 & grail == 1 ~ "biopsy_matched",
# 		   					   MSK == 1 & grail == 0 ~ "biopsy_only",
# 		   					   category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
# 		   					   category %in% c("germline", "germlineish") ~ "germline",
# 		   					   category %in% c("blood", "bloodier") ~ "WBC_matched",
# 		   					   category == "somatic" ~ "VUSo",
# 		   					   TRUE ~ "other"),
# 		   		  af = ifelse(is.na(af), 0, af),
# 		   		  af_nobaq = round(adnobaq / dpnobaq * 100, 2),
# 		   		  af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq))
# 
# variants = variants %>%
# 		   filter(is_nonsyn) %>%
#  		   filter(subj_type!="Control") %>%
#  		   filter(bio_source %in% c("VUSo","biopsy_matched"))
#  		   
# summ_per_patient_w_indels = variants %>%
#  				   group_by(patient_id) %>%
#  				   summarize(num_called = n()) %>%
#  				   mutate(TMB = num_called/total_bed_Mb) %>%
#  				   ungroup()
#  
# summ_per_patient_wo_indels = variants %>%
# 				   filter(is_snv) %>%
#  				   group_by(patient_id) %>%
#  				   summarize(num_called = n()) %>%
#  				   mutate(TMB = num_called/total_bed_Mb) %>%
#  				   ungroup()
 
# # mean TMB + 95% CI
# x95 = ci.mean(summ_per_patient_w_indels$TMB, alpha=.05)
# 
# # mean TMB + 99% CI
# x99 = ci.mean(summ_per_patient_w_indels$TMB, alpha=.01)
# 
# # Zehir et al. median TMB + 2*IQR
# xmed = median(summ_per_patient_w_indels$TMB) + 2*iqr(summ_per_patient_w_indels$TMB)
# 
# # 90% percentile
# p90 = quantile(summ_per_patient_w_indels$TMB, .90)
# 
# # 95% percentile
# p95 = quantile(summ_per_patient_w_indels$TMB, .95)
# 
# # 99% percentile
# p99 = quantile(summ_per_patient_w_indels$TMB, .99)
 
# #==================================================
# # extract mutation calls from wbc in all genes
# #==================================================
# load(all_vars_and_clinical)
# 
# #==================================================
# # default filters
# #==================================================
# all_vars = all_vars %>%
# 		   filter(is_patient_valid) %>%
# 		   filter(c_panel) %>%
# 		   filter(!is_hypermutator) %>%
# 		   filter(!is_lowdepth) %>%
# 		   filter(!is_lowqual) %>%
# 		   filter(!is_somatic)
# 
# #==================================================
# # save variants
# #==================================================
# save_vars = all_vars
# 
# #==================================================
# # < 5% recurrence | is_hotspot | Frame-shifting
# # in CH gene
# #==================================================
# n_samples = all_vars %>%
#  			distinct(patient_id) %>%
#   		    count()
# recurrence = all_vars %>%
#    			 group_by(loc_lng) %>%
#    			 count() %>%
#    			 ungroup() %>%
#    			 rename(n_recurrence=n) %>%
#    			 mutate(f_recurrence=n_recurrence/n_samples$n)
# all_vars = all_vars %>%
#    		   left_join(recurrence)
# all_vars = all_vars %>%
#  		   filter(f_recurrence < 0.05 | is_hotspot | (Variant_Classification %in% c("Frame_Shift_Ins", "Frame_Shift_Del", "Nonsense_Mutation", "Nonstop_Mutation") & (SYMBOL %in% chip_genes)))
#   		   
# #==================================================
# # co-occurring indels filter
# #==================================================
# recurrence = all_vars %>%
#    			 group_by(patient_id, loc_srt) %>%
#    			 count() %>%
#    			 ungroup() %>%
#    			 rename(n_indel=n)
# all_vars = all_vars %>%
# 		   left_join(recurrence)
# all_vars = all_vars %>%
# 		   filter(!(n_indel > 1 & indel))
# 		   
# #==================================================
# # variant classes filter
# #==================================================
# all_vars = all_vars %>%
# 		   filter(Variant_Classification!="3'Flank") %>%
#  		   filter(Variant_Classification!="3'UTR") %>%
#  		   filter(Variant_Classification!="5'Flank") %>%
#  		   filter(Variant_Classification!="5'UTR") %>%
#  		   filter(Variant_Classification!="In_Frame_Del") %>%
#  		   filter(Variant_Classification!="In_Frame_Ins") %>%
#  		   filter(Variant_Classification!="Intron") %>%
#  		   filter(Variant_Classification!="RNA") %>%
#  		   filter(Variant_Classification!="Silent")
#  		   
# #==================================================
# # HLA-A and IGR
# #==================================================
# all_vars = all_vars %>%
#  		   filter(SYMBOL!="HLA-A") %>%
#  		   filter(SYMBOL!="")
# 
# #==================================================
# # germ-line filter
# #==================================================
# all_vars = all_vars %>%
#  		   filter((adnobaq/dpnobaq)<=.3 | (Variant_Classification %in% c("Frame_Shift_Ins", "Frame_Shift_Del", "Nonsense_Mutation", "Nonstop_Mutation") & SYMBOL %in% chip_genes))
# 		   
# #==================================================
# # ExAC filter
# #==================================================
# all_vars = all_vars %>%
#  		   filter(!in_exac)
#  		   
# #==================================================
# # gnomAD filter
# #==================================================
# all_vars = all_vars %>%
#  		   filter(!in_gnomad)
#  		   
# burden_healthy = all_vars %>%
#   		   		 filter(subj_type=="Control") %>%
#   		   		 group_by(patient_id) %>%
#   		   		 summarize(num_called = n()) %>%
#   		   		 ungroup() %>%
#   		   		 left_join(clinical)
# patient_ids = (unique(save_vars$patient_id[save_vars$subj_type=="Control"]))[!unique(save_vars$patient_id[save_vars$subj_type=="Control"]) %in% burden_healthy$patient_id]
# index = clinical$patient_id %in% patient_ids
# tmp = clinical[index,,drop=FALSE]
# tmp = cbind(tmp, num_called=rep(0, sum(index)))
# burden_healthy = rbind(burden_healthy, tmp[,colnames(burden_healthy)])
# burden_cancer = all_vars %>%
#   		   		filter(subj_type!="Control") %>%
#   		   		group_by(patient_id) %>%
#   	 	   		summarize(num_called = n()) %>%
#   	 	   		ungroup() %>%
#   	 	   		left_join(clinical)
# patient_ids = (unique(save_vars$patient_id[save_vars$subj_type!="Control"]))[!unique(save_vars$patient_id[save_vars$subj_type!="Control"]) %in% burden_cancer$patient_id]
# index = clinical$patient_id %in% patient_ids
# tmp = clinical[index,,drop=FALSE]
# tmp = cbind(tmp, num_called=rep(0, sum(index)))
# burden_cancer = rbind(burden_cancer, tmp[,colnames(burden_cancer)])
# 
# data_.cancer_all = data.frame(burden_cancer[,c("age", "num_called", "study"),drop=FALSE], row.names=burden_cancer$patient_id)
# data_.healthy_all = data.frame(burden_healthy[,c("age", "num_called", "study"),drop=FALSE], row.names=burden_healthy$patient_id)
# burden_all = bind_rows(data_.cancer_all, data_.healthy_all)
# rownames(burden_all) = c(rownames(data_.cancer_all), rownames(data_.healthy_all))
# burden_all = burden_all[!is.na(burden_all$age),,drop=FALSE]
# 
# burden_healthy = all_vars %>%
# 				 filter(SYMBOL %in% chip_genes) %>%
#   		   		 filter(subj_type=="Control") %>%
#   		   		 group_by(patient_id) %>%
#   		   		 summarize(num_called = n()) %>%
#   		   		 ungroup() %>%
#   		   		 left_join(clinical)
# patient_ids = (unique(save_vars$patient_id[save_vars$subj_type=="Control"]))[!unique(save_vars$patient_id[save_vars$subj_type=="Control"]) %in% burden_healthy$patient_id]
# index = clinical$patient_id %in% patient_ids
# tmp = clinical[index,,drop=FALSE]
# tmp = cbind(tmp, num_called=rep(0, sum(index)))
# burden_healthy = rbind(burden_healthy, tmp[,colnames(burden_healthy)])
# burden_cancer = all_vars %>%
# 				filter(SYMBOL %in% chip_genes) %>%
#   		   		filter(subj_type!="Control") %>%
#   		   		group_by(patient_id) %>%
#   	 	   		summarize(num_called = n()) %>%
#   	 	   		ungroup() %>%
#   	 	   		left_join(clinical)
# patient_ids = (unique(save_vars$patient_id[save_vars$subj_type!="Control"]))[!unique(save_vars$patient_id[save_vars$subj_type!="Control"]) %in% burden_cancer$patient_id]
# index = clinical$patient_id %in% patient_ids
# tmp = clinical[index,,drop=FALSE]
# tmp = cbind(tmp, num_called=rep(0, sum(index)))
# burden_cancer = rbind(burden_cancer, tmp[,colnames(burden_cancer)])
# 
# data_.cancer_ch = data.frame(burden_cancer[,c("age", "num_called", "study"),drop=FALSE], row.names=burden_cancer$patient_id)
# data_.healthy_ch = data.frame(burden_healthy[,c("age", "num_called", "study"),drop=FALSE], row.names=burden_healthy$patient_id)
# burden_ch = bind_rows(data_.cancer_ch, data_.healthy_ch)
# rownames(burden_ch) = c(rownames(data_.cancer_ch), rownames(data_.healthy_ch))
# burden_ch = burden_ch[!is.na(burden_ch$age),,drop=FALSE]
# 
# patient_ids = intersect(rownames(burden_all), rownames(burden_ch))
# burden_ch = burden_ch[patient_ids,,drop=FALSE]
# burden_all = burden_all[patient_ids,,drop=FALSE]
# 
# save(burden_all, burden_ch, file="../res/etc/MSK_TechVal_Number_Mutations.RData")
# 
## Pedram's code
# require(ggplot2)
# require(sandwich)
# require(msm)
# require(cowplot)
# require(RColorBrewer)
# 
# load("../res/etc/MSK_TechVal_Number_Mutations.RData")
# burden_all$ID = row.names(burden_all)
# burden_ch$ID=row.names(burden_ch)
# 
# tx = read.csv(file="../res/etc/prior_tx_techval_0818.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# txburden_ch = merge(burden_ch, tx, by="ID", all=T)
# txburden_ch$cancer_stage = ifelse(grepl("EB|EL", txburden_ch$ID), "early",
#                                   ifelse(grepl("VB|VL|VP", txburden_ch$ID), "metastatic", "control"))
# txburden_ch$cancer_type = ifelse(grepl("EB|VB", txburden_ch$ID), "breast",
#                                  ifelse(grepl("EL|VL", txburden_ch$ID), "lung",
#                                         ifelse(grepl("VP", txburden_ch$ID), "prostate", "control")))
# txburden_ch$prior_rt[is.na(txburden_ch$prior_rt)] = 0
# txburden_ch$prior_chemo[is.na(txburden_ch$prior_chemo)] = 0
# txburden_ch$prior_horm_targeted[is.na(txburden_ch$prior_horm_targeted)] = 0
# txburden_ch$prior_immuno[is.na(txburden_ch$prior_immuno)] = 0
# txburden_ch$study[is.na(txburden_ch$study)] = "TechVal"
# txburden_ch$prior_chemo_rt=ifelse(txburden_ch$prior_chemo==1 & txburden_ch$prior_rt==1, "Chemo+RT",
#                                   ifelse(txburden_ch$prior_chemo==1 & txburden_ch$prior_rt==0, "Chemo-alone",
#                                          ifelse(txburden_ch$prior_chemo==0 & txburden_ch$prior_rt==1, "RT-alone", "None")))
# txburden_ch$prior_chemo_or_rt=ifelse(txburden_ch$prior_chemo==1 | txburden_ch$prior_rt==1, "Chemo/RT", "None")
# txburden_ch$prior_chemo_rt = factor(txburden_ch$prior_chemo_rt, levels = c("None","Chemo-alone","RT-alone", "Chemo+RT"))
# txburden_ch$prior_chemo_or_rt = factor(txburden_ch$prior_chemo_or_rt, levels = c("None","Chemo/RT"))
# 
# p1 = glm(num_called ~ age+cancer_type+study+prior_chemo*prior_rt+prior_immuno+prior_horm_targeted, family="poisson", data=txburden_ch)
# summary(p1)
# 
# p1 = glm(num_called ~ age+study+cancer_type+prior_chemo+prior_rt+prior_immuno+prior_horm_targeted, family="poisson", data=txburden_ch)
# summary(p1)
# 
# cov.p1 = vcovHC(p1, type="HC0")
# std.err = sqrt(diag(cov.p1))
# r.est = cbind(Estimate= coef(p1), "Robust SE" = std.err,
#               "Pr(>|z|)" = 2 * pnorm(abs(coef(p1)/std.err), lower.tail=FALSE),
#               LL = coef(p1) - 1.96 * std.err,
#               UL = coef(p1) + 1.96 * std.err)
# r.est
# ## update p1 model dropping the variables to determine their significane
# p2 = update(p1, . ~ . - cancer_type)
# ## test model differences with chi square test
# anova(p2, p1, test="Chisq")
# 
# ### Individual plots ###
# 
# ## Prior chemotherapy
# p1 = glm(num_called ~ age+prior_chemo, family="poisson", data=txburden_ch)
# summary(p1)
# 
# ## calculate and store predicted values
# p = txburden_ch[which(!is.na(txburden_ch$age)),] # remove NA age cases
# p$phat = predict(p1, type="response")
# p = p[with(p, order(as.factor(prior_chemo),age)), ]
# 
# ## create the plot
# ggplot(p, aes(x = age, y = phat, colour = as.factor(prior_chemo))) +
#        geom_point(aes(y = num_called), alpha=.5, position=position_jitter(h=.2)) +
#        geom_line(size = 1) +
#        scale_colour_brewer(palette = "Dark2") +
#        labs(x = "Age", y = "Expected number of CH mutations")
# 
# ## Prior radiation therapy
# p1 = glm(num_called ~ age+prior_rt, family="poisson", data=txburden_ch)
# summary(p1)
# 
# ## calculate and store predicted values
# p = txburden_ch[which(!is.na(txburden_ch$age)),]
# p$phat = predict(p1, type="response")
# p = p[with(p, order(as.factor(prior_rt),age)), ]
# 
# ## create the plot
# ggplot(p, aes(x = age, y = phat, colour = as.factor(prior_rt))) +
# 	   geom_point(aes(y = num_called), alpha=.5, position=position_jitter(h=.2)) +
# 	   geom_line(size = 1) +
# 	   scale_colour_brewer(palette = "Dark2") +
# 	   labs(x = "Age", y = "Expected number of CH mutations")
# 
# 
# ## Stage
# p1 = glm(num_called ~ age+cancer_stage, family="poisson", data=txburden_ch)
# summary(p1)
# 
# ## calculate and store predicted values
# p = txburden_ch[which(!is.na(txburden_ch$age)),]
# p$phat = predict(p1, type="response")
# p = p[with(p, order(as.factor(cancer_stage),age)), ]
# 
# ## create the plot
# ggplot(p, aes(x = age, y = phat, colour = as.factor(cancer_stage))) +
# 	   geom_point(aes(y = num_called), alpha=.5, position=position_jitter(h=.2)) +
# 	   geom_line(size = 1) +
# 	   scale_colour_brewer(palette = "Dark2") +
# 	   labs(x = "Age", y = "Expected number of CH mutations")
#
#
#
#
# 
# #==================================================
# # compile number is the manuscript
# #==================================================
# clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
# 		   type_convert() %>%
# 		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
# 			
# tracker_grail = read_csv(file=patient_tracker)
# 
# tracker_impact = read_csv(file=impact_tracker)
# 
# valid_patient_ids = tracker_grail %>%
# 					filter(study=="TechVal") %>%
# 					filter(patient_id %in% tracker_impact$patient_id) %>%
# 					.[["patient_id"]]
# 					
# valid_control_ids = tracker_grail %>%
# 					filter(tissue=="Healthy" & study=="Merlin") %>%
# 					.[["patient_id"]]
#  
# clinical = clinical %>%
# 	 	   filter(patient_id %in% valid_patient_ids | patient_id %in% valid_control_ids)
# 	 	   
# median(clinical$age[grep("VB", clinical$patient_id)])
# range(clinical$age[grep("VB", clinical$patient_id)])
# table(clinical$receptor_status[grep("VB", clinical$patient_id)])
# table(clinical$histology[grep("VB", clinical$patient_id)])
# 
# median(clinical$age[grep("VL", clinical$patient_id)])
# range(clinical$age[grep("VL", clinical$patient_id)])
# table(clinical$sex[grep("VL", clinical$patient_id)])
# table(clinical$histology[grep("VL", clinical$patient_id)])
# 
# median(clinical$age[grep("VP", clinical$patient_id)])
# range(clinical$age[grep("VP", clinical$patient_id)])
# table(clinical$cancer_type[grep("VP", clinical$patient_id)])
# 
# tx = read.csv(file="../res/etc/prior_tx_techval_0818.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE, row.names=1)
# tx = subset(tx, rownames(tx) %in% clinical$patient_id)
# 
# tx2 = tx %>%
# 	  filter(cancer_type == "breast")
# sum(tx2$prior_rt==1 | tx2$prior_chemo==1 | tx2$prior_horm_targeted==1 | tx2$prior_immuno==1)/nrow(tx2) * 100
# 
# tx2 = tx %>%
# 	  filter(cancer_type == "lung")
# sum(tx2$prior_rt==1 | tx2$prior_chemo==1 | tx2$prior_horm_targeted==1 | tx2$prior_immuno==1)/nrow(tx2) * 100
# 
# tx2 = tx %>%
# 	  filter(cancer_type == "prostate")
# sum(tx2$prior_rt==1 | tx2$prior_chemo==1 | tx2$prior_horm_targeted==1 | tx2$prior_immuno==1)/nrow(tx2) * 100
# 	  
# 
# sum(clinical$lntxcat1==">=3" & grepl("VB", clinical$patient_id)) / length(grep("VB", clinical$patient_id)) * 100
# sum(clinical$lntxcat1==">=3" & grepl("VL", clinical$patient_id)) / length(grep("VL", clinical$patient_id)) * 100
# sum(clinical$lntxcat1!="." & grepl("VP", clinical$patient_id)) / length(grep("VL", clinical$patient_id)) * 100
# 
# clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
# 		   type_convert() %>%
# 		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
# 		   
# snv_vars = read_tsv(snv_file$scored, col_types = cols(.default = col_character())) %>%
# 		   type_convert()
# 		   
# indel_vars = read_tsv(indel_file$scored, col_types = cols(.default = col_character())) %>%
# 			 type_convert()
# 
# wbc_stack = read_tsv(wbc_variants$scored, col_types = cols(.default = col_character())) %>%
# 			type_convert()
# 			
# msk_anno = read_tsv(msk_anno_joined, col_types = cols(.default = col_character())) %>%
#   		   type_convert()
# 			
# tracker_grail = read_csv(file=patient_tracker)
# 
# tracker_impact = read_csv(file=impact_tracker)
# 
# valid_patient_ids = tracker_grail %>%
# 					filter(patient_id %in% tracker_impact$patient_id) %>%
# 					.[["patient_id"]]
#   
# indel_vars = indel_vars %>%
# 			 mutate(filter = replace(filter,
#              		patient_id == "MSK-VB-0001" &
#              		gene == "GATA3" &
#              		filter == "PASS",
#              		"CSR_MATCH_ELIMINATED"),
#          	 		ccd = replace(ccd,
#              			   		  patient_id == "MSK-VB-0001" &
#                            		  gene == "GATA3" &
#                            		  filter == "CSR_MATCH_ELIMINATED",
#                            		  0))
# 
# snv_plasma = snv_vars %>%
#   			 filter(ccd == 1,
#          			(c_panel == 1 | panel == 1),
#          			study == "TechVal",
#          			grail == 1 | MSK == 1,
#          			patient_id %in% valid_patient_ids) %>%
# 			 mutate(vtype = "SNV")
# 
# indel_plasma = indel_vars %>%
# 			   filter(ccd == 1,
#          			  (c_panel == 1 | panel == 1),
#          			  study == "TechVal",
#          			  grail == 1 | MSK == 1,
#          			  patient_id %in% valid_patient_ids) %>%
#   			   mutate(vtype = "INDEL",
#          			  altenddistmedian = as.integer(altenddistmedian))
#          			  
# healthy_snv = snv_vars %>%
#   			  filter((c_panel == 1 | panel == 1),
#          			  subj_type == "Healthy",
#          			  grail == 1) %>%
#   			  mutate(vtype = "SNV")
# 
# healthy_indel = indel_vars %>%
#   				filter((c_panel == 1 | panel == 1),
#          			    subj_type == "Healthy",
#          				grail == 1) %>%
#   				mutate(vtype = "INDEL",
#          			   altenddistmedian = as.integer(altenddistmedian))
# 
# small_vars_plasma = full_join(snv_plasma, indel_plasma) %>%
# 					full_join(healthy_snv) %>%
# 					full_join(healthy_indel)
# small_vars_plasma = small_vars_plasma %>%
#   					mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
# 
# small_vars_plasma = small_vars_plasma %>%
# 					filter(!(patient_id %in% hypermutators$patient_id)) %>%
# 					filter(!(patient_id %in% msi_hypermutators$patient_id)) %>%
# 					mutate(loc = str_c(chrom, ":", position_orig, "_", ref_orig, ">", alt_orig))
# 
# all_patient_table = small_vars_plasma %>%
# 					distinct(subj_type, patient_id)
# all_patient_table = cbind.data.frame(subj_type = rep(all_patient_table$subj_type, 4),
#                                      patient_id = rep(all_patient_table$patient_id, 4),
#                                      bio_source = rep(c("WBC_matched",
#                                                         "VUSo",
#                                                         "biopsy_matched",
#                                                         "biopsy_only"),
#                                                        each = nrow(all_patient_table)))
# 
# variants = label_bio_source(small_vars_plasma)
# variants = left_join(variants, msk_anno %>% select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
# variants = variants %>%
# 		   mutate(bio_source = case_when(
# 		   					   MSK == 1 & grail == 1 ~ "biopsy_matched",
# 		   					   MSK == 1 & grail == 0 ~ "biopsy_only",
# 		   					   category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
# 		   					   category %in% c("germline", "germlineish") ~ "germline",
# 		   					   category %in% c("blood", "bloodier") ~ "WBC_matched",
# 		   					   category == "somatic" & `IM-T.alt_count` > bam_read_count_cutoff ~ "IMPACT-BAM_matched",
# 		   					   category == "somatic" ~ "VUSo",
# 		   					   TRUE ~ "other"),
# 		   		  af = ifelse(is.na(af), 0, af),
# 		   		  af_nobaq = round(adnobaq / dpnobaq * 100, 2),
# 		   		  af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq))
#  
# variants_ctrl = variants %>%
# 				filter(is_nonsyn) %>%
# 				filter(subj_type=="Control") %>%
# 				filter(bio_source=="VUSo")
# vaf_ctrl = 100*(variants_ctrl$adnobaq / variants_ctrl$dpnobaq)
# sum(vaf_ctrl>1)
# 
# variants_mbc = variants %>%
# 			   filter(is_nonsyn) %>%
# 			   filter(subj_type=="Breast") %>%
# 			   filter(bio_source %in% c("VUSo", "biopsy_matched", "IMPACT-BAM_matched"))
# vaf_mbc = 100*(variants_mbc$adnobaq / variants_mbc$dpnobaq)
# sum(vaf_mbc>=1)/length(vaf_mbc) * 100
# 
# variants_nsclc = variants %>%
# 				 filter(is_nonsyn) %>%
# 			   	 filter(subj_type=="Lung") %>%
# 			     filter(bio_source %in% c("VUSo", "biopsy_matched", "IMPACT-BAM_matched"))
# vaf_nsclc = 100*(variants_nsclc$adnobaq / variants_nsclc$dpnobaq)
# sum(vaf_nsclc>=1)/length(vaf_nsclc) * 100
# 
# variants_crpc = variants %>%
# 				filter(is_nonsyn) %>%
# 			   	filter(subj_type=="Prostate") %>%
# 			    filter(bio_source %in% c("VUSo", "biopsy_matched", "IMPACT-BAM_matched"))
# vaf_crpc = 100*(variants_crpc$adnobaq / variants_crpc$dpnobaq)
# sum(vaf_crpc>=1)/length(vaf_crpc) * 100
# 
# patient_ids = unique(variants %>%
# 		      filter(subj_type!="Control") %>%
# 		      .[["patient_id"]])
# 
# maxvaf_cat = NULL
# for (i in 1:length(patient_ids)) {
# 	variants_per_pat = variants %>%
# 					   filter(is_nonsyn) %>%
# 			   		   filter(patient_id==patient_ids[i]) %>%
# 			    	   filter(bio_source %in% c("VUSo", "biopsy_matched", "IMPACT-BAM_matched"))
# 	category = variants_per_pat$bio_source
# 	vaf = 100*(variants_per_pat$adnobaq / variants_per_pat$dpnobaq)
# 	index = which.max(vaf)
# 	maxvaf_cat = c(maxvaf_cat, category[index])
# }
# 
# small_vars_plasma = full_join(snv_plasma, indel_plasma) %>%
# 					full_join(healthy_snv) %>%
# 					full_join(healthy_indel)
# small_vars_plasma = small_vars_plasma %>%
#   					mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
# 
# small_vars_plasma = small_vars_plasma %>%
# 					mutate(loc = str_c(chrom, ":", position_orig, "_", ref_orig, ">", alt_orig))
# 
# variants = label_bio_source(small_vars_plasma)
# variants = left_join(variants, msk_anno %>% select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
# variants = variants %>%
# 		   mutate(bio_source = case_when(
# 		   					   MSK == 1 & grail == 1 ~ "biopsy_matched",
# 		   					   MSK == 1 & grail == 0 ~ "biopsy_only",
# 		   					   category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
# 		   					   category %in% c("germline", "germlineish") ~ "germline",
# 		   					   category %in% c("blood", "bloodier") ~ "WBC_matched",
# 		   					   category == "somatic" & `IM-T.alt_count` > bam_read_count_cutoff ~ "IMPACT-BAM_matched",
# 		   					   category == "somatic" ~ "VUSo",
# 		   					   TRUE ~ "other"),
# 		   		  af = ifelse(is.na(af), 0, af),
# 		   		  af_nobaq = round(adnobaq / dpnobaq * 100, 2),
# 		   		  af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq))
# 		   		  
# vars_hypermutators = variants %>%
# 					 filter(is_nonsyn) %>%
# 					 filter(patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id)) %>%
# 					 filter(bio_source %in% c("VUSo","biopsy_matched","IMPACT-BAM_matched"))
# 					 
# vars_nonhypermutators = variants %>%
# 						filter(is_nonsyn) %>%
# 						filter(!(patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id))) %>%
# 						filter(subj_type!="Control") %>%
# 					    filter(bio_source %in% c("VUSo","biopsy_matched","IMPACT-BAM_matched"))
# 
# 
# clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
# 		   type_convert() %>%
# 		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
# 		   
# snv_vars = read_tsv(snv_file$scored, col_types = cols(.default = col_character())) %>%
# 		   type_convert()
# 		   
# indel_vars = read_tsv(indel_file$scored, col_types = cols(.default = col_character())) %>%
# 			 type_convert()
# 
# wbc_stack = read_tsv(wbc_variants$scored, col_types = cols(.default = col_character())) %>%
# 			type_convert()
# 			
# tracker_grail = read_csv(file=patient_tracker)
# 
# tracker_impact = read_csv(file=impact_tracker)
# 
# bed_file = rtracklayer::import.bed(con=common_bed)
# bed_ranges = GenomicRanges::ranges(bed_file)
# total_bed_Mb = round(sum(GenomicRanges::width(bed_ranges)) / 1e6, 1)
# 
# valid_patient_ids = tracker_grail %>%
# 					filter(patient_id %in% tracker_impact$patient_id) %>%
# 					.[["patient_id"]]
#   
# indel_vars = indel_vars %>%
# 			 mutate(filter = replace(filter,
#              		patient_id == "MSK-VB-0001" &
#              		gene == "GATA3" &
#              		filter == "PASS",
#              		"CSR_MATCH_ELIMINATED"),
#          	 		ccd = replace(ccd,
#              			   		  patient_id == "MSK-VB-0001" &
#                            		  gene == "GATA3" &
#                            		  filter == "CSR_MATCH_ELIMINATED",
#                            		  0))
# 
# snv_plasma = snv_vars %>%
#   			 filter(ccd == 1,
#          			(c_panel == 1 | panel == 1),
#          			study == "TechVal",
#          			grail == 1,
#          			patient_id %in% valid_patient_ids) %>%
# 			 mutate(vtype = "SNV")
# 
# indel_plasma = indel_vars %>%
# 			   filter(ccd == 1,
#          			  (c_panel == 1 | panel == 1),
#          			  study == "TechVal",
#          			  grail == 1,
#          			  patient_id %in% valid_patient_ids) %>%
#   			   mutate(vtype = "INDEL",
#          			  altenddistmedian = as.integer(altenddistmedian))
#          			  
# healthy_snv = snv_vars %>%
#   			  filter((c_panel == 1 | panel == 1),
#          			  subj_type == "Healthy",
#          			  grail == 1) %>%
#   			  mutate(vtype = "SNV")
# 
# healthy_indel = indel_vars %>%
#   				filter((c_panel == 1 | panel == 1),
#          			    subj_type == "Healthy",
#          				grail == 1) %>%
#   				mutate(vtype = "INDEL",
#          			   altenddistmedian = as.integer(altenddistmedian))
# 
# small_vars_plasma = full_join(snv_plasma, indel_plasma) %>%
# 					full_join(healthy_snv) %>%
# 					full_join(healthy_indel)
# small_vars_plasma = small_vars_plasma %>%
#   					mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
#   					
# small_vars_plasma = small_vars_plasma %>%
# 					filter(!(patient_id %in% hypermutators$patient_id)) %>%
# 					filter(!(patient_id %in% msi_hypermutators$patient_id)) %>%
# 					mutate(loc = str_c(chrom, ":", position_orig, "_", ref_orig, ">", alt_orig))  					
# 
# variants = label_bio_source(small_vars_plasma)
# variants = left_join(variants, msk_anno %>% select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
# variants = variants %>%
# 		   mutate(bio_source = case_when(
# 		   					   MSK == 1 & grail == 1 ~ "biopsy_matched",
# 		   					   MSK == 1 & grail == 0 ~ "biopsy_only",
# 		   					   category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
# 		   					   category %in% c("germline", "germlineish") ~ "germline",
# 		   					   category %in% c("blood", "bloodier") ~ "WBC_matched",
# 		   					   category == "somatic" & `IM-T.alt_count` > bam_read_count_cutoff ~ "IMPACT-BAM_matched",
# 		   					   category == "somatic" ~ "VUSo",
# 		   					   TRUE ~ "other"),
# 		   		  af = ifelse(is.na(af), 0, af),
# 		   		  af_nobaq = round(adnobaq / dpnobaq * 100, 2),
# 		   		  af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq))
# 
# somatic_vars = variants %>%
#  			   filter(is_nonsyn | bio_source=="IMPACT-BAM_matched") %>%
#  			   filter(subj_type!="Control") %>%
#  			   filter(bio_source %in% c("biopsy_matched", "IMPACT-BAM_matched", "VUSo", "WBC_matched"))
#  
# tumor_vars = variants %>%
# 		     filter(is_nonsyn | bio_source=="IMPACT-BAM_matched") %>%
#  			 filter(bio_source %in% c("biopsy_matched", "IMPACT-BAM_matched"))
#  			   
# patient_ids = unique(somatic_vars$patient_id)
# pct_tumor_vars = NULL
# for (i in 1:length(patient_ids)) {
# 	som_vars = somatic_vars %>%
# 			   filter(patient_id == patient_ids[i])
# 	pat_vars = tumor_vars %>%
# 			   filter(patient_id == patient_ids[i])
# 	pct_tumor_vars = c(pct_tumor_vars, 100*nrow(pat_vars)/nrow(som_vars))
# }
# 
# ctrl_vars = variants %>%
# 			filter(is_nonsyn) %>%
# 			filter(subj_type=="Control") %>%
# 			filter(bio_source %in% c("VUSo", "WBC_matched"))
# 
# patient_ids = unique(ctrl_vars$patient_id)
# vars_per_mb = NULL
# for (i in 1:length(patient_ids)) {
# 	tmp_vars = ctrl_vars %>%
# 			   filter(patient_id == patient_ids[i])
# 	vars_per_mb = c(vars_per_mb, nrow(tmp_vars)/total_bed_Mb)
# }
# print(sum(ctrl_vars$bio_source=="WBC_matched")/nrow(ctrl_vars) * 100)
# 
# ctrl_vars = variants %>%
# 			filter(is_nonsyn) %>%
# 			filter(subj_type=="Control") %>%
# 			filter(bio_source %in% c("VUSo", "WBC_matched"))
# 
# patient_ids = unique(ctrl_vars$patient_id)
# vars_per_mb = NULL
# for (i in 1:length(patient_ids)) {
# 	tmp_vars = ctrl_vars %>%
# 			   filter(patient_id == patient_ids[i])
# 	vars_per_mb = c(vars_per_mb, nrow(tmp_vars)/total_bed_Mb)
# }
# print(sum(ctrl_vars$bio_source=="WBC_matched")/nrow(ctrl_vars) * 100)
# 
# som_vars = variants %>%
# 		   filter(is_nonsyn) %>%
# 		   filter(subj_type!="Control") %>%
# 		   filter(bio_source %in% c("VUSo", "WBC_matched", "biopsy_matched", "IMPACT-BAM_matched"))
# print(sum(som_vars$bio_source=="WBC_matched")/nrow(som_vars) * 100)
# print(sum(som_vars$bio_source %in% c("biopsy_matched", "IMPACT-BAM_matched"))/nrow(som_vars) * 100)
# 
# wbc_vars = variants %>%
# 		   filter(is_nonsyn) %>%
# 		   filter(subj_type=="Control") %>%
# 		   filter(bio_source=="WBC_matched")
# print(length(unique(wbc_vars$patient_id))/47 * 100)
# 
# wbc_vars = variants %>%
# 		   filter(is_nonsyn) %>%
# 		   filter(subj_type!="Control") %>%
# 		   filter(bio_source=="WBC_matched")
# print(length(unique(wbc_vars$patient_id))/124 * 100)
# 
# 
# 
# 
# 
# 
# clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
# 		   type_convert() %>%
# 		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
# 		   
# snv_vars = read_tsv(snv_file$scored, col_types = cols(.default = col_character())) %>%
# 		   type_convert()
# 		   
# indel_vars = read_tsv(indel_file$scored, col_types = cols(.default = col_character())) %>%
# 			 type_convert()
# 
# wbc_stack = read_tsv(wbc_variants$scored, col_types = cols(.default = col_character())) %>%
# 			type_convert()
# 			
# msk_anno = read_tsv(msk_anno_joined, col_types = cols(.default = col_character())) %>%
#   		   type_convert()
# 			
# tracker_grail = read_csv(file=patient_tracker)
# 
# tracker_impact = read_csv(file=impact_tracker)
# 
# valid_patient_ids = tracker_grail %>%
# 					filter(patient_id %in% tracker_impact$patient_id) %>%
# 					.[["patient_id"]]
#   
# indel_vars = indel_vars %>%
# 			 mutate(filter = replace(filter,
#              		patient_id == "MSK-VB-0001" &
#              		gene == "GATA3" &
#              		filter == "PASS",
#              		"CSR_MATCH_ELIMINATED"),
#          	 		ccd = replace(ccd,
#              			   		  patient_id == "MSK-VB-0001" &
#                            		  gene == "GATA3" &
#                            		  filter == "CSR_MATCH_ELIMINATED",
#                            		  0))
# 
# snv_plasma = snv_vars %>%
#   			 filter(ccd == 1,
#          			(c_panel == 1 | panel == 1),
#          			study == "TechVal",
#          			grail == 1 | MSK == 1,
#          			patient_id %in% valid_patient_ids) %>%
# 			 mutate(vtype = "SNV")
# 
# indel_plasma = indel_vars %>%
# 			   filter(ccd == 1,
#          			  (c_panel == 1 | panel == 1),
#          			  study == "TechVal",
#          			  grail == 1 | MSK == 1,
#          			  patient_id %in% valid_patient_ids) %>%
#   			   mutate(vtype = "INDEL",
#          			  altenddistmedian = as.integer(altenddistmedian))
#          			  
# healthy_snv = snv_vars %>%
#   			  filter((c_panel == 1 | panel == 1),
#          			  subj_type == "Healthy",
#          			  grail == 1) %>%
#   			  mutate(vtype = "SNV")
# 
# healthy_indel = indel_vars %>%
#   				filter((c_panel == 1 | panel == 1),
#          			    subj_type == "Healthy",
#          				grail == 1) %>%
#   				mutate(vtype = "INDEL",
#          			   altenddistmedian = as.integer(altenddistmedian))
# 
# small_vars_plasma = full_join(snv_plasma, indel_plasma) %>%
# 					full_join(healthy_snv) %>%
# 					full_join(healthy_indel)
# small_vars_plasma = small_vars_plasma %>%
#   					mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
# 
# small_vars_plasma = small_vars_plasma %>%
# 					filter(!(patient_id %in% hypermutators$patient_id)) %>%
# 					filter(!(patient_id %in% msi_hypermutators$patient_id)) %>%
# 					mutate(loc = str_c(chrom, ":", position_orig, "_", ref_orig, ">", alt_orig))
# 
# all_patient_table = small_vars_plasma %>%
# 					distinct(subj_type, patient_id)
# all_patient_table = cbind.data.frame(subj_type = rep(all_patient_table$subj_type, 4),
#                                      patient_id = rep(all_patient_table$patient_id, 4),
#                                      bio_source = rep(c("WBC_matched",
#                                                         "VUSo",
#                                                         "biopsy_matched",
#                                                         "biopsy_only"),
#                                                        each = nrow(all_patient_table)))
# 
# variants = label_bio_source(small_vars_plasma)
# variants = left_join(variants, msk_anno %>% select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
# variants = variants %>%
# 		   mutate(bio_source = case_when(
# 		   					   MSK == 1 & grail == 1 ~ "biopsy_matched",
# 		   					   MSK == 1 & grail == 0 ~ "biopsy_only",
# 		   					   category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
# 		   					   category %in% c("germline", "germlineish") ~ "germline",
# 		   					   category %in% c("blood", "bloodier") ~ "WBC_matched",
# 		   					   category == "somatic" & `IM-T.alt_count` > bam_read_count_cutoff ~ "IMPACT-BAM_matched",
# 		   					   category == "somatic" ~ "VUSo",
# 		   					   TRUE ~ "other"),
# 		   		  af = ifelse(is.na(af), 0, af),
# 		   		  af_nobaq = round(adnobaq / dpnobaq * 100, 2),
# 		   		  af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq))
#  
# variants_ctrl = variants %>%
# 				filter(is_nonsyn) %>%
# 				filter(subj_type=="Control") %>%
# 				filter(bio_source=="VUSo")
# length(unique(variants_ctrl$patient_id))
# 
# small_vars = variants %>%
# 			 filter(is_nonsyn) %>%
# 			 filter(subj_type=="Breast") %>%
# 			 filter(bio_source=="VUSo")
# length(unique(small_vars$patient_id))
# 
# small_vars = variants %>%
# 			 filter(is_nonsyn) %>%
# 			 filter(subj_type=="Lung") %>%
# 			 filter(bio_source=="VUSo")
# length(unique(small_vars$patient_id))
# 
# small_vars = variants %>%
# 			 filter(is_nonsyn) %>%
# 			 filter(subj_type=="Prostate") %>%
# 			 filter(bio_source=="VUSo")
# length(unique(small_vars$patient_id))
# 
# small_vars = variants %>%
# 			 filter(subj_type!="Control") %>%
# 			 filter(bio_source %in% c("biopsy_matched", "biopsy_only"))
# table(small_vars$bio_source)
# 
# small_vars = variants %>%
# 			 filter(is_nonsyn) %>%
# 			 filter(subj_type=="Control") %>%
# 			 filter(bio_source=="WBC_matched")
# length(unique(small_vars$patient_id))
# 
# small_vars = variants %>%
# 			 filter(is_nonsyn) %>%
# 			 filter(subj_type!="Control") %>%
# 			 filter(bio_source=="WBC_matched")
# length(unique(small_vars$patient_id))
# 
# clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
#  		   type_convert() %>%
#  		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
# tracker_grail = read_csv(file=patient_tracker)
# tracker_impact = read_csv(file=impact_tracker)
# valid_patient_ids = tracker_grail %>%
#  					filter(study=="TechVal") %>%
# 					filter(patient_id %in% tracker_impact$patient_id) %>%
# 					.[["patient_id"]]
# valid_control_ids = tracker_grail %>%
# 					filter(tissue=="Healthy" & study=="Merlin") %>%
# 					.[["patient_id"]]
# clinical = clinical %>%
#  	 	   filter(patient_id %in% valid_patient_ids | patient_id %in% valid_control_ids)
# 
# median(clinical$age[grep("VB", clinical$patient_id)])
# range(clinical$age[grep("VB", clinical$patient_id)])
# median(clinical$age[grep("VL", clinical$patient_id)])
# range(clinical$age[grep("VL", clinical$patient_id)])
# median(clinical$age[grep("VP", clinical$patient_id)])
# range(clinical$age[grep("VP", clinical$patient_id)])
# median(clinical$age[grep("W", clinical$patient_id)])
# range(clinical$age[grep("W", clinical$patient_id)])
# 
# sum(clinical$age[grep("VB", clinical$patient_id)]<=50)
# sum(clinical$age[grep("VB", clinical$patient_id)]>50 & clinical$age[grep("VB", clinical$patient_id)]<=60)
# sum(clinical$age[grep("VB", clinical$patient_id)]>60 & clinical$age[grep("VB", clinical$patient_id)]<=70)
# sum(clinical$age[grep("VB", clinical$patient_id)]>70)
# 
# sum(clinical$age[grep("VL", clinical$patient_id)]<=50)
# sum(clinical$age[grep("VL", clinical$patient_id)]>50 & clinical$age[grep("VL", clinical$patient_id)]<=60)
# sum(clinical$age[grep("VL", clinical$patient_id)]>60 & clinical$age[grep("VL", clinical$patient_id)]<=70)
# sum(clinical$age[grep("VL", clinical$patient_id)]>70)
# 
# sum(clinical$age[grep("VP", clinical$patient_id)]<=50)
# sum(clinical$age[grep("VP", clinical$patient_id)]>50 & clinical$age[grep("VP", clinical$patient_id)]<=60)
# sum(clinical$age[grep("VP", clinical$patient_id)]>60 & clinical$age[grep("VP", clinical$patient_id)]<=70)
# sum(clinical$age[grep("VP", clinical$patient_id)]>70)
# 
# sum(clinical$age[grep("W", clinical$patient_id)]<=50)
# sum(clinical$age[grep("W", clinical$patient_id)]>50 & clinical$age[grep("W", clinical$patient_id)]<=60)
# sum(clinical$age[grep("W", clinical$patient_id)]>60 & clinical$age[grep("W", clinical$patient_id)]<=70)
# sum(clinical$age[grep("W", clinical$patient_id)]>70)
# 
# table(clinical$sex[grep("VB", clinical$patient_id)])
# table(clinical$sex[grep("VL", clinical$patient_id)])
# table(clinical$sex[grep("VP", clinical$patient_id)])
# table(clinical$sex[grep("W", clinical$patient_id)])
# 
# table(clinical$receptor_status[grep("VB", clinical$patient_id)])
# table(clinical$histology[grep("VB", clinical$patient_id)])
# 
# table(clinical$histology[grep("VL", clinical$patient_id)])
# table(clinical$cancer_type[grep("VP", clinical$patient_id)])
# 
# tx = read.csv(file="../res/etc/prior_tx_techval_0818.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE, row.names=1)
# tx = subset(tx, rownames(tx) %in% clinical$patient_id)
# tx2 = tx %>%
#  	  filter(cancer_type == "breast")
# sum(tx2$prior_chemo==1)
# sum(tx2$prior_rt==1)
# tx2 = tx %>%
#  	  filter(cancer_type == "lung")
# sum(tx2$prior_chemo==1)
# sum(tx2$prior_rt==1)
# tx2 = tx %>%
#  	  filter(cancer_type == "prostate")
# sum(tx2$prior_chemo==1)
# sum(tx2$prior_rt==1)
# 
# sum(clinical$lntxcat1=="." & grepl("VB", clinical$patient_id))
# sum(clinical$lntxcat1=="=1" & grepl("VB", clinical$patient_id))
# sum(clinical$lntxcat1=="=2" & grepl("VB", clinical$patient_id))
# sum(clinical$lntxcat1==">=3" & grepl("VB", clinical$patient_id))
# 
# sum(clinical$lntxcat1=="." & grepl("VL", clinical$patient_id))
# sum(clinical$lntxcat1=="=1" & grepl("VL", clinical$patient_id))
# sum(clinical$lntxcat1=="=2" & grepl("VL", clinical$patient_id))
# sum(clinical$lntxcat1==">=3" & grepl("VL", clinical$patient_id))
# 
# sum(clinical$lntxcat1=="." & grepl("VP", clinical$patient_id))
# sum(clinical$lntxcat1=="=1" & grepl("VP", clinical$patient_id))
# sum(clinical$lntxcat1=="=2" & grepl("VP", clinical$patient_id))
# sum(clinical$lntxcat1==">=3" & grepl("VP", clinical$patient_id))
