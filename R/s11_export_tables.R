#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
source('config.R')

if (!dir.exists("../res/tables")) {
	dir.create("../res/tables")
}

#==================================================
# export vcf tables for vcf2maf
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
 		   		  
cfdna_variants = variants %>%
				 filter(is_nonsyn | bio_source=="biopsy_only") %>%
 				 filter(bio_source %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo"))
cfdna_variants = as.data.frame(cfdna_variants[,c("chrom","position","loc","ref_orig","alt_orig")])
cfdna_variants = cbind(cfdna_variants,
					   rep(".", nrow(cfdna_variants)),
					   rep("PASS", nrow(cfdna_variants)),
					   rep(".", nrow(cfdna_variants)))
cfdna_variants[,1] = gsub(pattern="chr", replacement="", x=cfdna_variants[,1], fixed=TRUE)
colnames(cfdna_variants) = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
write.table(cfdna_variants, file="../res/etc/export-cfdna-variants.vcf", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

impact_variants = variants %>%
				  filter(is_nonsyn | bio_source=="biopsy_only") %>%
 				  filter(bio_source %in% c("biopsy_matched", "biopsy_only", "IMPACT-BAM_matched"))
index = impact_variants$bio_source=="biopsy_only"
impact_variants$ref_orig[index] = impact_variants$c_ref_orig[index]
impact_variants$alt_orig[index] = impact_variants$c_alt_orig[index]
impact_variants = as.data.frame(impact_variants[,c("chrom","position","loc","ref_orig","alt_orig")])
impact_variants = cbind(impact_variants,
					   rep(".", nrow(impact_variants)),
					   rep("PASS", nrow(impact_variants)),
					   rep(".", nrow(impact_variants)))
impact_variants[,1] = gsub(pattern="chr", replacement="", x=impact_variants[,1], fixed=TRUE)
colnames(impact_variants) = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
index = c(786, 811, 819, 835, 856, 874, 880)-1
impact_variants = impact_variants[-index,,drop=FALSE]
write.table(impact_variants, file="../res/etc/export-tumor-variants.vcf", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

#==================================================
# Table S7: cfDNA variants
#==================================================
cfdna_variants = variants %>%
				 filter(is_nonsyn | bio_source=="biopsy_only") %>%
 				 filter(bio_source %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo"))
 				 
cfdna_variants_export = read.csv(file="../res/etc/export-cfdna-variants.maf", header=TRUE, sep="\t", comment.char="#", stringsAsFactors=FALSE) %>%
							 mutate(Tumor_Sample_Barcode = cfdna_variants$patient_id) %>%
							 mutate(Matched_Norm_Sample_Barcode = cfdna_variants$patient_id) %>%
							 mutate(t_depth = cfdna_variants$dpnobaq) %>%
							 mutate(t_ref_count = cfdna_variants$dpnobaq - cfdna_variants$adnobaq) %>%
							 mutate(t_alt_count = cfdna_variants$adnobaq) %>%
							 mutate(PHENO = cfdna_variants$bio_source) %>%
							 mutate(Hugo_Symbol = cfdna_variants$gene) %>%
							 mutate(SYMBOL = cfdna_variants$gene) %>%
							 mutate(Gene = ".") %>%
							 mutate(Entrez_Gene_Id = ".") %>%
							 mutate(Center = "GRAIL")
write.table(cfdna_variants_export, file="../res/tables/Table_S7.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

#==================================================
# Table S9: Tumor variants
#==================================================
impact_variants = variants %>%
				  filter(is_nonsyn | bio_source=="biopsy_only") %>%
 				  filter(bio_source %in% c("biopsy_matched", "biopsy_only", "IMPACT-BAM_matched"))
index = c(786, 811, 819, 835, 856, 874, 880)-1
impact_variants = impact_variants[-index,,drop=FALSE]
impact_variants_export = read.csv(file="../res/etc/export-tumor-variants.maf", header=TRUE, sep="\t", comment.char="#", stringsAsFactors=FALSE) %>%
							  mutate(Tumor_Sample_Barcode = impact_variants$patient_id) %>%
							  mutate(Matched_Norm_Sample_Barcode = impact_variants$patient_id) %>%
							  mutate(t_depth = impact_variants$`IM-T.ref_count` + impact_variants$`IM-T.alt_count`) %>%
							  mutate(t_ref_count = impact_variants$`IM-T.ref_count`) %>%
							  mutate(t_alt_count = impact_variants$`IM-T.alt_count`) %>%
							  mutate(PHENO = impact_variants$bio_source) %>%
							  mutate(Hugo_Symbol = impact_variants$gene) %>%
							  mutate(SYMBOL = impact_variants$gene) %>%
 							  mutate(Gene = ".") %>%
							  mutate(Entrez_Gene_Id = ".") %>%
							  mutate(Center = "MSKCC")
tmp = read.csv(file="../res/etc/export-tumor-variants-submission-na-filled-manual-msk-vb-0044 (1).txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
rownames(impact_variants_export) = paste0(impact_variants_export[,"Tumor_Sample_Barcode"], ":", impact_variants_export[,"Chromosome"], ":", impact_variants_export[,"Start_Position"], ":", impact_variants_export[,"End_Position"], ":", impact_variants_export[,"Tumor_Seq_Allele1"], ">", impact_variants_export[,"Tumor_Seq_Allele2"])
rownames(tmp) = paste0(tmp[,"Tumor_Sample_Barcode"], ":", tmp[,"Chromosome"], ":", tmp[,"Start_Position"], ":", tmp[,"End_Position"], ":", tmp[,"Tumor_Seq_Allele1"], ">", tmp[,"Tumor_Seq_Allele2"])
tmp = tmp[rownames(impact_variants_export),,drop=FALSE]
impact_variants_export[,"t_depth"] = tmp[,"t_depth"]
impact_variants_export[,"t_ref_count"] = tmp[,"t_ref_count"]
impact_variants_export[,"t_alt_count"] = tmp[,"t_alt_count"]
write.table(impact_variants_export, file="../res/tables/Table_S9.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

#==================================================
# Table S8: WBC variants
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
 		   af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq)) %>%
 		   mutate(ID_x = make_id_x(patient_id, chrom, position, ref, alt))
 		   
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
			   
som_vars_breast = read_tsv(file=somatic_vars_breast) %>%
				  type_convert()
som_vars_lung = read_tsv(file=somatic_vars_lung) %>%
				type_convert()
som_vars_prostate = read_tsv(file=somatic_vars_prostate) %>%
					type_convert()
som_vars = bind_rows(som_vars_breast, som_vars_lung, som_vars_prostate) %>%
		   mutate(ID_x = make_id_x(CASE, CHROM, POS, REF, ALT)) %>%
		   mutate(ID_y = make_id_y(CASE, vcf2maf_HGVSp_Short))
load(cosmic_file)
cosmic_db = cosmic_db %>%
			filter(SNP==0) %>%
			filter(as.numeric(n)>=2) %>%
			mutate(id = paste0("chr", CHROM, ":", POS, "_", REF, ">", ALT))
load(gnomad_file)
gnomad_db = gnomad_db %>%
			filter(AF>=0.01) %>%
			mutate(id = paste0("chr", CHROM, ":", POS, "_", REF, ">", ALT))
load(hotspot_file)
cancer_hotspot = cancer_hotspot %>%
				 mutate(id = paste0("chr", Chromosome, ":", Start_Position, "_", Reference_Allele, ">", Tumor_Seq_Allele2))
	 		   
all_vars = read_tsv(wbc_variants$scored, col_types = cols(.default = col_character()))  %>%
		   type_convert()
all_annotations = read_tsv(wbc_variants$annotations, comment="#", col_types = cols(.default = col_character()))  %>%
				  type_convert()
all_vars = bind_cols(all_vars, all_annotations) %>%
		   left_join(clinical, by=c("patient_id","study")) %>%
           mutate(subj_type = ifelse(subj_type=="Healthy" | is.na(subj_type), "Control", subj_type)) %>%
 		   mutate(is_patient_valid = (patient_id %in% valid_patient_ids)) %>%
 		   mutate(ID_x = make_id_x(patient_id, chrom, pos, ref, alt)) %>%
 		   mutate(ID_y = make_id_y(patient_id, HGVSp_Short)) %>%
 		   mutate(loc_lng = str_c(chrom, ":", pos, "_", ref, ">", alt)) %>%
 		   mutate(loc_srt = str_c(chrom, ":", pos)) %>%
 		   mutate(c_panel = in_common_bed(chrom, pos)) %>%
 		   mutate(is_hypermutator = (patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id))) %>%
 		   mutate(is_lowdepth = (dpnobaq < 500)) %>%
 		   mutate(is_lowqual = (qualnobaq < 60)) %>%
 		   mutate(is_tumor_matched = (ID_x %in% som_vars$ID_x) | (ID_y %in% som_vars$ID_y)) %>%
 		   mutate(is_cfdna_matched = (ID_x %in% variants$ID_x[variants$bio_source %in% c("biopsy_matched", "biopsy_only", "germline", "IMPACT-BAM_matched")])) %>%
 		   mutate(in_cosmic = loc_lng %in% cosmic_db$id) %>%
  		   mutate(in_exac = !is.na(ExAC_AF) & ExAC_AF>=0.01) %>%
  		   mutate(in_gnomad = loc_lng %in% gnomad_db$id) %>%
  		   mutate(is_hotspot = loc_lng %in% cancer_hotspot$id)

all_vars = all_vars %>%
 		   filter(is_patient_valid) %>%
 		   filter(c_panel) %>%
 		   filter(!is_hypermutator) %>%
 		   filter(!is_lowdepth) %>%
 		   filter(!is_lowqual) %>%
 		   filter(!is_tumor_matched) %>%
 		   filter(!is_cfdna_matched)
 
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
  		   
recurrence = all_vars %>%
   			 group_by(patient_id, loc_srt) %>%
   			 count() %>%
   			 ungroup() %>%
   			 rename(n_indel=n)
all_vars = all_vars %>%
		   left_join(recurrence)
all_vars = all_vars %>%
		   filter(!(n_indel > 1 & indel))
		   
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
  		   
all_vars = all_vars %>%
 		   filter(SYMBOL!="HLA-A")

all_vars = all_vars %>%
 		   filter((adnobaq/dpnobaq)<=.3 | (Variant_Classification %in% c("Frame_Shift_Ins", "Frame_Shift_Del", "Nonsense_Mutation", "Nonstop_Mutation") & SYMBOL %in% chip_genes))
		   
all_vars = all_vars %>%
 		   filter(!in_exac)
 		   
all_vars = all_vars %>%
 		   filter(!in_gnomad)

all_vars = all_vars %>%
		   mutate(Tumor_Sample_Barcode = all_vars$patient_id) %>%
		   mutate(Matched_Norm_Sample_Barcode = all_vars$patient_id) %>%
		   mutate(n_depth = all_vars$dpnobaq) %>%
		   mutate(n_ref_count = all_vars$dpnobaq - all_vars$adnobaq) %>%
		   mutate(n_alt_count = all_vars$adnobaq) %>%
		   mutate(PHENO = "CH_derived") %>%
		   mutate(Gene = ".") %>%
		   mutate(Entrez_Gene_Id = ".") %>%
		   mutate(Center = "GRAIL")

all_vars = all_vars[,colnames(impact_variants_export),drop=FALSE]
write.table(all_vars, file="../res/tables/Table_S8.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

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
oncotree_code = rep("", length(valid_patient_ids))
oncotree_code[grep("VB", valid_patient_ids, fixed=TRUE)] = "BRCA"
oncotree_code[grep("VL", valid_patient_ids, fixed=TRUE)] = "NSCLC"
oncotree_code[grep("VP", valid_patient_ids, fixed=TRUE)] = "PRAD"
oncotree = cbind(SAMPLE_ID=valid_patient_ids,
				 ONCOTREE_CODE=oncotree_code)
write.table(oncotree, file="../res/tables/clinical.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

