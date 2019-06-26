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
# Pie-charts of the number of variants classified
# by source across control, cancer non-hypermutators
# and cancer hypermutators
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
total_bed_Mb = sum(GenomicRanges::width(bed_ranges)) / 1e6

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
 		   
guardant_g360 = read_tsv(url_guardant_g360, col_types = cols(.default = col_character())) %>%
				type_convert()

variants_control = variants %>%
		   filter(gene %in% guardant_g360$gene_id) %>%
		   filter(subj_type=="Control") %>%
		   filter(is_nonsyn | bio_source=="biopsy_only") %>%
		   filter(bio_source %in% c("VUSo","biopsy_matched", "biopsy_only", "IMPACT-BAM_matched", "WBC_matched"))

variants_nohyper = variants %>%
		   filter(gene %in% guardant_g360$gene_id) %>%
		   filter(subj_type!="Control") %>%
		   filter(!patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id)) %>%
		   filter(is_nonsyn | bio_source=="biopsy_only") %>%
		   filter(bio_source %in% c("VUSo","biopsy_matched", "biopsy_only", "IMPACT-BAM_matched", "WBC_matched"))
 
variants_hyper = variants %>%
		   filter(gene %in% guardant_g360$gene_id) %>%
		   filter(subj_type!="Control") %>%
		   filter(patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id)) %>%
		   filter(is_nonsyn | bio_source=="biopsy_only") %>%
		   filter(bio_source %in% c("VUSo","biopsy_matched", "biopsy_only", "IMPACT-BAM_matched", "WBC_matched"))
 		   
M2 = matrix(NA, nrow=5, ncol=5)
M2[1,1] = sum(variants_control$bio_source=="biopsy_only")
M2[2,1] = sum(variants_nohyper$bio_source=="biopsy_only")
M2[3,1] = sum(variants_nohyper$bio_source=="biopsy_only" & variants_nohyper$subj_type=="Breast")
M2[4,1] = sum(variants_nohyper$bio_source=="biopsy_only" & variants_nohyper$subj_type=="Lung")
M2[5,1] = sum(variants_nohyper$bio_source=="biopsy_only" & variants_nohyper$subj_type=="Prostate")
M2[1,2] = sum(variants_control$bio_source=="biopsy_matched")
M2[2,2] = sum(variants_nohyper$bio_source=="biopsy_matched")
M2[3,2] = sum(variants_nohyper$bio_source=="biopsy_matched" & variants_nohyper$subj_type=="Breast")
M2[4,2] = sum(variants_nohyper$bio_source=="biopsy_matched" & variants_nohyper$subj_type=="Lung")
M2[5,2] = sum(variants_nohyper$bio_source=="biopsy_matched" & variants_nohyper$subj_type=="Prostate")
M2[1,3] = sum(variants_control$bio_source=="IMPACT-BAM_matched")
M2[2,3] = sum(variants_nohyper$bio_source=="IMPACT-BAM_matched")
M2[3,3] = sum(variants_nohyper$bio_source=="IMPACT-BAM_matched" & variants_nohyper$subj_type=="Breast")
M2[4,3] = sum(variants_nohyper$bio_source=="IMPACT-BAM_matched" & variants_nohyper$subj_type=="Lung")
M2[5,3] = sum(variants_nohyper$bio_source=="IMPACT-BAM_matched" & variants_nohyper$subj_type=="Prostate")
M2[1,4] = sum(variants_control$bio_source=="WBC_matched")
M2[2,4] = sum(variants_nohyper$bio_source=="WBC_matched")
M2[3,4] = sum(variants_nohyper$bio_source=="WBC_matched" & variants_nohyper$subj_type=="Breast")
M2[4,4] = sum(variants_nohyper$bio_source=="WBC_matched" & variants_nohyper$subj_type=="Lung")
M2[5,4] = sum(variants_nohyper$bio_source=="WBC_matched" & variants_nohyper$subj_type=="Prostate")
M2[1,5] = sum(variants_control$bio_source=="VUSo")
M2[2,5] = sum(variants_nohyper$bio_source=="VUSo")
M2[3,5] = sum(variants_nohyper$bio_source=="VUSo" & variants_nohyper$subj_type=="Breast")
M2[4,5] = sum(variants_nohyper$bio_source=="VUSo" & variants_nohyper$subj_type=="Lung")
M2[5,5] = sum(variants_nohyper$bio_source=="VUSo" & variants_nohyper$subj_type=="Prostate")
 
M3 = matrix(NA, nrow=5, ncol=5)
M3[1,1] = 0
M3[2,1] = sum(variants_hyper$bio_source=="biopsy_only")
M3[3,1] = sum(variants_hyper$bio_source=="biopsy_only" & variants_hyper$subj_type=="Breast")
M3[4,1] = sum(variants_hyper$bio_source=="biopsy_only" & variants_hyper$subj_type=="Lung")
M3[5,1] = sum(variants_hyper$bio_source=="biopsy_only" & variants_hyper$subj_type=="Prostate")
M3[1,2] = 0
M3[2,2] = sum(variants_hyper$bio_source=="biopsy_matched")
M3[3,2] = sum(variants_hyper$bio_source=="biopsy_matched" & variants_hyper$subj_type=="Breast")
M3[4,2] = sum(variants_hyper$bio_source=="biopsy_matched" & variants_hyper$subj_type=="Lung")
M3[5,2] = sum(variants_hyper$bio_source=="biopsy_matched" & variants_hyper$subj_type=="Prostate")
M3[1,3] = 0
M3[2,3] = sum(variants_hyper$bio_source=="IMPACT-BAM_matched")
M3[3,3] = sum(variants_hyper$bio_source=="IMPACT-BAM_matched" & variants_hyper$subj_type=="Breast")
M3[4,3] = sum(variants_hyper$bio_source=="IMPACT-BAM_matched" & variants_hyper$subj_type=="Lung")
M3[5,3] = sum(variants_hyper$bio_source=="IMPACT-BAM_matched" & variants_hyper$subj_type=="Prostate")
M3[1,4] = 0
M3[2,4] = sum(variants_hyper$bio_source=="WBC_matched")
M3[3,4] = sum(variants_hyper$bio_source=="WBC_matched" & variants_hyper$subj_type=="Breast")
M3[4,4] = sum(variants_hyper$bio_source=="WBC_matched" & variants_hyper$subj_type=="Lung")
M3[5,4] = sum(variants_hyper$bio_source=="WBC_matched" & variants_hyper$subj_type=="Prostate")
M3[1,5] = 0
M3[2,5] = sum(variants_hyper$bio_source=="VUSo")
M3[3,5] = sum(variants_hyper$bio_source=="VUSo" & variants_hyper$subj_type=="Breast")
M3[4,5] = sum(variants_hyper$bio_source=="VUSo" & variants_hyper$subj_type=="Lung")
M3[5,5] = sum(variants_hyper$bio_source=="VUSo" & variants_hyper$subj_type=="Prostate")
 
colnames(M2) = colnames(M3) = c("biopsy_only", "biopsy_matched", "biopsy_subthreshold", "WBC_mathed", "VUSo")
rownames(M2) = rownames(M3) = c("Control", "Cancer", "Breast", "Lung", "Prostate")
 
pdf(file="../res/rebuttal/pie_bio_sources.pdf", height=7, width=15)
par(mar = c(6.1, 6.5, 4.1, 1.1))
z = split.screen(figs=matrix(c(0, 0.33, 0, 1, .33, .67, 0, 1, .67, 1, 0, 1), nrow=3, ncol=4, byrow=TRUE))
screen(z[1])
pie(x=M2[1,4:5], col=c("#EA7180", "#4E9B3C"), labels=M2[1,4:5])
screen(z[2])
pie(x=M2[2,2:5], col=c("#2C80C3", "#F2B151", "#EA7180", "#4E9B3C"), labels=M2[2,2:5])
screen(z[3])
pie(x=M3[2,2:5], col=c("#2C80C3", "#F2B151", "#EA7180", "#4E9B3C"), labels=M3[2,2:5])
close.screen(all.screens=TRUE)
dev.off()
