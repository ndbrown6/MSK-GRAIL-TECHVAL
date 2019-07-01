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
# CH-derived mutations in cfDNA
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
 		   
tracker_grail = read_csv(file=patient_tracker, col_types = cols(.default = col_character()))  %>%
				type_convert()
tracker_impact = read_csv(impact_tracker, col_types = cols(.default = col_character()))  %>%
				 type_convert()
valid_patient_ids = tracker_grail %>%
  				    filter(patient_id %in% tracker_impact$patient_id | tissue == "Healthy") %>%
  				    filter(!(tissue %in% c("Breast", "Lung", "Prostate") & study=="Merlin")) %>%
  				    .[["patient_id"]]
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
burden_healthy = variants %>%
 		  		 filter(bio_source %in% c("WBC_matched", "VUSo", "biopsy_matched"), is_nonsyn) %>%
 		  		 group_by(subj_type, patient_id, bio_source) %>%
 		  		 summarize(num_called = n()) %>%
 		  		 ungroup() %>%
 		  		 left_join(clinical) %>%
 		  		 mutate(group = case_when(grepl("Control", subj_type) ~ "Control", TRUE ~ "Cancer"), group = factor(group, levels = c("Control", "Cancer"))) %>%
 		  		 filter(subj_type=="Control")
burden_cancer = variants %>%
 		  		filter(bio_source %in% c("WBC_matched", "VUSo", "biopsy_matched", "IMPACT-BAM_matched"), is_nonsyn) %>%
 		  		group_by(subj_type, patient_id, bio_source) %>%
 		  		summarize(num_called = n()) %>%
 		  		ungroup() %>%
 		  		left_join(clinical) %>%
 		  		mutate(group = case_when(grepl("Control", subj_type) ~ "Control", TRUE ~ "Cancer"), group = factor(group, levels = c("Control", "Cancer"))) %>%
 		  		filter(subj_type!="Control") %>%
 		  		filter(!patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id))

valid_patient_ids = tracker_grail %>%
  				    filter(patient_id %in% tracker_impact$patient_id | tissue == "Healthy") %>%
  				    filter(!(tissue %in% c("Breast", "Lung", "Prostate") & study=="Merlin")) %>%
  				    .[["patient_id"]]
valid_patient_ids = intersect(valid_patient_ids, clinical$patient_id)		   

burden = bind_rows(burden_cancer, burden_healthy)
burden_cfdna = data.frame(burden[,c("patient_id", "age", "num_called", "subj_type", "bio_source"),drop=FALSE]) %>%
			   filter(bio_source=="WBC_matched") %>%
	 		   filter(patient_id %in% valid_patient_ids) %>%
	 		   dplyr::select(patient_id, age, num_called, subj_type)

indx = which(!(valid_patient_ids %in% burden_cfdna$patient_id))

if (length(indx)!=0) {
	zzz = data.frame(clinical[clinical$patient_id %in% valid_patient_ids[indx],c("patient_id","age","subj_type")])
	zzz = cbind(zzz, num_called = rep(0, length(indx)))
	zzz = zzz[,colnames(burden_cfdna)]
	burden_cfdna = bind_rows(burden_cfdna, zzz)
}
		   
#==================================================
# CH-derived mutations in WBC
#==================================================
all_vars = read_tsv(file=wbc_scored_annotated_and_clinical$scored, col_types = cols(.default = col_character())) %>%
		   type_convert()
clinical = read_tsv(file=wbc_scored_annotated_and_clinical$clinical, col_types = cols(.default = col_character())) %>%
 		   type_convert()
  		   
save_vars = all_vars
 
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
 		   
burden_healthy = all_vars %>%
  		   		 filter(subj_type=="Control") %>%
  		   		 group_by(patient_id) %>%
   		   		 summarize(num_called = n()) %>%
  		   		 ungroup()
  		   		 
burden_cancer = all_vars %>%
  		   		filter(subj_type!="Control") %>%
 		   		group_by(patient_id) %>%
   	 	   		summarize(num_called = n()) %>%
   	 	   		ungroup()
   	 	   		
burden_wbc = bind_rows(burden_healthy, burden_cancer)

data = left_join(burden_cfdna, burden_wbc, by="patient_id") %>%
 	   dplyr::rename(num_called_cfdna = num_called.x, num_called_wbc = num_called.y) %>%
 	   mutate(subj_type = ifelse(subj_type=="Control", "Healthy", "Cancer")) %>%
 	   mutate(subj_type = ifelse(grepl("W", patient_id), "Healthy", "Cancer")) %>%
 	   mutate(num_called_cfdna = ifelse(is.na(num_called_cfdna), 0, num_called_cfdna)) %>%
 	   mutate(num_called_wbc = ifelse(is.na(num_called_wbc), 0, num_called_wbc))

cfdna_fraction = read_csv(file=url_ctdna_frac, col_types = cols(.default = col_character())) %>%
 		   		 type_convert() %>%
 		   		 mutate(ctdna_frac = ifelse(is.na(ctdna_frac), 0, ctdna_frac)) %>%
 		   		 dplyr::rename(patient_id = ID)
 		   		 
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

qc_metrics_cfdna = read.csv(file=url_qc_metrics_cfdna, header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
			 	   dplyr::select(sample_id, patient_id, sample_type, tissue, volume_of_blood_mL, volume_of_DNA_source_mL, DNA_extraction_yield_ng, DNA_input_concentration_ng_uL, Library_preparation_input_ng, raw.MEAN_BAIT_COVERAGE, collapsed.MEAN_BAIT_COVERAGE, collapsed_fragment.MEAN_BAIT_COVERAGE, readErrorRate, readSubstErrorRate, Study) %>%
			 	   filter(sample_type=="cfDNA")
 
data = left_join(data, cfdna_fraction, by="patient_id") %>%
	   left_join(qc_metrics_cfdna %>% dplyr::select(patient_id, uncollapsed_cfdna_coverage=raw.MEAN_BAIT_COVERAGE, collapsed_cfdna_coverage=collapsed.MEAN_BAIT_COVERAGE), by="patient_id") %>%
 	   mutate(ctdna_frac = ifelse(is.na(ctdna_frac), 0, ctdna_frac)) %>%
 	   mutate(bio_source = "Somatic CH-derived mutations / Mb")
 
fit.null = lm(num_called_cfdna ~ num_called_wbc, data=data %>% filter(subj_type=="Cancer"))
fit.wage = lm(num_called_cfdna ~ num_called_wbc + age, data=data %>% filter(subj_type=="Cancer"))
lrt = lrtest(fit.wage, fit.null)

fit.wucov = lm(num_called_cfdna ~ num_called_wbc + uncollapsed_cfdna_coverage, data=data %>% filter(subj_type=="Cancer"))
lrt = lrtest(fit.wucov, fit.null)

fit.wccov = lm(num_called_cfdna ~ num_called_wbc + collapsed_cfdna_coverage, data=data %>% filter(subj_type=="Cancer"))
lrt = lrtest(fit.wccov, fit.null)

fit.wctdna = lm(num_called_cfdna ~ num_called_wbc + ctdna_frac, data=data %>% filter(subj_type=="Cancer"))
lrt = lrtest(fit.wctdna, fit.null)

fit.alt = lm(num_called_cfdna ~ num_called_wbc + age + uncollapsed_cfdna_coverage + collapsed_cfdna_coverage + ctdna_frac, data=data %>% filter(subj_type=="Cancer"))

plot.0 = ggplot(data %>% filter(subj_type=="Healthy"), aes(x = num_called_wbc, y = num_called_cfdna)) +
		 geom_point(alpha=.85, shape = 24, size=2.5, color = "#231F20", fill = "#FDAE61") +
		 geom_point(alpha=.85, shape = 21, color = "#231F20", fill = "salmon", aes(x = num_called_wbc, y = num_called_cfdna, fill = subj_type, size=ctdna_frac), data = data %>% filter(subj_type=="Cancer"), inherit.aes=FALSE) +
		 geom_abline(slope = 1, intercept = 0, color = "goldenrod3", size = .75) +
		 geom_smooth(formula = y ~ x +0, method="lm", aes(x = num_called_wbc, y = num_called_cfdna), data = data %>% filter(subj_type=="Cancer") , inherit.aes=FALSE, color="grey50", fill="grey50", fullrange=TRUE) +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\nWBC\n", y="cfDNA\n") +
		 coord_cartesian(xlim = c(0, 40), ylim = c(0, 40)) +
		 guides(size=guide_legend(title=c("ctDNA fraction")))
		 
pdf(file="../res/rebuttal/cfDNA_matched_vs_WBC_Combined.pdf", width=5.5, height=6)
print(plot.0)
dev.off()

plot.0 = ggplot(data %>% filter(subj_type=="Cancer"), aes(x = num_called_wbc, y = num_called_cfdna, size=ctdna_frac)) +
		 geom_point(alpha=.85, shape = 21, color = "#231F20", fill="salmon") +
		 geom_abline(slope = 1, intercept = 0, color = "goldenrod3", size = .75) +
		 geom_smooth(formula = y ~ x +0, method="lm", aes(x = num_called_wbc, y = num_called_cfdna), data = data %>% filter(subj_type=="Cancer"), inherit.aes=FALSE, color="grey50", fill="grey50", fullrange=TRUE) +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\nWBC\n", y="cfDNA\n") +
		 coord_cartesian(xlim = c(0, 40), ylim = c(0, 40)) +
		 guides(size=guide_legend(title=c("ctDNA fraction")))
	  	
pdf(file="../res/rebuttal/cfDNA_matched_vs_WBC_Cancer.pdf", width=5.5, height=6)
print(plot.0)
dev.off()

plot.0 = ggplot(data %>% filter(subj_type=="Healthy"), aes(x = num_called_wbc, y = num_called_cfdna)) +
		 geom_point(alpha=.85, shape = 21, color = "#231F20", fill="#FDAE61", size=2.5) +
		 geom_abline(slope = 1, intercept = 0, color = "goldenrod3", size = .75) +
		 geom_smooth(formula = y ~ x +0, method="lm", aes(x = num_called_wbc, y = num_called_cfdna), data = data %>% filter(subj_type=="Healthy") %>% filter(num_called_wbc<20), inherit.aes=FALSE, color="grey50", fill="grey50", fullrange=TRUE) +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\nWBC\n", y="cfDNA\n") +
		 coord_cartesian(xlim = c(0, 40), ylim = c(0, 40)) +
		 guides(size=guide_legend(title=c("ctDNA fraction")))
	  	
pdf(file="../res/rebuttal/cfDNA_matched_vs_WBC_Healthy.pdf", width=5.5, height=6)
print(plot.0)
dev.off()

'gatherpairs' <- function (data, ..., xkey = '.xkey', xvalue = '.xvalue', ykey = '.ykey', yvalue = '.yvalue', na.rm = FALSE, convert = FALSE, factor_key = FALSE)
{
  vars <- quos(...)
  xkey <- enquo(xkey)
  xvalue <- enquo(xvalue)
  ykey <- enquo(ykey)
  yvalue <- enquo(yvalue)

  data %>% {
    cbind(gather(., key = !!xkey, value = !!xvalue, !!!vars,
                 na.rm = na.rm, convert = convert, factor_key = factor_key),
          dplyr::select(., !!!vars)) 
  } %>% gather(., key = !!ykey, value = !!yvalue, !!!vars,
               na.rm = na.rm, convert = convert, factor_key = factor_key)
}
