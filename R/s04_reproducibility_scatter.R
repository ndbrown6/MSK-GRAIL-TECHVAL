#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/figureS4")) {
	dir.create("../res/figureS4")
}

#==================================================
# 2-by-2 scatterplots of technical replicates
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

patient_ids = c("MSK-VB-0050", "MSK-VB-0041", "MSK-VL-0028", "MSK-VL-0042", "MSK-VB-0023", "MSK-VL-0038")
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
				  dplyr::select(-t1, -t2, -t3, -t4, -g1, -g2, -g3, -g4, -p1, -p2, -p3, -p4) %>%
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
				  dplyr::select(-min_p) %>%
				  mutate(loc = str_c(chrom, ":", pos, "_", ref, ">", alt)) %>%
				  mutate(position_orig = pos, ref_orig = ref, alt_orig = alt, position = pos) %>%
				  mutate(gratio = (adnobaq+2)*(dpgdna+4)/((adgdna+2)*(dpnobaq+4)),
				  		 gzero = adgdna/sqrt(adgdna+2)) %>%
				  mutate(start = position_orig,
                         end = position_orig + 1) %>%
                  genome_left_join(clean_target_region, by = c("chrom", "start", "end")) %>%
                  mutate(chrom = chrom.x) %>%
                  dplyr::select(-c(chrom.x, chrom.y, start.x, end.x, start.y, end.y)) %>%
				  left_join(variants %>% dplyr::select(patient_id, chrom, position_orig, ref_orig, alt_orig, MSK, grail), by=c("patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
				  mutate(MSK = ifelse(is.na(MSK), 0, MSK)) %>%
				  mutate(grail = 1) %>%
				  mutate(study = "TechVal") %>%
				  filter(in_target)
				  
feature_names = intersect(colnames(small_vars_plasma), colnames(techval_repeats))
repeat_variants = bind_rows(small_vars_plasma[,feature_names,drop=FALSE] %>%
							filter(!(patient_id %in% patient_ids)),
							techval_repeats[,feature_names,drop=FALSE])

repeat_variants = label_bio_source(repeat_variants)
repeat_variants = left_join(repeat_variants, msk_anno %>% dplyr::select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
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
			
patient_ids = c("MSK-VL-0028", "MSK-VL-0042", "MSK-VB-0023")
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
				  dplyr::select(-t1, -t2, -t3, -t4, -g1, -g2, -g3, -g4, -p1, -p2, -p3, -p4) %>%
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
				  dplyr::select(-min_p) %>%
				  mutate(loc = str_c(chrom, ":", pos, "_", ref, ">", alt)) %>%
				  mutate(position_orig = pos, ref_orig = ref, alt_orig = alt, position = pos) %>%
				  mutate(gratio = (adnobaq+2)*(dpgdna+4)/((adgdna+2)*(dpnobaq+4)),
				  		 gzero = adgdna/sqrt(adgdna+2)) %>%
				  mutate(start = position_orig,
                         end = position_orig + 1) %>%
                  genome_left_join(clean_target_region, by = c("chrom", "start", "end")) %>%
                  mutate(chrom = chrom.x) %>%
                  dplyr::select(-c(chrom.x, chrom.y, start.x, end.x, start.y, end.y)) %>%
				  left_join(variants %>% dplyr::select(patient_id, chrom, position_orig, ref_orig, alt_orig, MSK, grail), by=c("patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
				  mutate(MSK = ifelse(is.na(MSK), 0, MSK)) %>%
				  mutate(grail = 1) %>%
				  mutate(study = "TechVal") %>%
				  filter(in_target)
				  
feature_names = intersect(colnames(small_vars_plasma), colnames(techval_repeats))
repeat_variants = bind_rows(small_vars_plasma[,feature_names,drop=FALSE] %>%
							filter(!(patient_id %in% patient_ids)),
							techval_repeats[,feature_names,drop=FALSE])

repeat_variants = label_bio_source(repeat_variants)
repeat_variants = left_join(repeat_variants, msk_anno %>% dplyr::select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
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
			

patient_ids = c("MSK-VB-0050", "MSK-VB-0041", "MSK-VL-0028", "MSK-VL-0042", "MSK-VB-0023", "MSK-VL-0038")			
all_vars = full_join(vars_rep0 %>% mutate(replicate = 1),
					 vars_rep1 %>% mutate(replicate = 2), by=c("study", "subj_type", "patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
		   dplyr::select(patient_id, adnobaq.x, adnobaq.y, dpnobaq.x, dpnobaq.y, bio_source.x, bio_source.y) %>%
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


cols = c("Not detected in one replicate"="#D7191C",
		 "Not called in one replicate due\nto low quality"="#FDAE61",
		 "Incorrect assignment between replicates"="#ABDDA4",
		 "Called in both replicates"="#2B83BA")
		 
for (i in 1:length(patient_ids)) {
	print(patient_ids[i])
	tmp_vars = all_vars %>% filter(patient_id == patient_ids[i])
	
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
			   mutate(fill = ifelse(bio_source.x=="VUSo" & bio_source.y == "WBC_matched", "Incorrect assignment between replicates", fill)) %>%
			   mutate(fill = ifelse(bio_source.y=="VUSo" & bio_source.x == "WBC_matched", "Incorrect assignment between replicates", fill))
	
	pdf(file=paste0("../res/figureS4/", patient_ids[i], "_R1_R2.pdf"), width=6.5, height=7)
	par(mar = c(6.1, 6, 4.1, 1))
	epsilon = 0
	shapes = c("Biopsy matched"=24,
			   "Biopsy unmatched"=21)
	cols = c("Not detected in one replicate"="#D7191C",
			 "Not called in one replicate due\nto low quality"="#FDAE61",
			 "Incorrect assignment between replicates"="#ABDDA4",
			 "Called in both replicates"="#2B83BA")
	plot(1, 1, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, xlim=c(0.01, 100), ylim=c(0.01, 100), log="xy")
	axis(1, at=c(0.01, 0.1, 1, 10, 100), labels=c("0", "0.1", "1", "10", "100"), cex.axis = 1.75, padj = 0.25, lwd=1.85, lwd.ticks=1.75)
	axis(2, at=c(0.01, 0.1, 1, 10, 100), labels=c("0", "0.1", "1", "10", "100"), cex.axis = 1.75, las = 1, lwd=1.85, lwd.ticks=1.75)
	mtext(side = 1, text = "Replicate 1 (%)", line = 4, cex = 1.75)
	mtext(side = 2, text = "Replicate 2 (%)", line = 4, cex = 1.75)
	points(c(.01,100), c(.01,100), type="l", lty=1, lwd=2, col="goldenrod3")
	x = tmp_vars$afnobaq.x
	y = tmp_vars$afnobaq.y
	z1 = as.character(tmp_vars$shape)
	z2 = as.character(tmp_vars$fill)
	points(x, y, pch=shapes[z1], col="black", bg=cols[z2], cex=1.65)
	log10_axis(side=1, at=c(0.01, 0.1, 1, 10, 100), lwd=0, lwd.ticks=1)
	log10_axis(side=2, at=c(0.01, 0.1, 1, 10, 100), lwd=0, lwd.ticks=1)
	dev.off()
	
}

patient_ids = c("MSK-VL-0028", "MSK-VL-0042", "MSK-VB-0023")			
all_vars = full_join(vars_rep0 %>% mutate(replicate = 1),
					 vars_rep2 %>% mutate(replicate = 3), by=c("study", "subj_type", "patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
		   dplyr::select(patient_id, adnobaq.x, adnobaq.y, dpnobaq.x, dpnobaq.y, bio_source.x, bio_source.y) %>%
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


for (i in 1:length(patient_ids)) {
	print(patient_ids[i])
	tmp_vars = all_vars %>% filter(patient_id == patient_ids[i])
	
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
			   mutate(fill = ifelse(bio_source.x=="VUSo" & bio_source.y == "WBC_matched", "Incorrect assignment between replicates", fill)) %>%
			   mutate(fill = ifelse(bio_source.y=="VUSo" & bio_source.x == "WBC_matched", "Incorrect assignment between replicates", fill))
			   	
	pdf(file=paste0("../res/figureS4/", patient_ids[i], "_R1_R3.pdf"), width=6.5, height=7)
	par(mar = c(6.1, 6, 4.1, 1))
	epsilon = 0
	shapes = c("Biopsy matched"=24,
			   "Biopsy unmatched"=21)
	cols = c("Not detected in one replicate"="#D7191C",
			 "Not called in one replicate due\nto low quality"="#FDAE61",
			 "Incorrect assignment between replicates"="#ABDDA4",
			 "Called in both replicates"="#2B83BA")
	plot(1, 1, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, xlim=c(0.01, 100), ylim=c(0.01, 100), log="xy")
	axis(1, at=c(0.01, 0.1, 1, 10, 100), labels=c("0", "0.1", "1", "10", "100"), cex.axis = 1.75, padj = 0.25, lwd=1.85, lwd.ticks=1.75)
	axis(2, at=c(0.01, 0.1, 1, 10, 100), labels=c("0", "0.1", "1", "10", "100"), cex.axis = 1.75, las = 1, lwd=1.85, lwd.ticks=1.75)
	mtext(side = 1, text = "Replicate 1 (%)", line = 4, cex = 1.75)
	mtext(side = 2, text = "Replicate 3 (%)", line = 4, cex = 1.75)
	points(c(.01,100), c(.01,100), type="l", lty=1, lwd=2, col="goldenrod3")
	x = tmp_vars$afnobaq.x
	y = tmp_vars$afnobaq.y
	z1 = as.character(tmp_vars$shape)
	z2 = as.character(tmp_vars$fill)
	points(x, y, pch=shapes[z1], col="black", bg=cols[z2], cex=1.65)
	log10_axis(side=1, at=c(0.01, 0.1, 1, 10, 100), lwd=0, lwd.ticks=1)
	log10_axis(side=2, at=c(0.01, 0.1, 1, 10, 100), lwd=0, lwd.ticks=1)
	dev.off()
	
}

patient_ids = c("MSK-VL-0028", "MSK-VL-0042", "MSK-VB-0023")			
all_vars = full_join(vars_rep1 %>% mutate(replicate = 2),
					 vars_rep2 %>% mutate(replicate = 3), by=c("study", "subj_type", "patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
		   dplyr::select(patient_id, adnobaq.x, adnobaq.y, dpnobaq.x, dpnobaq.y, bio_source.x, bio_source.y) %>%
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


for (i in 1:length(patient_ids)) {
	print(patient_ids[i])
	tmp_vars = all_vars %>% filter(patient_id == patient_ids[i])
	
	# fix biopsy_matched in y
	index = tmp_vars$bio_source.x == "biopsy_matched" & (tmp_vars$bio_source.y !="biopsy_matched" & tmp_vars$bio_source.y != "unmatched")
	if (sum(index)!=0) {
		tmp_vars$bio_source.y[index] = "biopsy_matched"
	}
	
	# fix IMPACT-BAM_matched in y
	index = tmp_vars$bio_source.x == "IMPACT-BAM_matched" & (tmp_vars$bio_source.y !="IMPACT-BAM_matched" & tmp_vars$bio_source.y != "unmatched")
	if (sum(index)!=0) {
		tmp_vars$bio_source.y[index] = "IMPACT-BAM_matched"
	}
	
	tmp_vars = tmp_vars %>%
			   mutate(shape = ifelse(bio_source.x=="biopsy_matched" | bio_source.x=="IMPACT-BAM_matched" | bio_source.y=="biopsy_matched" | bio_source.y=="IMPACT-BAM_matched", "Biopsy matched", "Biopsy unmatched")) %>%
			   mutate(fill = ifelse(bio_source.x!=bio_source.y , "Incorrect assignment between replicates", "Called in both replicates")) %>%
			   mutate(fill = ifelse((bio_source.x=="unmatched" & bio_source.y!="unmatched") | (bio_source.x!="unmatched" & bio_source.y=="unmatched"), "Not detected in one replicate", "Called in both replicates")) %>%
			   mutate(fill = ifelse((bio_source.x=="noise" & (bio_source.y!="noise" & bio_source.y!="other")) | (bio_source.y=="noise" & (bio_source.x!="noise" & bio_source.x!="other")), "Not called in one replicate due\nto low quality", fill)) %>%
			   mutate(fill = ifelse(bio_source.x=="VUSo" & bio_source.y == "WBC_matched", "Incorrect assignment between replicates", fill)) %>%
			   mutate(fill = ifelse(bio_source.y=="VUSo" & bio_source.x == "WBC_matched", "Incorrect assignment between replicates", fill))
			   
	pdf(file=paste0("../res/figureS4/", patient_ids[i], "_R2_R3.pdf"), width=6.5, height=7)
	par(mar = c(6.1, 6, 4.1, 1))
	epsilon = 0
	shapes = c("Biopsy matched"=24,
			   "Biopsy unmatched"=21)
	cols = c("Not detected in one replicate"="#D7191C",
			 "Not called in one replicate due\nto low quality"="#FDAE61",
			 "Incorrect assignment between replicates"="#ABDDA4",
			 "Called in both replicates"="#2B83BA")
	plot(1, 1, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, xlim=c(0.01, 100), ylim=c(0.01, 100), log="xy")
	axis(1, at=c(0.01, 0.1, 1, 10, 100), labels=c("0", "0.1", "1", "10", "100"), cex.axis = 1.75, padj = 0.25, lwd=1.85, lwd.ticks=1.75)
	axis(2, at=c(0.01, 0.1, 1, 10, 100), labels=c("0", "0.1", "1", "10", "100"), cex.axis = 1.75, las = 1, lwd=1.85, lwd.ticks=1.75)
	mtext(side = 1, text = "Replicate 2 (%)", line = 4, cex = 1.75)
	mtext(side = 2, text = "Replicate 3 (%)", line = 4, cex = 1.75)
	points(c(.01,100), c(.01,100), type="l", lty=1, lwd=2, col="goldenrod3")
	x = tmp_vars$afnobaq.x
	y = tmp_vars$afnobaq.y
	z1 = as.character(tmp_vars$shape)
	z2 = as.character(tmp_vars$fill)
	points(x, y, pch=shapes[z1], col="black", bg=cols[z2], cex=1.65)
	log10_axis(side=1, at=c(0.01, 0.1, 1, 10, 100), lwd=0, lwd.ticks=1)
	log10_axis(side=2, at=c(0.01, 0.1, 1, 10, 100), lwd=0, lwd.ticks=1)
	dev.off()
	
}

#==================================================
# tabulate %ppa for 6 patients used to test
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
				  dplyr::select(-t1, -t2, -t3, -t4, -g1, -g2, -g3, -g4, -p1, -p2, -p3, -p4) %>%
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
				  dplyr::select(-min_p) %>%
				  mutate(loc = str_c(chrom, ":", pos, "_", ref, ">", alt)) %>%
				  mutate(position_orig = pos, ref_orig = ref, alt_orig = alt, position = pos) %>%
				  mutate(gratio = (adnobaq+2)*(dpgdna+4)/((adgdna+2)*(dpnobaq+4)),
				  		 gzero = adgdna/sqrt(adgdna+2)) %>%
				  mutate(start = position_orig,
                         end = position_orig + 1) %>%
                  genome_left_join(clean_target_region, by = c("chrom", "start", "end")) %>%
                  mutate(chrom = chrom.x) %>%
                  dplyr::select(-c(chrom.x, chrom.y, start.x, end.x, start.y, end.y)) %>%
				  left_join(variants %>% dplyr::select(patient_id, chrom, position_orig, ref_orig, alt_orig, MSK, grail), by=c("patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
				  mutate(MSK = ifelse(is.na(MSK), 0, MSK)) %>%
				  mutate(grail = 1) %>%
				  mutate(study = "TechVal")
				  
feature_names = intersect(colnames(small_vars_plasma), colnames(techval_repeats))
repeat_variants = bind_rows(small_vars_plasma[,feature_names,drop=FALSE] %>%
							filter(!(patient_id %in% patient_ids)),
							techval_repeats[,feature_names,drop=FALSE])

repeat_variants = label_bio_source(repeat_variants)
repeat_variants = left_join(repeat_variants, msk_anno %>% dplyr::select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
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
# variants with vaf <1%
#==================================================
for (i in 1:length(patient_ids)) {
	tmp_r1 = vars_rep1 %>%
			 filter(patient_id == patient_ids[i]) %>%
			 mutate(uuid = paste(patient_id, chrom, position_orig, ref_orig, alt_orig, sep="_"))
	tmp_r0 = vars_rep0 %>%
		  	 filter(patient_id == patient_ids[i]) %>%
		  	 filter(is_nonsyn) %>%
		  	 mutate(uuid = paste(patient_id, chrom, position_orig, ref_orig, alt_orig, sep="_")) %>%
		  	 filter(af_nobaq<1) %>%
		  	 mutate(r = uuid %in% tmp_r1$uuid)
	pander(table(tmp_r0$r, tmp_r0$bio_source))
}
				
#==================================================
# tabulate %ppa for 3 patients used to test
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
				  dplyr::select(-t1, -t2, -t3, -t4, -g1, -g2, -g3, -g4, -p1, -p2, -p3, -p4) %>%
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
				  dplyr::select(-min_p) %>%
				  mutate(loc = str_c(chrom, ":", pos, "_", ref, ">", alt)) %>%
				  mutate(position_orig = pos, ref_orig = ref, alt_orig = alt, position = pos) %>%
				  mutate(gratio = (adnobaq+2)*(dpgdna+4)/((adgdna+2)*(dpnobaq+4)),
				  		 gzero = adgdna/sqrt(adgdna+2)) %>%
				  mutate(start = position_orig,
                         end = position_orig + 1) %>%
                  genome_left_join(clean_target_region, by = c("chrom", "start", "end")) %>%
                  mutate(chrom = chrom.x) %>%
                  dplyr::select(-c(chrom.x, chrom.y, start.x, end.x, start.y, end.y)) %>%
				  left_join(variants %>% dplyr::select(patient_id, chrom, position_orig, ref_orig, alt_orig, MSK, grail), by=c("patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
				  mutate(MSK = ifelse(is.na(MSK), 0, MSK)) %>%
				  mutate(grail = 1) %>%
				  mutate(study = "TechVal")
				  
feature_names = intersect(colnames(small_vars_plasma), colnames(techval_repeats))
repeat_variants = bind_rows(small_vars_plasma[,feature_names,drop=FALSE] %>%
							filter(!(patient_id %in% patient_ids)),
							techval_repeats[,feature_names,drop=FALSE])

repeat_variants = label_bio_source(repeat_variants)
repeat_variants = left_join(repeat_variants, msk_anno %>% dplyr::select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
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
# variants with vaf < 1%
#==================================================
for (i in 1:length(patient_ids)) {
	tmp_r1 = vars_rep1 %>%
			 filter(patient_id == patient_ids[i]) %>%
			 mutate(uuid = paste(patient_id, chrom, position_orig, ref_orig, alt_orig, sep="_"))
	tmp_r0 = vars_rep0 %>%
		  	 filter(patient_id == patient_ids[i]) %>%
		  	 filter(is_nonsyn) %>%
		  	 mutate(uuid = paste(patient_id, chrom, position_orig, ref_orig, alt_orig, sep="_")) %>%
		  	 filter(af_nobaq<1) %>%
		  	 mutate(r = uuid %in% tmp_r1$uuid)
	pander(table(tmp_r0$r, tmp_r0$bio_source))
}

#==================================================
# tabulate variants for 6 patients used to test
# assay reproducibility where relabelling required
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
