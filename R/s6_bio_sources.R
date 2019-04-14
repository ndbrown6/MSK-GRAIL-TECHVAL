#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
source('config.R')

if (!dir.exists("../res/figureS6")) {
	dir.create("../res/figureS6")
}

#==================================================
# Barplot of mutation burden and sources of mutation
# per patient and cohort
#==================================================
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
		   
snv_vars = read_tsv(snv_file$scored, col_types = cols(.default = col_character())) %>%
		   type_convert()
		   
indel_vars = read_tsv(indel_file$scored, col_types = cols(.default = col_character())) %>%
			 type_convert()

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

per_pat_smry = variants %>%
			   filter(!(bio_source %in% c("noise","germline","other")), is_nonsyn) %>%
			   group_by(subj_type, patient_id, bio_source) %>%
			   summarise(num = n()) %>%
			   ungroup
  
per_pat_smry_all = variants %>%
				   filter(!(bio_source %in% c("noise","germline","other")), is_nonsyn) %>%
				   group_by(subj_type, patient_id) %>%
				   summarise(all = n()) %>%
				   ungroup

per_pat_smry = left_join(per_pat_smry, per_pat_smry_all)

all_patient_table = small_vars_plasma %>%
					distinct(subj_type, patient_id)
all_patient_table = cbind.data.frame(subj_type = rep(all_patient_table$subj_type, 4),
                                     patient_id = rep(all_patient_table$patient_id, 4),
                                     bio_source = rep(c("WBC_matched",
                                                        "VUSo",
                                                        "biopsy_matched",
                                                        "IMPACT-BAM_matched"),
                                                      each = nrow(all_patient_table)))
per_pat_smry = left_join(all_patient_table, per_pat_smry)
per_pat_smry[is.na(per_pat_smry)] = 0
per_pat_smry = per_pat_smry %>%
			   mutate(orig_num = num, num = round(num / total_bed_Mb, 0))
cohort_stats = per_pat_smry %>%
			   group_by(subj_type, bio_source) %>%
			   summarise(cohort_size = n(), positive_scored = sum(num > 0), percent_of_samples = 100 * sum(num > 0) / n()) %>%
			   ungroup()

pdf(file="../res/figureS6/barplot_bio_sources_hyper.pdf", width=12, height=9)
par(mar = c(6.1, 6, 4.1, 1))
cols = c("WBC_matched"="#EA7180", "biopsy_matched" = "#297AA3", "IMPACT-BAM_matched" = "#F4AC33", "VUSo" = "#3FBC45")
zz = split.screen(figs=matrix(c(0.05,.55,.43,.93, .46,.96,.43,.93, 0.05,.55,0.09,.59, .46,.96,0.09,.59, 0,1,0,1), nrow=5, ncol=4, byrow=TRUE))
subj = c("Control", "Breast", "Lung", "Prostate")
for (i in 1:length(subj)) {
	screen(zz[i])
  	subj_smry = per_pat_smry %>%
  				filter(subj_type == subj[i])
  
  	wbc_pats = subj_smry %>%
  			   filter(bio_source == "WBC_matched") %>%
  			   dplyr::select(patient_id, num, all)
  
  	ordered_patient_ids = wbc_pats %>%
  						  arrange(num, all) %>%
  						  .[["patient_id"]]
  
  	subj_smry = subj_smry %>%
  			    mutate(num = ifelse(bio_source == "WBC_matched", -num, num),
  			  		   patient_id = factor(patient_id, levels = ordered_patient_ids),
  			  		   bio_source = factor(bio_source, levels = c("WBC_matched", "biopsy_matched", "VUSo", "IMPACT-BAM_matched")))
  
  	wbc_matched = subj_smry %>%
  				  filter(bio_source=="WBC_matched")
  	vuso = subj_smry %>%
  		   filter(bio_source=="VUSo")
  	biopsy_matched = subj_smry %>%
  		   		     filter(bio_source=="biopsy_matched")
  	biopsy_subthreshold = subj_smry %>%
  		   		     	  filter(bio_source=="IMPACT-BAM_matched")
  	
  	index = order(biopsy_matched[,"num"]+vuso[,"num"])
  	wbc_matched = wbc_matched[index,,drop=FALSE]
  	vuso = vuso[index,,drop=FALSE]
  	biopsy_matched = biopsy_matched[index,,drop=FALSE]
  	biopsy_subthreshold = biopsy_subthreshold[index,,drop=FALSE]
  	
  	index = order(wbc_matched[,"num"], decreasing=TRUE)
  	wbc_matched = wbc_matched[index,,drop=FALSE]
  	vuso = vuso[index,,drop=FALSE]
  	biopsy_matched = biopsy_matched[index,,drop=FALSE]
  	biopsy_subthreshold = biopsy_subthreshold[index,,drop=FALSE]
  	
  	zzz = barplot(wbc_matched[,"num"]-.9, col=transparentRgb(cols["WBC_matched"], 175), border="black", space=.16, axes=FALSE, ylim=c(-30,150), lwd=.01)
  	abline(h=0, col="white", lwd=4)
  	barplot(ifelse((vuso[,"num"]+biopsy_matched[,"num"]+biopsy_subthreshold[,"num"])<110, (vuso[,"num"]+biopsy_matched[,"num"]+biopsy_subthreshold[,"num"]), 120), col=cols["VUSo"], border="black", space=.16, add=TRUE, axes=FALSE, lwd=.01)
  	barplot(biopsy_matched[,"num"]+biopsy_subthreshold[,"num"], col=cols["IMPACT-BAM_matched"], border="black", space=.16, add=TRUE, axes=FALSE, lwd=.01)
  	barplot(biopsy_matched[,"num"], col=cols["biopsy_matched"], border="black", space=.16, add=TRUE, axes=FALSE, lwd=.01)
  	if (i==1 | i==3) {
  		axis(2, at = c(-30,-15,0,25,50,75,100,115), labels=c(30,15,0,25,50,75,100,110), cex.axis = 1, las = 1)
  		axis(2, at = c(130, 150), labels=c(550,600), cex.axis = 1, las = 1)
  	} else {
		axis(2, at = c(-30,-15,0,25,50,75,100,115), labels=rep("",8), cex.axis = 1, las = 1)
		axis(2, at = c(130, 150), labels=c("",""), cex.axis = 1.4, las = 1)
	}
	axis(2, at = c(100, 150), labels=c("",""), cex.axis = 1, las = 1, lwd.ticks=-1)
  	axis.break(axis=2, breakpos=123, pos=NA,bgcol="white",breakcol="black",style="slash",brw=0.02)
	if (i==2) {
		rect(xleft=zzz[36]-.565, xright=zzz[36]+.565, ybottom=130, ytop=130+(587-550)-20, col=cols["VUSo"], border="black")
	}
	title(main=subj[i])
}
screen(zz[5])
plot(0,0, type="n", xlim=c(0,10), ylim=c(-1,1), xlab="", ylab="", axes=FALSE, frame.plot=FALSE)
mtext(side = 2, text = expression("Somatic cfDNA variants / Mb"), line = 1.5, cex = 1.4)
close.screen(all.screens=TRUE)
dev.off()