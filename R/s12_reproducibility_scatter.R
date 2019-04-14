#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/rebuttal")) {
	dir.create("../res/rebuttal")
}

ver_tidy <- "../modified_v11/"
s3_target <- paste0(ver_tidy, "Resources/Bed_files/pan_v2_wo_decoy_wo_iSNP_wo_CNV_IMPACT_common_regions.merged.bed")
target <- read_tsv(file=s3_target, col_names = c("chrom", "start", "end"))
target$is_clean <- TRUE

s3_tracker <- paste0(ver_tidy, "Tracker_files/s3_sample_tracker_TechVal_Merlin.csv")
s3_tracker_retest <- paste0(ver_tidy, "Tracker_files/s3_sample_tracker_TechVal_retest.csv")
tracker <- read_csv(file=s3_tracker, col_types = cols(.default = col_character())) %>%
		   type_convert()
tracker_retest <- read_csv(file=s3_tracker_retest, col_types = cols(.default = col_character())) %>%
		   		  type_convert()

s3_snv <- paste0(ver_tidy, "Variants_Calls/Joined_cfDNA_IMPACT_variants/scored_merged_snvs_20171115.tsv")
snv <- read_tsv(file=s3_snv, col_types = cols(.default = col_character())) %>%
		   		  type_convert() 

s3_retest <- paste0(ver_tidy, "Variants_Calls/Stacked_annotated_retest/TechVal_retest_annotated_stack.tsv")
retest <- read_tsv(file=s3_retest, col_types = cols(.default = col_character())) %>%
		  type_convert()

retest_subj <- unique(tracker_retest$patient_id)

snv_original <- snv %>%
  				filter(patient_id %in% retest_subj) %>%
  				mutate(replicate = "rep1", start = position, end = position +1) %>%
  				genome_left_join(target, by = c("chrom", "start", "end")) %>%
  				mutate(chrom = chrom.x) %>%
  				select(-c(chrom.x, chrom.y, start.x, end.x, start.y, end.y)) %>%
  				filter(is_clean)
  				
kGDNAParams <- data_frame(
  subj_type = c("Healthy", "Breast", "Lung", "Prostate"),
  min_p = c(0.8, 0.79, 0.82, 0.79))
  
'update_filter' <- function(filter, qualnobaq, pgtkxgdna, is_edge, min_p, min_qual = 60, sep = ";")
{
  high_qual <- which(qualnobaq >= min_qual)
  low_qual <- which(qualnobaq < min_qual)
  gdna <- which(pgtkxgdna < min_p)
  edge <- which(is_edge)
  pass <- high_qual %>%
    setdiff(gdna) %>%
    setdiff(edge)
  
  low_qual_string <- sprintf("QUAL_LT_%d", min_qual)
  gdna_strings <- sprintf("PGTKXGDNA_LT_%0.2f", min_p)
  edge_string <- "IS_EDGE"
  
  ids <- data_frame(i = seq_along(filter))
  filters <- bind_rows(
    data_frame(i = pass, filter = "PASS"),
    data_frame(i = low_qual, filter = low_qual_string),
    data_frame(i = gdna, filter = gdna_strings[gdna]),
    data_frame(i = edge, filter = edge_string)
  )
  filters <- filters %>%
    group_by(i) %>%
    summarize(filter = str_c(filter, collapse = sep)) %>%
    ungroup() %>%
    full_join(ids) %>%
    arrange(i)
  filters$filter[filters$i]
}

snv_retest <- retest %>%
  filter(snv) %>%
  mutate_at(vars(matches("qual")), funs(if_else(is.na(.), 0, .))) %>%
  mutate(pedge = ifelse(is.na(pedge), 0, pedge),
         pgtkxgdna = ifelse(is.na(pgtkxgdna), 0, pgtkxgdna)) %>%
  mutate(is_nonsyn = ifelse(
    grepl("missense", annotation) | grepl("stop", annotation) | grepl("frameshift", annotation),
    TRUE, FALSE)) %>%
  mutate(gene = genename,
         hgvs_p = hgvsp,
         position = pos,
         is_snv = snv,
         is_indel = indel,
         grail = 1L,
         study = "TechVal_repeat") %>%
  mutate(subj_type = case_when(
    grepl("VB", .$patient_id) ~ "Breast",
    grepl("VL", .$patient_id) ~ "Lung",
    grepl("VP", .$patient_id) ~ "Prostate"),
    replicate = ifelse(grepl("_B", sample_id), "rep3", "rep2")) %>%
  mutate(filter = "") %>%
  left_join(kGDNAParams) %>%
  mutate_(filter = ~update_filter(
    filter = filter,
    qualnobaq = qualnobaq,
    pgtkxgdna = pgtkxgdna,
    is_edge = isedge,
    min_p = min_p)) %>%
  mutate(start = position,
         end = position +1) %>%
  genome_left_join(target, by = c("chrom", "start", "end")) %>%
  mutate(chrom = chrom.x) %>%
  select(-genename, -pos, -sample_id, -snv, -indel, -min_p,
         -chrom.x, -chrom.y, -start.x, -end.x, -start.y, -end.y) %>%
  filter(is_clean)
  
nm_common <- intersect(names(snv_original), names(snv_retest))
nm_diff <- setdiff(names(snv_original), names(snv_retest))

snv_cfdna <- snv_original %>%
  dplyr::select(nm_common) %>%
  full_join(snv_retest) %>%
  mutate(af_cfdna = adnobaq/dpnobaq) %>%
  filter(is_nonsyn)

snv_tissue <- snv_original %>%
  filter(MSK == "1") %>%
  dplyr::select(patient_id, subj_type, chrom, position, is_snv, is_indel, is_nonsyn, gene, hgvs_p, nm_diff)

snv_cfdna %>%
  group_by(patient_id, replicate) %>%
  summarise(n_pass_nonsyn = length(filter[filter == "PASS"])) %>%
  pandoc.table()
  
cols = c("Not called in one replicate"="#D7191C",
		 "Called in both replicates"="#2B83BA")
  
PID = c("MSK-VB-0023", "MSK-VB-0041", "MSK-VB-0050", "MSK-VL-0028", "MSK-VL-0038", "MSK-VL-0042")
for (i in 1:length(PID)) {

	snv_cfdna_r1 <- snv_cfdna %>%
					filter(patient_id == PID[i]) %>%
	  				filter(replicate == "rep1") %>%
	  				dplyr::select(patient_id, subj_type, chrom, position, gene, is_snv, is_indel, is_nonsyn, grail_rep1 = grail, adnobaq_rep1 = adnobaq, dpnobaq_rep1 = dpnobaq, af_cfdna_rep1 = af_cfdna, filter_rep1 = filter)

	snv_cfdna_r2 <- snv_cfdna %>% 
					filter(patient_id == PID[i]) %>%
					filter(replicate == "rep2") %>%
					dplyr::select(patient_id, subj_type, chrom, position, gene, is_snv, is_indel, is_nonsyn, grail_rep2 = grail, adnobaq_rep2 = adnobaq, dpnobaq_rep2 = dpnobaq, af_cfdna_rep2 = af_cfdna, filter_rep2 = filter)
         
	snv_biopsy <- snv_tissue %>%
				  filter(patient_id == PID[i]) %>%
				  dplyr::select(patient_id, subj_type, chrom, position, hgvs_p, is_snv, is_indel, is_nonsyn, MSK, af_tissue = c_af) %>%
				  mutate(af_tissue = af_tissue/100)

	all_snv <- full_join(snv_cfdna_r1, snv_cfdna_r2, by = c("patient_id", "subj_type", "chrom", "position", "gene", "is_snv", "is_indel", "is_nonsyn")) %>%
			   filter(filter_rep1 == "PASS" | filter_rep2 == "PASS") %>%
			   full_join(snv_biopsy) %>%
	  		   mutate_at(vars(matches("af_cfdna")), funs(if_else(is.na(.), 0, .))) %>%
			   filter(filter_rep1 == "PASS" | filter_rep2 == "PASS") %>%
			   mutate(called_in_cfDNA = case_when(
			  			filter_rep1 == "PASS" & filter_rep2 == "PASS" ~ "Called in both replicates",
			  			TRUE ~ "Not called in one replicate")) %>%
	  		   mutate(called_in_tissue = case_when(
	    				MSK == 1 ~ "Biopsy matched",
	    				is.na(MSK) ~ "Biopsy unmatched")) %>%
	    	   mutate(af_cfdna_rep1 = af_cfdna_rep1*100) %>%
	    	   mutate(af_cfdna_rep2 = af_cfdna_rep2*100) %>%
	    	   mutate(af_cfdna_rep1 = ifelse(af_cfdna_rep1==0, 0.01, af_cfdna_rep1)) %>%
	    	   mutate(af_cfdna_rep2 = ifelse(af_cfdna_rep2==0, 0.01, af_cfdna_rep2))

	plot.0 = ggplot(all_snv, aes(x = af_cfdna_rep1, y = af_cfdna_rep2, shape = called_in_tissue, fill = called_in_cfDNA)) +
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
		 
	pdf(file=paste0("../res/rebuttal/", PID[i],"_R1_R2.pdf"), width=5.5, height=6.5)
	print(plot.0)
	dev.off()
}

PID = c("MSK-VB-0023", "MSK-VL-0028", "MSK-VL-0042")
for (i in 1:length(PID)) {

	snv_cfdna_r1 <- snv_cfdna %>%
					filter(patient_id == PID[i]) %>%
	  				filter(replicate == "rep1") %>%
	  				dplyr::select(patient_id, subj_type, chrom, position, gene, is_snv, is_indel, is_nonsyn, grail_rep1 = grail, adnobaq_rep1 = adnobaq, dpnobaq_rep1 = dpnobaq, af_cfdna_rep1 = af_cfdna, filter_rep1 = filter)

	snv_cfdna_r2 <- snv_cfdna %>% 
					filter(patient_id == PID[i]) %>%
					filter(replicate == "rep2") %>%
					dplyr::select(patient_id, subj_type, chrom, position, gene, is_snv, is_indel, is_nonsyn, grail_rep2 = grail, adnobaq_rep2 = adnobaq, dpnobaq_rep2 = dpnobaq, af_cfdna_rep2 = af_cfdna, filter_rep2 = filter)
					
	snv_cfdna_r3 <- snv_cfdna %>% 
					filter(patient_id == PID[i]) %>%
					filter(replicate == "rep3") %>%
					dplyr::select(patient_id, subj_type, chrom, position, gene, is_snv, is_indel, is_nonsyn, grail_rep3 = grail, adnobaq_rep3 = adnobaq, dpnobaq_rep3 = dpnobaq, af_cfdna_rep3 = af_cfdna, filter_rep3 = filter)
         
	snv_biopsy <- snv_tissue %>%
				  filter(patient_id == PID[i]) %>%
				  dplyr::select(patient_id, subj_type, chrom, position, hgvs_p, is_snv, is_indel, is_nonsyn, MSK, af_tissue = c_af) %>%
				  mutate(af_tissue = af_tissue/100)

	all_snv <- full_join(snv_cfdna_r1, snv_cfdna_r3, by = c("patient_id", "subj_type", "chrom", "position", "gene", "is_snv", "is_indel", "is_nonsyn")) %>%
			   filter(filter_rep1 == "PASS" | filter_rep3 == "PASS") %>%
			   full_join(snv_biopsy) %>%
	  		   mutate_at(vars(matches("af_cfdna")), funs(if_else(is.na(.), 0, .))) %>%
			   filter(filter_rep1 == "PASS" | filter_rep3 == "PASS") %>%
			   mutate(called_in_cfDNA = case_when(
			  			filter_rep1 == "PASS" & filter_rep3 == "PASS" ~ "Called in both replicates",
			  			TRUE ~ "Not called in one replicate")) %>%
	  		   mutate(called_in_tissue = case_when(
	    				MSK == 1 ~ "Biopsy matched",
	    				is.na(MSK) ~ "Biopsy unmatched")) %>%
	    	   mutate(af_cfdna_rep1 = af_cfdna_rep1*100) %>%
	    	   mutate(af_cfdna_rep3 = af_cfdna_rep3*100) %>%
	    	   mutate(af_cfdna_rep1 = ifelse(af_cfdna_rep1==0, 0.01, af_cfdna_rep1)) %>%
	    	   mutate(af_cfdna_rep3 = ifelse(af_cfdna_rep3==0, 0.01, af_cfdna_rep3))

	plot.0 = ggplot(all_snv, aes(x = af_cfdna_rep1, y = af_cfdna_rep3, shape = called_in_tissue, fill = called_in_cfDNA)) +
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
		 
	pdf(file=paste0("../res/rebuttal/", PID[i],"_R1_R3.pdf"), width=5.5, height=6.5)
	print(plot.0)
	dev.off()
	
	all_snv <- full_join(snv_cfdna_r2, snv_cfdna_r3, by = c("patient_id", "subj_type", "chrom", "position", "gene", "is_snv", "is_indel", "is_nonsyn")) %>%
			   filter(filter_rep2 == "PASS" | filter_rep3 == "PASS") %>%
			   full_join(snv_biopsy) %>%
	  		   mutate_at(vars(matches("af_cfdna")), funs(if_else(is.na(.), 0, .))) %>%
			   filter(filter_rep2 == "PASS" | filter_rep3 == "PASS") %>%
			   mutate(called_in_cfDNA = case_when(
			  			filter_rep2 == "PASS" & filter_rep3 == "PASS" ~ "Called in both replicates",
			  			TRUE ~ "Not called in one replicate")) %>%
	  		   mutate(called_in_tissue = case_when(
	    				MSK == 1 ~ "Biopsy matched",
	    				is.na(MSK) ~ "Biopsy unmatched")) %>%
	    	   mutate(af_cfdna_rep2 = af_cfdna_rep2*100) %>%
	    	   mutate(af_cfdna_rep3 = af_cfdna_rep3*100) %>%
	    	   mutate(af_cfdna_rep2 = ifelse(af_cfdna_rep2==0, 0.01, af_cfdna_rep2)) %>%
	    	   mutate(af_cfdna_rep3 = ifelse(af_cfdna_rep3==0, 0.01, af_cfdna_rep3))

	plot.0 = ggplot(all_snv, aes(x = af_cfdna_rep2, y = af_cfdna_rep3, shape = called_in_tissue, fill = called_in_cfDNA)) +
			 geom_abline(linetype = 1, color = "goldenrod3") +
			 geom_point(alpha=1, size=3.5) +
			 scale_fill_manual(values = cols) +
			 scale_shape_manual(values = c(24, 21)) +
			 facet_wrap(~patient_id) +
			 theme_bw(base_size=15) +
			 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.3, 0.8), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
			 labs(x="\nReplicate 2 (%)\n", y="Replicate 3 (%)\n") +
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
		 
	pdf(file=paste0("../res/rebuttal/", PID[i],"_R2_R3.pdf"), width=5.5, height=6.5)
	print(plot.0)
	dev.off()
	
}

#==================================================
# Updated assay reproducibility
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
			

patient_ids = c("MSK-VB-0050", "MSK-VB-0041", "MSK-VL-0028", "MSK-VL-0042", "MSK-VB-0023", "MSK-VL-0038")			
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


cols = c("Not detected in one replicate"="#D7191C",
		 "Not called in one replicate due\nto low quality"="#FDAE61",
		 "Incorrect assignment between replicates"="#ABDDA4",
		 "Called in both replicates"="#2B83BA"
		 )
		 
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
			   mutate(fill = ifelse((bio_source.x=="noise" & (bio_source.y!="noise" & bio_source.y!="other")) | (bio_source.y=="noise" & (bio_source.x!="noise" & bio_source.x!="other")), "Not called in one replicate due\nto low quality", fill))
			   
	
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
		 
	pdf(file=paste0("../res/rebuttal/", patient_ids[i], "_R1_R2.pdf"), width=5.5, height=6.5)
	print(plot.0)
	dev.off()
	
}

patient_ids = c("MSK-VL-0028", "MSK-VL-0042", "MSK-VB-0023")			
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
			   mutate(fill = ifelse((bio_source.x=="noise" & (bio_source.y!="noise" & bio_source.y!="other")) | (bio_source.y=="noise" & (bio_source.x!="noise" & bio_source.x!="other")), "Not called in one replicate due\nto low quality", fill))
			   	
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
		 
	pdf(file=paste0("../res/rebuttal/", patient_ids[i], "_R1_R3.pdf"), width=5.5, height=6.5)
	print(plot.0)
	dev.off()
	
}

patient_ids = c("MSK-VL-0028", "MSK-VL-0042", "MSK-VB-0023")			
all_vars = full_join(vars_rep0 %>% mutate(replicate = 2),
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
			   mutate(fill = ifelse((bio_source.x=="noise" & (bio_source.y!="noise" & bio_source.y!="other")) | (bio_source.y=="noise" & (bio_source.x!="noise" & bio_source.x!="other")), "Not called in one replicate due\nto low quality", fill))
			   	
	plot.0 = ggplot(tmp_vars, aes(x = afnobaq.x, y = afnobaq.y, shape = shape, fill = fill)) +
			 geom_abline(linetype = 1, color = "goldenrod3") +
			 geom_point(alpha=1, size=3.5) +
			 scale_fill_manual(values = cols) +
			 scale_shape_manual(values = c(24, 21)) +
			 facet_wrap(~patient_id) +
			 theme_bw(base_size=15) +
			 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.3, 0.8), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
			 labs(x="\nReplicate 2 (%)\n", y="Replicate 3 (%)\n") +
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
		 
	pdf(file=paste0("../res/rebuttal/", patient_ids[i], "_R2_R3.pdf"), width=5.5, height=6.5)
	print(plot.0)
	dev.off()
	
}
