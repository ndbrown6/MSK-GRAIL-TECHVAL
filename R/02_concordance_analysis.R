#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/figure2")) {
	dir.create("../res/figure2")
}

#==================================================
# Barplot of mutation burden and sources of mutation
# per patient and cohort
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

ordered_pat_list = list()

cols = c("#97C2E4", "#297AA3")

pdf(file="../res/figure2/allele_frequency_combined.pdf", width=28, height=10)
par(mar = c(6.1, 6, 4.1, 1))
zz = split.screen(figs=matrix(c(0+.05,1/3+.05,0.39,1, 1/3,2/3,0.39,1, 2/3-.05,1-.05,0.39,1, 0+.05,1/3+.05,0,.6, 1/3,2/3,0,.6, 2/3-.05,1-.05,0,.6), nrow=6, ncol=4, byrow=TRUE))
cancer_types = c("Breast", "Lung", "Prostate")
for (i in 1:length(cancer_types)) {
	subj_small_vars = variants %>%
 					  filter(subj_type == cancer_types[i],
 					  bio_source %in% c("biopsy_matched", "biopsy_only") | (bio_source %in% c("IMPACT-BAM_matched", "VUSo") & is_nonsyn))
 
	ordered_patient_ids = subj_small_vars %>%
 						  group_by(patient_id) %>%
 						  mutate(max_af = max(af_nobaq)) %>%
 						  dplyr::select(patient_id, max_af) %>%
 						  unique() %>%
 						  arrange(max_af) %>%
 						  dplyr::select(patient_id) %>%
 						  unlist()
   
 	if (cancer_types[i] == "Control") {
     	sample_ids = tracker_grail %>%
     				 filter(study == "Merlin", tissue == "Healthy") %>% .[["patient_id"]]
   	} else {
     	sample_ids = tracker_grail %>%
     				 filter(patient_id %in% unique(tracker_impact$patient_id)) %>%
     				 filter(study == "TechVal", tissue == cancer_types[i]) %>% .[["patient_id"]]
   	}
   	zero_ids = sample_ids[!sample_ids %in% ordered_patient_ids]
   
 	if (length(zero_ids) > 0) {
     	sample_id_table = data.frame(patient_id = zero_ids,
                                   	 bio_source = "VUSo",
                                   	 af_nobaq = NA)
     	subj_small_vars = bind_rows(sample_id_table, subj_small_vars %>%
                                     dplyr::select(patient_id, bio_source, af_nobaq))
     	ordered_patient_ids = c(zero_ids, ordered_patient_ids)
 	}
   	subj_small_vars = subj_small_vars %>%
     				  mutate(patient_id = factor(patient_id, levels = ordered_patient_ids),
         			  		 bio_source = factor(bio_source, levels = c("VUSo", "IMPACT-BAM_matched", "biopsy_matched", "biopsy_only")))

	ordered_pat_list[[cancer_types[i]]] = ordered_patient_ids
 	
 	cancer_vars_tissue_ref = variants %>%
  						     filter(subj_type != "Control",
  						     bio_source %in% c("biopsy_matched", "biopsy_only", "VUSo"))
 	per_pat_smry = cancer_vars_tissue_ref %>%
   				   filter(bio_source != "VUSo") %>%
   				   group_by(subj_type, patient_id, bio_source) %>%
   				   summarise(num_variants = n()) %>%
   				   ungroup()
 
 	subj_smry = per_pat_smry %>%
   			    filter(subj_type == cancer_types[i])
   			  
 	# add zero-mutation samples
 	pat_ids = unique(subj_smry$patient_id)
  	zero_ids = ordered_pat_list[[cancer_types[i]]][!ordered_pat_list[[cancer_types[i]]] %in% pat_ids]
   	if (length(zero_ids) != 0) {
   		zero_id_table = data.frame(subj_type = cancer_types[i],
     							   patient_id = zero_ids,
     							   bio_source = "biopsy_matched",
     							   num_variants = NA)
     	subj_smry = bind_rows(zero_id_table, subj_smry)
   	}
   
   	plt_smry = data.frame(matrix(0, nrow=length(ordered_patient_ids), ncol=2))
   	rownames(plt_smry) = ordered_patient_ids
   	colnames(plt_smry) = c("biopsy_matched","biopsy_only")
   	tmp = subj_smry$num_variants[subj_smry$patient_id %in% ordered_patient_ids & subj_smry$bio_source=="biopsy_matched"]
   	names(tmp) = subj_smry$patient_id[subj_smry$patient_id %in% ordered_patient_ids & subj_smry$bio_source=="biopsy_matched"]
   	plt_smry[names(tmp),"biopsy_matched"] = tmp
 	tmp = subj_smry$num_variants[subj_smry$patient_id %in% ordered_patient_ids & subj_smry$bio_source=="biopsy_only"]
   	names(tmp) = subj_smry$patient_id[subj_smry$patient_id %in% ordered_patient_ids & subj_smry$bio_source=="biopsy_only"]
   	plt_smry[names(tmp),"biopsy_only"] = tmp
   	plt_smry[is.na(plt_smry)] = 0
   	screen(zz[i+3])
   	plot(1, 1, type="n", xlim=c(1,length(ordered_patient_ids)+5), ylim=c(35,0), xlab="", ylab="", axes=FALSE, frame.plot=FALSE)
   	zzz = barplot(plt_smry[,"biopsy_matched"]+plt_smry[,"biopsy_only"], space=.08, add=TRUE, axes=FALSE, ylab="", col=cols[1])
   	barplot(plt_smry[,"biopsy_matched"], space=.08, add=TRUE, axes=FALSE, ylab="", col=cols[2])
   	if (i==1) {
 		axis(2, at = seq(0,35,by=5), labels=c("", seq(5,35,by=5)), cex.axis = 1.95, las = 1, line = 0, lwd=2)
 		mtext(side = 2, text = "Number of variants", line = 6, cex = 1.95)
   	} else {
   		axis(2, at = seq(0,35,by=5), labels=rep("",8), cex.axis = 1.25, las = 1, line=0, lwd=2)
   	}
   	
   	screen(zz[i])
 	dot_cols = c("VUSo" = "#231F20", "biopsy_matched" = "#231F20", "IMPACT-BAM_matched" = "#231F20", "biopsy_only" = NA)
 	dot_shapes = c("VUSo" = 21, "biopsy_matched" = 21, "IMPACT-BAM_matched" = 25, "biopsy_only" = NA)
 	dot_bg = c("VUSo" = "#3FBC45", "biopsy_matched" = "#297AA3", "IMPACT-BAM_matched" = "#F4AC33", "biopsy_only" = NA)
 	plot(1, 1, type="n", xlim=c(1,length(ordered_patient_ids)+5), ylim=c(.02, 100), log="y", xlab="", ylab="", axes=FALSE, frame.plot=FALSE)
 	for (j in 1:length(ordered_patient_ids)) {
 		points(rep(zzz[j,1],sum(subj_small_vars[,"patient_id"]==ordered_patient_ids[j])), subj_small_vars[subj_small_vars[,"patient_id"]==ordered_patient_ids[j],"af_nobaq",drop=TRUE],
   		   	   pch=dot_shapes[as.character(subj_small_vars[subj_small_vars[,"patient_id"]==ordered_patient_ids[j],"bio_source",drop=TRUE])],
   		   	   col=dot_cols[as.character(subj_small_vars[subj_small_vars[,"patient_id"]==ordered_patient_ids[j],"bio_source",drop=TRUE])],
   		   	   bg=dot_bg[as.character(subj_small_vars[subj_small_vars[,"patient_id"]==ordered_patient_ids[j],"bio_source",drop=TRUE])], lwd=.5)
   	}
   	if (i==1) {
 		axis(2, at = c(.02,.1,.5,1,5,10,50,100), labels = c("0.02","0.1","0.5","1","5","10","50", "100"), cex.axis = 1.95, las = 1, line = 0, lwd=2)
 		mtext(side = 2, text = "VAF (%)", line = 6, cex = 1.95)
 	} else {
   	  	axis(2, at = c(.02,.1,.5,1,5,10,50,100), labels = rep("",8), cex.axis = 1.25, las = 1, line=0, lwd=2)
    }
    title(main=paste0("\n", cancer_types[i]), cex.main=1.85)
  	
}
close.screen(all.screens=TRUE)
dev.off()
 
pdf(file="../res/figure2/allele_frequency_combined_bis.pdf", width=28, height=10)
par(mar = c(6.1, 6, 4.1, 1))
zz = split.screen(figs=matrix(c(0+.05,1/3+.05,0.39,1, 1/3,2/3,0.39,1, 2/3-.05,1-.05,0.39,1, 0+.05,1/3+.05,0,.6, 1/3,2/3,0,.6, 2/3-.05,1-.05,0,.6, 0,1,0,1), nrow=7, ncol=4, byrow=TRUE))
cancer_types = c("Control")
for (i in 1:length(cancer_types)) {
	subj_small_vars = variants %>%
 					  filter(subj_type == cancer_types[i],
 					  bio_source %in% c("biopsy_matched", "VUSo", "biopsy_only"))
 
 	ordered_patient_ids = subj_small_vars %>%
 						  group_by(patient_id) %>%
 						  mutate(max_af = max(af_nobaq)) %>%
 						  dplyr::select(patient_id, max_af) %>%
 						  unique() %>%
 						  arrange(max_af) %>%
 						  dplyr::select(patient_id) %>%
 						  unlist()
   
 	if (cancer_types[i] == "Control") {
     	sample_ids = tracker_grail %>%
     				 filter(study == "Merlin", tissue == "Healthy") %>% .[["patient_id"]]
   	} else {
     	sample_ids = tracker_grail %>%
     				 filter(patient_id %in% unique(tracker_impact$patient_id)) %>%
     				 filter(study == "TechVal", tissue == cancer_types[i]) %>% .[["patient_id"]]
   	}
   	zero_ids = sample_ids[!sample_ids %in% ordered_patient_ids]
   
 	if (length(zero_ids) > 0) {
     	sample_id_table = data.frame(patient_id = zero_ids,
                                   	 bio_source = "VUSo",
                                   	 af_nobaq = NA)
     	subj_small_vars = bind_rows(sample_id_table, subj_small_vars %>%
                                    dplyr::select(patient_id, bio_source, af_nobaq))
    	ordered_patient_ids = c(zero_ids, ordered_patient_ids)
 	}
   	subj_small_vars = subj_small_vars %>%
     				  mutate(patient_id = factor(patient_id, levels = ordered_patient_ids),
         			  bio_source = factor(bio_source, levels = c("VUSo", "biopsy_matched", "biopsy_only")))
   
 	ordered_pat_list[[cancer_types[i]]] = ordered_patient_ids
 	
 	cancer_vars_tissue_ref = variants %>%
   						     filter(subj_type != "Control",
   						     bio_source %in% c("biopsy_matched", "biopsy_only", "VUSo"))
 	per_pat_smry = cancer_vars_tissue_ref %>%
   				   filter(bio_source != "VUSo") %>%
   				   group_by(subj_type, patient_id, bio_source) %>%
   				   summarise(num_variants = n()) %>%
   				   ungroup()
 
 	subj_smry = per_pat_smry %>%
   			    filter(subj_type == cancer_types[i])
   			  
 	# add zero-mutation samples
 	pat_ids = unique(subj_smry$patient_id)
  	zero_ids = ordered_pat_list[[cancer_types[i]]][!ordered_pat_list[[cancer_types[i]]] %in% pat_ids]
   	if (length(zero_ids) != 0) {
   		zero_id_table = data.frame(subj_type = cancer_types[i],
     							   patient_id = zero_ids,
     							   bio_source = "biopsy_matched",
     							   num_variants = NA)
     	subj_smry = bind_rows(zero_id_table, subj_smry)
   	}
   
   	plt_smry = data.frame(matrix(0, nrow=length(ordered_patient_ids), ncol=2))
   	rownames(plt_smry) = ordered_patient_ids
   	colnames(plt_smry) = c("biopsy_matched","biopsy_only")
   	tmp = subj_smry$num_variants[subj_smry$patient_id %in% ordered_patient_ids & subj_smry$bio_source=="biopsy_matched"]
   	names(tmp) = subj_smry$patient_id[subj_smry$patient_id %in% ordered_patient_ids & subj_smry$bio_source=="biopsy_matched"]
   	plt_smry[names(tmp),"biopsy_matched"] = tmp
 	tmp = subj_smry$num_variants[subj_smry$patient_id %in% ordered_patient_ids & subj_smry$bio_source=="biopsy_only"]
   	names(tmp) = subj_smry$patient_id[subj_smry$patient_id %in% ordered_patient_ids & subj_smry$bio_source=="biopsy_only"]
   	plt_smry[names(tmp),"biopsy_only"] = tmp
   	plt_smry[is.na(plt_smry)] = 0
   	screen(zz[i+3])
   	plot(1, 1, type="n", xlim=c(1,length(ordered_patient_ids)+5), ylim=c(35,0), xlab="", ylab="", axes=FALSE, frame.plot=FALSE)
   	zzz = barplot(plt_smry[,"biopsy_matched"]+plt_smry[,"biopsy_only"], space=.08, add=TRUE, axes=FALSE, ylab="", col=cols[1], plot=FALSE)
   	
   	screen(zz[i])
 	dot_cols = c("VUSo" = "#231F20", "biopsy_matched" = "#231F20", "IMPACT-BAM_matched" = "#231F20", "biopsy_only" = NA)
 	dot_shapes = c("VUSo" = 21, "biopsy_matched" = 21, "IMPACT-BAM_matched" = 25, "biopsy_only" = NA)
 	dot_bg = c("VUSo" = "#3FBC45", "biopsy_matched" = "#297AA3", "IMPACT-BAM_matched" = "#F4AC33", "biopsy_only" = NA)
 	plot(1, 1, type="n", xlim=c(1,length(ordered_patient_ids)+5), ylim=c(.02, 100), log="y", xlab="", ylab="", axes=FALSE, frame.plot=FALSE)
 	for (j in 1:length(ordered_patient_ids)) {
 		points(rep(zzz[j,1],sum(subj_small_vars[,"patient_id"]==ordered_patient_ids[j])), subj_small_vars[subj_small_vars[,"patient_id"]==ordered_patient_ids[j],"af_nobaq"],
   		   	   pch=dot_shapes[subj_small_vars[subj_small_vars[,"patient_id"]==ordered_patient_ids[j],"bio_source"]],
   		   	   col=dot_cols[subj_small_vars[subj_small_vars[,"patient_id"]==ordered_patient_ids[j],"bio_source"]],
   		   	   bg=dot_bg[as.character(subj_small_vars[subj_small_vars[,"patient_id"]==ordered_patient_ids[j],"bio_source",drop=TRUE])], lwd=.5)
   	}
   	if (i==1) {
 		axis(2, at = c(.02,.1,.5,1,5,10,50,100), labels = c("0.02","0.1","0.5","1","5","10","50", "100"), cex.axis = 1.95, las = 1, line = 0, lwd=2)
 		mtext(side = 2, text = "VAF (%)", line = 6, cex = 1.95)
 		axis(1, at = zzz[,1], labels=rep("",length(ordered_patient_ids)))
 	} else {
   	  	axis(2, at = c(.02,.1,.5,1,5,10,50,100), labels = rep("",8), cex.axis = 1.25, las = 1, line=0, lwd=2)
    }
    title(main=paste0("\n", cancer_types[i]), cex.main=1.85)
   	
}
close.screen(all.screens=TRUE)
dev.off()

#==================================================
# Barplot of recurrent genes
#==================================================
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
 						 filter(bio_source %in% c("biopsy_matched", "IMPACT-BAM_matched", "VUSo"), gene %in% gene_list, subj_type != "Control") %>%
 						 dplyr::select(subj_type, gene, percent_patient, bio_source)
   
top_cancer_genes_table = top_cancer_genes_table %>%
 						 mutate(gene = factor(gene, levels = top_cancer_genes_ordered),
 						 subj_type = factor(subj_type, levels = c("Breast", "Lung", "Prostate")),
 						 bio_source = factor(bio_source, levels = c("VUSo", "biopsy_matched", "IMPACT-BAM_matched")))

bar_fill_cols = c("biopsy_matched" = "#297AA3", "IMPACT-BAM_matched" = "#F4AC33", "VUSo" = "#3FBC45")
pdf(file="../res/figure2/recurrent_genes_combined.pdf", width=10, height=6.7)
par(mar = c(6.1, 6, 4.1, 1))
zz = split.screen(figs=matrix(c(0,1,2/3-.16,1, 0,1,1/3-.08,2/3+.08, 0,1,0,1/3+.16, 0,1,0,1), nrow=4, ncol=4, byrow=TRUE))
cancer_types = c("Breast", "Lung", "Prostate")
for (i in 1:length(cancer_types)) {
	screen(zz[i])
	x.1 = top_cancer_genes_table %>%
		  filter(subj_type==cancer_types[i], bio_source=="biopsy_matched")
	x.2 = top_cancer_genes_table %>%
		  filter(subj_type==cancer_types[i], bio_source=="IMPACT-BAM_matched")
	x.3 = top_cancer_genes_table %>%
		  filter(subj_type==cancer_types[i], bio_source=="VUSo")
	x = matrix(0, nrow=length(top_cancer_genes_ordered), ncol=3)
 	rownames(x) = top_cancer_genes_ordered
 	colnames(x) = c("biopsy_matched","IMPACT-BAM_matched","VUSo")
 	if (nrow(x.1)!=0) {
 		x[as.character(x.1$gene),"biopsy_matched"] = as.numeric(x.1$percent_patient)
 	}
 	if (nrow(x.2)!=0) {
 		x[as.character(x.2$gene),"IMPACT-BAM_matched"] = as.numeric(x.2$percent_patient)
 	}
 	if (nrow(x.3)!=0) {
 		x[as.character(x.3$gene),"VUSo"] = as.numeric(x.3$percent_patient)
 	}
 	zzz = barplot(t(x), beside=FALSE, names.arg=rep("",nrow(x)), col=bar_fill_cols, axes=FALSE, ylim=c(0,50))
 	axis(2, at = seq(0,50,l=6), labels = seq(0,50,l=6), cex.axis = 1.05, las = 1, line = 0, lwd=1)
 	if (i==3) {
 		axis(side=1, at=zzz, labels=top_cancer_genes_ordered, las=2, font=3, cex=.75, lwd=-1)
 	}
 	title(main=paste0("\n\n",cancer_types[i]), cex.main=1.25)
}
screen(zz[4])
plot(0, 0, type="n", xlim=c(0,10), ylim=c(0,10), xlab="", ylab="", axes=FALSE, frame.plot=FALSE)
mtext(side = 2, text = "% of patients", line = 3, cex = 1.15)
close.screen(all.screens=TRUE)
dev.off()
 
#==================================================
# Violin plot of cfDNA fraction by cancer type
#==================================================
cfdna_frac = read.csv(file=url_ctdna_frac, header=TRUE, sep=",", stringsAsFactors=FALSE) %>%
			 filter(!is.na(ctdna_frac)) %>%
			 mutate(index = order_samples(ID)) %>%
 			 arrange(desc(ctdna_frac)) %>%
			 arrange(index)
cancer_types = c("Breast", "Lung", "Prostate")

pdf(file="../res/figure2/cfdna_fraction_by_cancer_type.pdf", width=7, height=7)
par(mar = c(6.1, 6, 4.1, 1))
plot(0, 0, type="n", xlim=c(0.5,3.5), ylim=c(0.01,1), axes=FALSE, frame.plot=FALSE, xlab="", ylab="")
vioplot(cfdna_frac[cfdna_frac[,3]==1,2], col=awtools::mpalette[3], add=TRUE, at=1)
vioplot(cfdna_frac[cfdna_frac[,3]==2,2], col=awtools::mpalette[2], add=TRUE, at=2)
vioplot(cfdna_frac[cfdna_frac[,3]==3,2], col=awtools::mpalette[4], add=TRUE, at=3)
axis(1, at = 1:3, labels=cancer_types, cex.axis = 1.5, padj = 0.25)
axis(2, at = NULL, cex.axis = 1.5, las = 1)
text(x=2, y=.975, "p = 0.0046" , cex=1.05)
mtext(side = 2, text = "Estimated cfDNA fraction", line = 4, cex = 1.5)
dev.off()

#==================================================
# Box plot plot of cfDNA fraction by number of
# metastatic sites
#==================================================
clinical = read_tsv(file=clinical_file_updated, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
cfdna_frac = read.csv(file=url_ctdna_frac, header=TRUE, sep=",", stringsAsFactors=FALSE) %>%
			 filter(!is.na(ctdna_frac)) %>%
			 mutate(index = order_samples(ID)) %>%
 			 arrange(desc(ctdna_frac)) %>%
			 arrange(index)
number_metastatic_sites = unlist(lapply(strsplit(clinical[,"metastatic_sites",drop=TRUE], split=",", fixed=TRUE), function(x) {length(x)}))
number_metastatic_sites[is.na(clinical[,"metastatic_sites",drop=TRUE])] = 0
names(number_metastatic_sites) = clinical[,"patient_id",drop=TRUE]
number_metastatic_sites = number_metastatic_sites[cfdna_frac[,1]]

pdf(file="../res/figure2/cfdna_fraction_by_number_of_site.pdf", width=7, height=7)
par(mar = c(6.1, 6, 4.1, 1))
x = list(cfdna_frac[number_metastatic_sites==1 | number_metastatic_sites==2,2],
		 cfdna_frac[number_metastatic_sites==3,2],
		 cfdna_frac[number_metastatic_sites>4,2])
y = list(cfdna_frac[number_metastatic_sites==1 | number_metastatic_sites==2,1],
		 cfdna_frac[number_metastatic_sites==3,1],
		 cfdna_frac[number_metastatic_sites>4,1])
for (i in 1:3) {
	y[[i]][grep("VB", y[[i]])] = awtools::mpalette[3]
	y[[i]][grep("VL", y[[i]])] = awtools::mpalette[2]
	y[[i]][grep("VP", y[[i]])] = awtools::mpalette[4]
}
plot(0, 0, type="n", xlim=c(0.5,3.5), ylim=c(0.01,1), axes=FALSE, frame.plot=FALSE, xlab="", ylab="")
boxplot(x[[1]], add=TRUE, at=1, outline=FALSE, axes=FALSE)
boxplot(x[[2]], add=TRUE, at=2, outline=FALSE, axes=FALSE)
boxplot(x[[3]], add=TRUE, at=3, outline=FALSE, axes=FALSE)
for (i in 1:3) {
	points(jitter(rep(i, length(x[[i]])), amount = 0.20), x[[i]], pch = 21, col = "black", bg = y[[i]], cex = 1.5, lwd=.5)
}
axis(1, at = 1:3, labels=c("1-2","3","4+"), cex.axis = 1.5, padj = 0.25)
axis(2, at = NULL, cex.axis = 1.5, las = 1)
rect(xleft=0, ybottom=0.8, xright=2, ytop=1.5, col="white", border="white")
box(lwd=1)
legend("topleft", legend=c("Breast, p = 0.002", "Lung, p = 0.018", "Prostate, p = 0.18"), pch=21, col="black", pt.bg=c(awtools::mpalette[3], awtools::mpalette[2], awtools::mpalette[4]), box.lwd=-1, pt.cex=1.5)
mtext(side = 2, text = "Estimated cfDNA fraction", line = 4, cex = 1.5)
dev.off()

#==================================================
# Trend test ctDNA fraction versus number of sites
#==================================================
x = cfdna_frac[,"ctdna_frac"]
y = number_metastatic_sites

print(cor.test(x=x, y=y, method="kendall", alternative="two.sided", exact=FALSE)$p.value)
print(jonckheere.test(x=x, g=y, alternative="two.sided", nperm=10000)$p.value)
 
## jonckheere-terpstra test by cancer type
for (i in c("VB", "VL", "VP")) {
 	index = grepl(i, names(number_metastatic_sites))
 	x = cfdna_frac[index,"ctdna_frac"]
 	y = number_metastatic_sites[index]
 	print(jonckheere.test(x=x, g=y, alternative="two.sided", nperm=10000)$p.value)
}
for (i in c("VB", "VL", "VP")) {
 	index = grepl(i, names(number_metastatic_sites))
 	x = cfdna_frac[index,"ctdna_frac"]
 	y = number_metastatic_sites[index]
 	y0 = rep(1, length(y))
 	y0[y==3] = 2
 	y0[y>3] = 3
	print(jonckheere.test(x=x, g=y0, alternative="increasing", nperm=10000)$p.value)
}

## kendall correlation test by cancer type
for (i in c("VB", "VL", "VP")) {
 	index = grepl(i, names(number_metastatic_sites))
 	x = cfdna_frac[index,"ctdna_frac"]
 	y = number_metastatic_sites[index]
 	print(cor.test(x=x, y=y, method="kendall", alternative="two.sided", exact = FALSE)$p.value)
}
for (i in c("VB", "VL", "VP")) {
 	index = grepl(i, names(number_metastatic_sites))
 	x = cfdna_frac[index,"ctdna_frac"]
 	y = number_metastatic_sites[index]
 	y0 = rep(1, length(y))
 	y0[y==3] = 2
 	y0[y>3] = 3
 	print(cor.test(x=x, y=y0, method="kendall", alternative="greater", exact = FALSE)$p.value)
}

## mblm association test by cancer type
for (i in c("VB", "VL", "VP")) {
 	index = grepl(i, names(number_metastatic_sites))
 	x = cfdna_frac[index,"ctdna_frac"]
 	y = number_metastatic_sites[index]
	data = data.frame(ctDNA_fraction=x, n_sites=y) 
 
 	z = mblm(ctDNA_fraction ~ n_sites, dataframe=data, repeated=FALSE)
 	p = summary(z)$coefficients[2,4]
 	print(p)
}

#==================================================
# Trend test ctDNA fraction by cancer type
#==================================================
x = cfdna_frac[,"ctdna_frac"]
y = rep(1, length(x))
y[grepl("VL", cfdna_frac[,1])] = 2
y[grepl("VP", cfdna_frac[,1])] = 3

print(cor.test(x=x, y=y, method="kendall", alternative="two.sided", exact=FALSE)$p.value)
print(jonckheere.test(x=x, g=y, alternative="two.sided", nperm=10000)$p.value)
print(kruskal.test(x=x, g=y)$p.value)
 	
#==================================================
# P-values of detection rate
#==================================================
N = c(39, 41, 44)
n = c(37, 31, 36)
 
# Breast versus lung
m1 = matrix(NA, 2, 2)
m1[1,1] = N[1]-n[1]
m1[2,1] = N[2]-n[2]
m1[1,2] = n[1]
m1[2,2] = n[2]
p1 = fisher.test(m1)$p.value
  
# Breast versus prostate
m2 = matrix(NA, 2, 2)
m2[1,1] = N[1]-n[1]
m2[2,1] = N[3]-n[3]
m2[1,2] = n[1]
m2[2,2] = n[3]
p2 = fisher.test(m2)$p.value
  
# Lung versus prostate
m3 = matrix(NA, 2, 2)
m3[1,1] = N[2]-n[2]
m3[2,1] = N[3]-n[3]
m3[1,2] = n[2]
m3[2,2] = n[3]
p3 = fisher.test(m3)$p.value
  
p = p.adjust(c(p1, p2, p3), method="bonferroni")
print(p)
