#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/figure2")) {
	dir.create("../res/figure2")
}

if (!dir.exists("../res/etc/Source_Data_Fig_2")) {
	dir.create("../res/etc/Source_Data_Fig_2")
}

#==================================================
# barplot of recurrent genes including
# hypermutators
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
		   		  
# sum the number of times a gene occurs in a bio_source in each tissue type
subj_num_smry = variants %>%
				group_by(subj_type) %>%
    			dplyr::summarize(subj_num = length(unique(patient_id))) %>%
    			ungroup()
  
gene_recurrences = variants %>%
				   filter(!is.na(gene), bio_source %in% c("biopsy_matched", "biopsy_only") | (bio_source %in% c("IMPACT-BAM_matched", "VUSo") & is_nonsyn)) %>%
 				   group_by(bio_source, subj_type, gene) %>%
 				   dplyr::summarize(number = n(), num_patient = length(unique(patient_id))) %>%
 				   ungroup %>%
 				   left_join(subj_num_smry) %>%
 				   mutate(percent_patient = round(num_patient / subj_num * 100, 0))
				   
gene_list = c()
for (subj in subj_num_smry$subj_type[subj_num_smry$subj_type != "Control"]) {
	top_cancer_gene_ids = gene_recurrences %>%
						  filter(subj_type == subj,
						  bio_source %in% c("biopsy_matched", "biopsy_only")) %>%
 						  group_by(gene) %>%
 						  dplyr::summarize(all = sum(percent_patient)) %>%
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
 						   dplyr::summarize(perc = sum(percent_patient)) %>%
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
 	zzz = barplot(t(x), beside=FALSE, names.arg=rep("",nrow(x)), col=variant_cols[colnames(x)], axes=FALSE, ylim=c(0,50))
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

export_x = top_cancer_genes_table %>%
		   dplyr::select(`tissue` = `subj_type`,
		   				 `gene_id` = `gene`,
		   				 `percent_patient` = `percent_patient`,
		   				 `bio_source` = `bio_source`)
write_tsv(export_x, path="../res/etc/Source_Data_Fig_2/Fig_2b.tsv", append=FALSE, col_names=TRUE)

#==================================================
# barplot of mutation burden and sources of mutation
# per patient and cohort including hypermutators
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
     	sample_id_table = data_frame(patient_id = zero_ids,
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
   				   dplyr::summarize(num_variants = n()) %>%
   				   ungroup()
  
 	subj_smry = per_pat_smry %>%
   			    filter(subj_type == cancer_types[i])
   			  
 	# add zero-mutation samples
  	pat_ids = unique(subj_smry$patient_id)
  	zero_ids = ordered_pat_list[[cancer_types[i]]][!ordered_pat_list[[cancer_types[i]]] %in% pat_ids]
   	if (length(zero_ids) != 0) {
   		zero_id_table = data_frame(subj_type = cancer_types[i],
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
   	zzz = barplot(plt_smry[,"biopsy_matched"]+plt_smry[,"biopsy_only"], space=.08, add=TRUE, axes=FALSE, ylab="", col=transparent_rgb(variant_cols["biopsy_only"], 155))
   	barplot(plt_smry[,"biopsy_matched"], space=.08, add=TRUE, axes=FALSE, ylab="", col=variant_cols["biopsy_matched"])
   	if (i==1) {
 		axis(2, at = seq(0,35,by=5), labels=c("", seq(5,35,by=5)), cex.axis = 1.95, las = 1, line = 0, lwd=2)
 		mtext(side = 2, text = "Number of variants", line = 6, cex = 1.95)
   	} else {
   		axis(2, at = seq(0,35,by=5), labels=rep("",8), cex.axis = 1.25, las = 1, line=0, lwd=2)
   	}
   	
   	screen(zz[i])
  	dot_cols = c("VUSo" = "#231F20", "biopsy_matched" = "#231F20", "IMPACT-BAM_matched" = "#231F20", "biopsy_only" = NA)
  	dot_shapes = c("VUSo" = 21, "biopsy_matched" = 21, "IMPACT-BAM_matched" = 25, "biopsy_only" = NA)
  	dot_bg = c(variant_cols[c("VUSo", "biopsy_matched", "IMPACT-BAM_matched")], "biopsy_only" = NA)
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
    
    export_x = plt_smry %>%
    		   mutate(`patient_id` = rownames(plt_smry)) %>%
    		   dplyr::select(`patient_id`,
    		   				 `biopsy_matched`,
    		   				 `biopsy_only`)
    export_y = subj_small_vars %>%
    		   dplyr::select(`patient_id` = `patient_id`,
    		   				 `bio_source` = `bio_source`,
    		   				 `af_cfdna` = `af_nobaq`)
    write_tsv(export_x, path=paste0("../res/etc/Source_Data_Fig_2/Fig_2d_1_", i, ".tsv"), append=FALSE, col_names=TRUE)
    write_tsv(export_y, path=paste0("../res/etc/Source_Data_Fig_2/Fig_2d_2_", i, ".tsv"), append=FALSE, col_names=TRUE)
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
   				   dplyr::summarize(num_variants = n()) %>%
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
   	zzz = barplot(plt_smry[,"biopsy_matched"]+plt_smry[,"biopsy_only"], space=.08, add=TRUE, axes=FALSE, ylab="", col=transparent_rgb(variant_cols["biopsy_only"]), plot=FALSE)
   	
   	screen(zz[i])
 	dot_cols = c("VUSo" = "#231F20", "biopsy_matched" = "#231F20", "IMPACT-BAM_matched" = "#231F20", "biopsy_only" = NA)
 	dot_shapes = c("VUSo" = 21, "biopsy_matched" = 21, "IMPACT-BAM_matched" = 25, "biopsy_only" = NA)
 	dot_bg = c(variant_cols[c("VUSo", "biopsy_matched", "IMPACT-BAM_matched")], "biopsy_only" = NA)
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
    
    export_x = subj_small_vars %>%
    		   dplyr::select(`patient_id` = `patient_id`,
    		   				 `bio_source` = `bio_source`,
    		   				 `af_cfdna` = `af_nobaq`)
    write_tsv(export_x, path="../res/etc/Source_Data_Fig_2/Fig_2c.tsv", append=FALSE, col_names=TRUE)
}
close.screen(all.screens=TRUE)
dev.off()

#==================================================
# violin plot of ctdna fraction by cancer type
#==================================================
cfdna_frac = read.csv(file=url_ctdna_frac, header=TRUE, sep=",", stringsAsFactors=FALSE) %>%
			 filter(!is.na(ctdna_frac)) %>%
			 mutate(index = order_samples(ID)) %>%
 			 arrange(desc(ctdna_frac)) %>%
			 arrange(index)
cancer_types = c("Breast", "Lung", "Prostate")

pdf(file="../res/figure2/cfdna_fraction_by_cancer_type.pdf", width=5.5, height=7)
par(mar = c(6.1, 6, 4.1, 1))
plot(0, 0, type="n", xlim=c(0.5,3.5), ylim=c(0.01,1), axes=FALSE, frame.plot=FALSE, xlab="", ylab="")
beeswarm(cfdna_frac[cfdna_frac[,3]==1,2], method="swarm", col=cohort_cols["Breast"], add=TRUE, pch=19, at=1, cex=1.75, corral="random", corralWidth=.9)
beeswarm(cfdna_frac[cfdna_frac[,3]==2,2], method="swarm", col=cohort_cols["Lung"], add=TRUE, at=2, pch=19, cex=1.75, corral="random", corralWidth=.85)
beeswarm(cfdna_frac[cfdna_frac[,3]==3,2], method="swarm", col=cohort_cols["Prostate"], add=TRUE, at=3, pch=19, cex=1.75, corral="random", corralWidth=.8)
axis(1, at = 1:3, labels=cancer_types, cex.axis = 1.5, padj = 0.25)
axis(2, at = NULL, cex.axis = 1.5, las = 1)
text(x=2, y=.995, "p = 0.0046" , cex=1.05)
box()
mtext(side = 2, text = "ctDNA fraction", line = 4, cex = 1.5)
dev.off()

pdf(file="../res/figure2/cfdna_fraction_by_cancer_type_bis.pdf", width=5.5, height=7)
par(mar = c(6.1, 6, 4.1, 1))
plot(0, 0, type="n", xlim=c(0.5,3.5), ylim=c(0.01,1), axes=FALSE, frame.plot=FALSE, xlab="", ylab="")
vioplot(cfdna_frac[cfdna_frac[,3]==1,2], at = 1, add = TRUE, col=cohort_cols["Breast"], border = "black", rectCol = "black", lineCol = "black", pchMed = 19, colMed = "white")
vioplot(cfdna_frac[cfdna_frac[,3]==2,2], at = 2, add = TRUE, col=cohort_cols["Lung"], border = "black", rectCol = "black", lineCol = "black", pchMed = 19, colMed = "white")
vioplot(cfdna_frac[cfdna_frac[,3]==3,2], at = 3, add = TRUE, col=cohort_cols["Prostate"], border = "black", rectCol = "black", lineCol = "black", pchMed = 19, colMed = "white")
axis(1, at = 1:3, labels=cancer_types, cex.axis = 1.5, padj = 0.25)
axis(2, at = NULL, cex.axis = 1.5, las = 1)
text(x=2, y=.995, "p = 0.0046" , cex=1.05)
mtext(side = 2, text = "ctDNA fraction", line = 4, cex = 1.5)
dev.off()

export_x = cfdna_frac %>%
		   dplyr::select(`patient_id` = `ID`,
		   				 `tissue` = `index`,
		   				 `ctdna_fraction` = `ctdna_frac`) %>%
		   mutate(tissue = case_when(
		   			tissue == 1 ~ "Breast",
		   			tissue == 2 ~ "Lung",
		   			tissue == 3 ~ "Prostate"))
write_tsv(export_x, path="../res/etc/Source_Data_Fig_2/Fig_2f.tsv", append=FALSE, col_names=TRUE)
 
#==================================================
# box plot of ctdna fraction by
# number of metastatic sites
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

x = list(cfdna_frac[number_metastatic_sites==1 | number_metastatic_sites==2,2],
		 cfdna_frac[number_metastatic_sites==3,2],
		 cfdna_frac[number_metastatic_sites>4,2])
y = list(cfdna_frac[number_metastatic_sites==1 | number_metastatic_sites==2,1],
		 cfdna_frac[number_metastatic_sites==3,1],
		 cfdna_frac[number_metastatic_sites>4,1])
for (i in 1:3) {
	y[[i]][grep("VB", y[[i]])] = "Breast"
	y[[i]][grep("VL", y[[i]])] = "Lung"
	y[[i]][grep("VP", y[[i]])] = "Prostate"
}
z = list()
for (i in 1:3) {
	z[[i]] = rep(i, length(x[[i]]))
}

tmp.0 = data_frame(
			ctdna = as.numeric(unlist(x)),
			tissue = factor(unlist(y), levels=c("Breast", "Lung", "Prostate"), ordered=TRUE),
			cat3 = factor(unlist(z), levels=c("1", "2", "3"), ordered=TRUE))
				 
plot.0 = ggplot(tmp.0, aes(x = tissue, y = ctdna, color = cat3, group=interaction(tissue, cat3))) + 
	     geom_boxplot(alpha=1, outlier.size=NA, outlier.shape=NA, fill="white", width=.8) +
	     scale_colour_manual(values=c("1"="black", "2"="black", "3"="black")) +
	     geom_jitter(
		 	aes(x = tissue, y = ctdna, group = cat3, fill=tissue),
  		    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
  		    alpha=1, size=3, shape=21,
  		 ) +
  		 scale_fill_manual(values=cohort_cols) +
	     theme_classic(base_size=15) +
 	     theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10), legend.position="none") +
 	     labs(x="", y="ctDNA fraction\n") +
 	     coord_cartesian(ylim = c(0,1))
 	     
pdf(file="../res/figure2/cfdna_fraction_by_number_of_site.pdf", width=8.5, height=5)
print(plot.0)
dev.off()

#==================================================
# box plot of ctdna fraction by
# disease volume
#==================================================
ctdna_fraction = read_csv(file=url_ctdna_frac, col_types = cols(.default = col_character()))  %>%
 				 type_convert() %>%
 				 rename(GRAIL_ID = ID)

vol_breast = read.xls(xls=url_volumetric_data, sheet=1, stringsAsFactors=FALSE) %>%
			 type_convert() %>%
			 mutate(Date_Tissue = as.character(Date_Tissue)) %>%
			 mutate(Date_Most_Recent_Scan = as.character(Date_Most_Recent_Scan)) %>%
			 mutate(Date_of_Last_Treatment = as.character(Date_of_Last_Treatment))

vol_lung = read.xls(xls=url_volumetric_data, sheet=2, stringsAsFactors=FALSE) %>%
		   type_convert() %>%
		   mutate(Date_Tissue = as.character(Date_Tissue)) %>%
		   mutate(Date_Most_Recent_Scan = as.character(Date_Most_Recent_Scan)) %>%
		   mutate(Date_of_Last_Treatment = as.character(Date_of_Last_Treatment)) %>%
		   mutate(Date_Dx = as.character(Date_Dx)) %>%
		   mutate(Date_Metastatic = as.character(Date_Metastatic)) %>%
		   mutate(Date_Tissue = as.character(Date_Tissue))
		   
vol_prostate = read.xls(xls=url_volumetric_data, sheet=3, stringsAsFactors=FALSE) %>%
		       type_convert() %>%
		       mutate(BSI_Value = ifelse(BSI_Value==">13", 14, BSI_Value)) %>%
		       mutate(BSI_Value = ifelse(BSI_Value=="n/a", NA, BSI_Value)) %>%
		       mutate(BSI_Value = as.numeric(BSI_Value)) %>%
		       rename(BSI_Volume = BSI_Value) %>%
		       mutate(Date_Dx = as.character(Date_Dx)) %>%
		       mutate(Date_Metastatic = as.character(Date_Metastatic)) %>%
		       mutate(Date_Tissue = as.character(Date_Tissue))

volumetric_data = bind_rows(vol_breast, vol_lung) %>%
				  bind_rows(vol_prostate) %>%
				  left_join(ctdna_fraction, by="GRAIL_ID")
				  
volumetric_data = volumetric_data %>%
				  mutate(tot_volu = apply( volumetric_data[,which(grepl('volume', colnames(volumetric_data), ignore.case = T))] , 1 , function(x) sum(as.numeric(x), na.rm=T) )) %>%
				  mutate(tot_mets = apply( volumetric_data[,which(!grepl('volume', colnames(volumetric_data), ignore.case = T) & grepl('Lesions|Lymph_Nodes|Pleural_Disease', colnames(volumetric_data), ignore.case = T))] , 1 , function(x) sum(as.numeric(x), na.rm=T))) %>%
				  filter(!(GRAIL_ID %in% c('MSK-VL-0057','MSK-VL-0003','MSK-VB-0046')))
				  
pattern = c("VB", "VL", "VP")
tissue = c("Breast", "Lung", "Prostate")

tmp.1 = p.1 = NULL
for (i in 1:length(pattern)) {
	tmp = volumetric_data %>%
		  mutate(tot_volu = ifelse(tot_volu==0, 0.1, tot_volu)) %>%
	 	  mutate(ln_frac = log(ctdna_frac)) %>%
	 	  mutate(ln_vol = log(tot_volu)) %>%
		  filter(grepl(pattern[i], GRAIL_ID)) %>%
		  mutate(Tissue = tissue[i]) %>%
		  filter(!(is.infinite(ln_frac) | is.na(ln_frac))) %>%
		  filter(!(is.infinite(ln_vol) | is.na(ln_vol)))
		
	if (tissue[i]=="Breast") {
		tmp = tmp %>%
			  mutate(shape = ifelse(is.na(Bone_Lesions) & is.na(Pleural_Disease), 21, 24))
	} else if (tissue[i]=="Lung") {
		tmp = tmp %>%
			  mutate(shape = ifelse(is.na(Bone_Lesions) & is.na(Pleural_Disease), 21, 24)) %>%
	  		  # Large left forearm bone and soft-tissue mass 10 by 7 by 17 cm not included in PET-scan and measurements
	  		  mutate(shape = ifelse(GRAIL_ID=="MSK-VL-0008", 24, shape))
	} else if (tissue[i]=="Prostate") {
		tmp = tmp %>%
			  mutate(shape = 24)
	}
	

	tmp.0 = tmp %>%
			mutate(cat = case_when(
				tot_volu >= 0 & tot_volu <= quantile(tot_volu, 1/3) ~ 1,
				tot_volu > quantile(tot_volu, 1/3) & tot_volu <= quantile(tot_volu, 2/3) ~ 2,
				tot_volu > quantile(tot_volu, 2/3) ~ 3,
			)) %>%
 			mutate(cat = as.factor(cat))
 	p = signif(jonckheere.test(x=tmp.0$ctdna_frac, g=as.numeric(tmp.0$cat), alternative = "increasing")$p.value, 3)
 	tmp.1 = rbind(tmp.1, tmp.0)
 	p.1 = c(p.1, p)
}

plot.0 = ggplot(tmp.1, aes(x = Tissue, y = ctdna_frac, color = cat, group=interaction(Tissue, cat))) + 
		 geom_boxplot(alpha=1, outlier.size=NA, outlier.shape=NA, fill="white", width=.8) +
		 scale_colour_manual(values=c("1"="black", "2"="black", "3"="black")) +
		 geom_jitter(
		 	aes(x = Tissue, y = ctdna_frac, fill = factor(Tissue), shape = factor(shape)),
  			position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
  			alpha=1, size=3
  			) +
  		 scale_shape_manual(values=c("21"=21, "24"=24)) +
  		 scale_fill_manual(values=cohort_cols) +
		 theme_classic(base_size=15) +
 		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10), legend.position="none") +
 		 labs(x="", y="ctDNA fraction\n") +
 		 coord_cartesian(ylim = c(0,1))

pdf(file="../res/figure2/cfdna_fraction_by_total_dx_vol.pdf", width=8.5, height=5)
print(plot.0)
dev.off()

export_x = tmp.1 %>%
		   dplyr::select(`patient_id` = `GRAIL_ID`,
		   				 `tissue` = `Tissue`,
		   				 `ctdna_fraction` = `ctdna_frac`,
		   				 `tertiles_disease_volume` = `cat`)
write_tsv(export_x, path="../res/etc/Source_Data_Fig_2/Fig_2g.tsv", append=FALSE, col_names=TRUE)

#==================================================
# trend test ctdna fraction by cancer type
#==================================================
x = cfdna_frac[,"ctdna_frac"]
y = rep(1, length(x))
y[grepl("VL", cfdna_frac[,1])] = 2
y[grepl("VP", cfdna_frac[,1])] = 3
 
pander(cor.test(x=x, y=y, method="kendall", alternative="two.sided", exact=FALSE)$p.value)
pander(kruskal.test(x=x, g=y)$p.value)


#==================================================
# trend test ctdna fraction versus number of sites
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

x = cfdna_frac[,"ctdna_frac"]
y = number_metastatic_sites

pander(cor.test(x=x, y=y, method="kendall", alternative="greater", exact=FALSE)$p.value)
pander(jonckheere.test(x=x, g=y, alternative="increasing", nperm=10000)$p.value)

# jonckheere-terpstra test by cancer type
for (i in c("VB", "VL", "VP")) {
  	index = grepl(i, names(number_metastatic_sites))
  	x = cfdna_frac[index,"ctdna_frac"]
  	y = number_metastatic_sites[index]
  	pander(jonckheere.test(x=x, g=y, alternative="increasing", nperm=10000)$p.value)
}
for (i in c("VB", "VL", "VP")) {
  	index = grepl(i, names(number_metastatic_sites))
  	x = cfdna_frac[index,"ctdna_frac"]
  	y = number_metastatic_sites[index]
  	y0 = rep(1, length(y))
  	y0[y==3] = 2
  	y0[y>3] = 3
 	pander(jonckheere.test(x=x, g=y0, alternative="increasing", nperm=10000)$p.value)
}
 
# kendall correlation test by cancer type
for (i in c("VB", "VL", "VP")) {
 	index = grepl(i, names(number_metastatic_sites))
 	x = cfdna_frac[index,"ctdna_frac"]
 	y = number_metastatic_sites[index]
  	pander(cor.test(x=x, y=y, method="kendall", alternative="greater", exact = FALSE)$p.value)
}
for (i in c("VB", "VL", "VP")) {
 	index = grepl(i, names(number_metastatic_sites))
 	x = cfdna_frac[index,"ctdna_frac"]
 	y = number_metastatic_sites[index]
 	y0 = rep(1, length(y))
 	y0[y==3] = 2
  	y0[y>3] = 3
 	pander(cor.test(x=x, y=y0, method="kendall", alternative="greater", exact = FALSE)$p.value)
}
 
# mblm association test by cancer type
for (i in c("VB", "VL", "VP")) {
 	index = grepl(i, names(number_metastatic_sites))
 	x = cfdna_frac[index,"ctdna_frac"]
 	y = number_metastatic_sites[index]
	data = data.frame(ctDNA_fraction=x, n_sites=y) 
 
 	z = mblm(ctDNA_fraction ~ n_sites, dataframe=data, repeated=FALSE)
 	pander(z)
}
  	
#==================================================
# p-values of detection rate
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
pander(p)


#==================================================
# cfDNA detection rate by tumor CCF
#==================================================
tumor_vars = read_tsv(file="../modified_v11/Resources/MSK_additional_data/20180405_table_1f_grail_paper_filters.called_in_both_and_impact.tsv",
			 col_types = cols(.default = col_character())) %>%
		     type_convert()
tumor_vars = tumor_vars %>%
			 mutate("tissue" = case_when(
			 	COHORT == "validation_breast" ~ "Breast",
			 	COHORT == "validation_lung" ~ "Lung",
			 	COHORT == "validation_prostate" ~ "Prostate"
			 )) %>%
			 mutate(is_snv = ifelse(nchar(ALT)==1 & nchar(REF)==1, TRUE, FALSE)) %>%
			 filter(is_snv)

cohort = c("Breast", "Lung", "Prostate")
p = rep(NA, length(cohort))
names(p) = cohort
dfm = list()
for (i in 1:length(cohort)) {
	tmp_vars = tumor_vars %>%
			   filter(tissue == cohort[i] & CALLED_IN %in% c("both","impact"))
	ccf = tmp_vars %>%
		  .[['IM-T.ccf']]
	ccf_cat = cut(ccf, c(0,0.5,0.8,1))
	called_in = tmp_vars %>%
		  	    .[['CALLED_IN']]
	dr_sum = by(called_in, ccf_cat, function(x) { sum(x=="both") })
	tt_sum = table(ccf_cat)
	dr = dr_sum/tt_sum
	ci = binconf(x=dr_sum, n=tt_sum)
	p[i] = prop.trend.test(dr_sum, tt_sum)$p.value

	x = c(sum(tumor_vars$`IM-T.clonality` %in% c("clonal","likely_clonal") & tumor_vars$CALLED_IN=="both"),
		  sum(tumor_vars$`IM-T.clonality` %in% c("clonal","likely_clonal") & tumor_vars$`IM-T.clonality`!=""),
		  sum(tumor_vars$`IM-T.clonality` %in% c("subclonal") & tumor_vars$CALLED_IN=="both"),
		  sum(tumor_vars$`IM-T.clonality` %in% c("subclonal") & tumor_vars$`IM-T.clonality`!=""))

	dfm[[i]] = data.frame(intervals = names(tt_sum),
					 	  dr = as.vector(dr),
					 	  lci = as.vector(ci[,2]),
					 	  uci = as.vector(ci[,3])) %>%
		  				  mutate(tissue = cohort[i])
}
dfm = do.call(rbind, dfm)

export_x = tumor_vars %>%
		   dplyr::select(
		   		`patient_id` = `CASE`,
		   		`tissue` = `tissue`,
		   		`gene_id` = `GENE`,
		   		`chromosome` = `CHROM`,
		   		`position` = `POS`,
		   		`called_in` = `CALLED_IN`,
		   		`t_ccf` = `IM-T.ccf`,
		   		`t_clonality` = `IM-T.clonality`) %>%
		   filter(called_in %in% c("both", "impact"))
write_tsv(export_x, path="../res/etc/Source_Data_Fig_2/Fig_2e.tsv", append=FALSE, col_names=TRUE)
