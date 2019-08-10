#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/figureS8")) {
	dir.create("../res/figureS8")
}

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

pdf(file="../res/figureS8/barplot_bio_sources_hyper.pdf", width=12, height=9)
par(mar = c(6.1, 6, 4.1, 1))
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
  	
  	zzz = barplot(wbc_matched[,"num"]-.9, col=variant_cols["WBC_matched"], border="black", space=.16, axes=FALSE, ylim=c(-30,150), lwd=.01)
  	abline(h=0, col="white", lwd=4)
  	barplot(ifelse((vuso[,"num"]+biopsy_matched[,"num"]+biopsy_subthreshold[,"num"])<110, (vuso[,"num"]+biopsy_matched[,"num"]+biopsy_subthreshold[,"num"]), 120), col=variant_cols["VUSo"], border="black", space=.16, add=TRUE, axes=FALSE, lwd=.01)
  	barplot(biopsy_matched[,"num"]+biopsy_subthreshold[,"num"], col=variant_cols["IMPACT-BAM_matched"], border="black", space=.16, add=TRUE, axes=FALSE, lwd=.01)
  	barplot(biopsy_matched[,"num"], col=variant_cols["biopsy_matched"], border="black", space=.16, add=TRUE, axes=FALSE, lwd=.01)
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
		rect(xleft=zzz[36]-.565, xright=zzz[36]+.565, ybottom=130, ytop=130+(587-550)-20, col=variant_cols["VUSo"], border="black")
	}
	title(main=subj[i])
}
screen(zz[5])
plot(0,0, type="n", xlim=c(0,10), ylim=c(-1,1), xlab="", ylab="", axes=FALSE, frame.plot=FALSE)
mtext(side = 2, text = expression("Somatic cfDNA variants / Mb"), line = 1.5, cex = 1.4)
close.screen(all.screens=TRUE)
dev.off()

#==================================================
# heatmap of wbc-matched mutations with hypermutators
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

variants = left_join(variants, msk_anno %>% dplyr::select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
variants = variants %>%
		   mutate(bio_source = case_when(
		   					   MSK == 1 & grail == 1 ~ "biopsy_matched",
		   					   MSK == 1 & grail == 0 ~ "biopsy_only",
		   					   category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
		   					   category %in% c("germline", "germlineish") ~ "germline",
		   					   category %in% c("blood", "bloodier") ~ "WBC_matched",
		   					   category == "somatic" ~ "VUSo",
		   					   TRUE ~ "other"),
		   		  af = ifelse(is.na(af), 0, af),
		   		  af_nobaq = round(adnobaq / dpnobaq * 100, 2),
		   		  af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq))
		   		  
variants = variants %>%
 		   filter(bio_source=="WBC_matched", is_nonsyn)
unique_genes = na.omit(unique(variants$gene))
subj_types = c("Control", "Breast", "Lung", "Prostate")
n = matrix(NA, nrow=length(unique_genes), ncol=4, dimnames=list(unique_genes, subj_types))
for (i in 1:length(unique_genes)) {
	for (j in 1:length(subj_types)) {
		index = which(variants$gene==unique_genes[i] & variants$subj_type==subj_types[j])
		n[i,j] = length(index)
	}
}
index = order(apply(n, 1, sum), decreasing=TRUE)
n = n[index,,drop=FALSE]
pdf(file="../res/figureS8/heatmap_wbc_matched_by_variants.pdf", width=5, height=21)
corrplot(corr=n[1:40,,drop=FALSE], method="color", type="full", add=FALSE,
		 col=colorRampPalette(c(rep("#ffffff",2), "#19547b"))(100), bg=NULL, is.corr=FALSE,
		 addgrid.col="white", diag=TRUE, outline=FALSE, mar=c(1, 0, 1, 0), cl.pos="n",
 		 addCoef.col = "black", number.font = 1, number.cex = 1.15, tl.cex = 1.5, tl.col = "black", tl.offset = 0.5)
dev.off()

subj_num_smry = variants %>%
				group_by(subj_type) %>%
				summarise(subj_num = length(unique(patient_id))) %>%
				ungroup() %>%
				mutate(subj_type_num = paste0(subj_type, "(", subj_num, ")"))

gene_recurrences = variants %>%
 				   filter(!is.na(gene)) %>%
				   group_by(bio_source, subj_type, gene, is_nonsyn) %>%
 				   summarize(number = n(), num_patient = length(unique(patient_id))) %>%
 				   ungroup %>%
 				   left_join(subj_num_smry) %>%
 				   mutate(percent_patient = round(num_patient / subj_num * 100, 0))
   
gene_list = c()
gene_stats_table_list = list()

for (subj in subj_num_smry$subj_type) {
   top_wbc_gene_ids = gene_recurrences %>%
  					  filter(subj_type == subj, bio_source == "WBC_matched", is_nonsyn) %>%
   					  group_by(gene) %>%
    				  summarise(all = sum(percent_patient)) %>%
   					  ungroup() %>%
   					  arrange(-all) %>%
   					  dplyr::slice(1:20) %>%
   					  .[["gene"]]
   
   gene_list = c(gene_list, top_wbc_gene_ids)
}
 
for (subj in c("Control", "Breast", "Lung", "Prostate")) {
  subj_smry = gene_recurrences %>%
   			  filter(subj_type == subj, bio_source == "WBC_matched", is_nonsyn, gene %in% unique(gene_list)) %>%
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

gene_stats_matrix = as.matrix(gene_stats_smry %>%
					dplyr::select(-gene, -Mean))
rownames(gene_stats_matrix) = gene_stats_smry$gene
idx_names = match(colnames(gene_stats_matrix), subj_num_smry$subj_type)
colnames(gene_stats_matrix) = subj_num_smry$subj_type_num[idx_names]

max_num = max(gene_stats_matrix) + 3
my_palette = colorRampPalette(c("white","#1F6784"))(n = 29)
col_breaks = c(seq(0, 5, length = 10), seq(5.1, max_num, length = 20))
colnames(gene_stats_matrix) = unlist(lapply(strsplit(subj_num_smry$subj_type_num[idx_names], "(", fixed=TRUE), function(x) {x[1]}))

pdf(file="../res/figureS8/heatmap_wbc_matched_by_gene.pdf", width=5, height=21)
corrplot(corr=gene_stats_matrix, method="color", type="full", add=FALSE,
		 col=colorRampPalette(c(rep("#ffffff",2), "#19547b"))(100), bg=NULL, is.corr=FALSE,
		 addgrid.col="white", diag=TRUE, outline=FALSE, mar=c(1, 0, 1, 0), cl.pos="n",
 		 addCoef.col = "black", number.font = 1, number.cex = 1.15, tl.cex = 1.5, tl.col = "black", tl.offset = 0.5)
dev.off()

#==================================================
# scatter plot of vaf in wbc and cfdna
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
						 mutate(facets_2 = "VAF with pseudocounts without BAQ") %>%
						 mutate(facets_3 = "VAF without pseudocounts without BAQ") %>%
						 mutate(bio_source = case_when(
						 	bio_source == "biopsy_matched" ~ "Biopsy matched",
						 	bio_source == "IMPACT-BAM_matched" ~ "Biopsy subthreshold",
						 	bio_source == "VUSo" ~ "VUSo",
						 	bio_source == "WBC_matched" ~ "WBC matched",
						 ))

cols = variant_cols
names(cols) = c("VUSo", "Biopsy matched", "Biopsy subthreshold", "WBC matched", "Biopsy only")
         
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
			 
pdf(file="../res/figureS8/VAF_VAF_posterior_mean.pdf", width=7, height=7)
par(mar = c(6.1, 6, 4.1, 1))
plot(1,1, type="n", xlim=c(0.01,15), ylim=c(0.01,15), xlab="", ylab="", axes=FALSE, frame.plot=FALSE, log="xy")
points(x=c(0.01, 15), y=c(0.01,15), type="l", lty=1, lwd=2, col="goldenrod3")
x = cfdna_vs_gdna_vaf_plot$afmeancfdna
y = cfdna_vs_gdna_vaf_plot$afmeangdna
x[x>15] = NA
y[y>15] = NA
for (i in c("VUSo", "WBC matched", "Biopsy matched", "Biopsy subthreshold")) {
	index = cfdna_vs_gdna_vaf_plot$bio_source==i
	points(x[index], y[index], type="p", col="black", bg=cols[i], pch=21, lwd=.9, cex=1.15)
}
axis(1, at = c(.01, .10, 1.0, 10.0, 15.0), labels = c(expression(10^-2), expression(10^-1), "1", "10", ""), cex.axis = 1.5, padj = 0.25, lwd=1.75, lwd.ticks=1.65)
axis(1, at = 15, labels = "15", cex.axis = 1.5, padj = 0.25, lwd=-1, las=1)
axis(2, at = c(.01, .10, 1.0, 10.0, 15.0), labels = c(expression(10^-2), expression(10^-1), "1", "10", ""), cex.axis = 1.5, las = 1, lwd=1.75, lwd.ticks=1.65)
axis(2, at = 15, labels = "15", cex.axis = 1.5, padj = 0.25, lwd=-1, las=1)
log10_axis(side=1, at=c(.01, .1, 1, 10), lwd=0, lwd.ticks=1)
log10_axis(side=2, at=c(.01, .1, 1, 10), lwd=0, lwd.ticks=1)
mtext(side = 1, text = "VAF in cfDNA (%)", line = 4, cex = 1.85)
mtext(side = 2, text = "VAF in WBC (%)", line = 4, cex = 1.85)
legend(x=0.009, y=19, pch=21, col="black", pt.bg=cols, pt.cex=1.55, pt.lwd=.5, legend=c("Biopsy matched", "Biopsy subthreshold", "VUSo", "WBC matched"), box.lwd=-1, cex=1.15)
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
		 
pdf(file="../res/figureS8/VAF_VAF_pseudo_no_BAQ.pdf", width=7, height=7)
par(mar = c(6.1, 6, 4.1, 1))
plot(1,1, type="n", xlim=c(0.01,15), ylim=c(0.01,15), xlab="", ylab="", axes=FALSE, frame.plot=FALSE, log="xy")
points(x=c(0.01, 15), y=c(0.01,15), type="l", lty=1, lwd=2, col="goldenrod3")
x = cfdna_vs_gdna_vaf_plot$afcfdna_nobaq
y = cfdna_vs_gdna_vaf_plot$afgdna_nobaq
x[x>15] = NA
y[y>15] = NA
for (i in c("VUSo", "WBC matched", "Biopsy matched", "Biopsy subthreshold")) {
	index = cfdna_vs_gdna_vaf_plot$bio_source==i
	points(x[index], y[index], type="p", col="black", bg=cols[i], pch=21, lwd=.9, cex=1.15)
}
axis(1, at = c(.01, .10, 1.0, 10.0, 15.0), labels = c(expression(10^-2), expression(10^-1), "1", "10", ""), cex.axis = 1.5, padj = 0.25, lwd=1.75, lwd.ticks=1.65)
axis(1, at = 15, labels = "15", cex.axis = 1.5, padj = 0.25, lwd=-1, las=1)
axis(2, at = c(.01, .10, 1.0, 10.0, 15.0), labels = c(expression(10^-2), expression(10^-1), "1", "10", ""), cex.axis = 1.5, las = 1, lwd=1.75, lwd.ticks=1.65)
axis(2, at = 15, labels = "15", cex.axis = 1.5, padj = 0.25, lwd=-1, las=1)
log10_axis(side=1, at=c(.01, .1, 1, 10), lwd=0, lwd.ticks=1)
log10_axis(side=2, at=c(.01, .1, 1, 10), lwd=0, lwd.ticks=1)
mtext(side = 1, text = "VAF in cfDNA (%)", line = 4, cex = 1.85)
mtext(side = 2, text = "VAF in WBC (%)", line = 4, cex = 1.85)
legend(x=0.009, y=19, pch=21, col="black", pt.bg=cols, pt.cex=1.55, pt.lwd=.5, legend=c("Biopsy matched", "Biopsy subthreshold", "VUSo", "WBC matched"), box.lwd=-1, cex=1.15)
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
  			 	labels = function(x) { c("0.01", "0.1", "1", "10", "100") }
  			 ) + 
  			 scale_y_log10(
  			 	breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
  			 	labels = function(x) { c("0.01", "0.1", "1", "10", "100") }
  			 ) +
 			 coord_cartesian(xlim = c(0.01, 100), ylim = c(0.01,100)) +
 			 annotation_logticks() +
			 guides(fill=guide_legend(title=c("Variant category")))
		 
pdf(file="../res/figureS8/VAF_VAF_nopsedo_no_BAQ.pdf", width=7, height=7)
par(mar = c(6.1, 6, 4.1, 1))
plot(1,1, type="n", xlim=c(0.01,15), ylim=c(0.01,15), xlab="", ylab="", axes=FALSE, frame.plot=FALSE, log="xy")
points(x=c(0.01, 15), y=c(0.01,15), type="l", lty=1, lwd=2, col="goldenrod3")
x = cfdna_vs_gdna_vaf_plot$afcfdna_nobaq_nos
y = cfdna_vs_gdna_vaf_plot$afgdna_nobaq_nos
x[x>15] = NA
y[y>15] = NA
for (i in c("VUSo", "WBC matched", "Biopsy matched", "Biopsy subthreshold")) {
	index = cfdna_vs_gdna_vaf_plot$bio_source==i
	points(x[index], y[index], type="p", col="black", bg=cols[i], pch=21, lwd=.9, cex=1.15)
}
axis(1, at = c(.01, .10, 1.0, 10.0, 15.0), labels = c(expression(10^-2), expression(10^-1), "1", "10", ""), cex.axis = 1.5, padj = 0.25, lwd=1.75, lwd.ticks=1.65)
axis(1, at = 15, labels = "15", cex.axis = 1.5, padj = 0.25, lwd=-1, las=1)
axis(2, at = c(.01, .10, 1.0, 10.0, 15.0), labels = c(expression(10^-2), expression(10^-1), "1", "10", ""), cex.axis = 1.5, las = 1, lwd=1.75, lwd.ticks=1.65)
axis(2, at = 15, labels = "15", cex.axis = 1.5, padj = 0.25, lwd=-1, las=1)
log10_axis(side=1, at=c(.01, .1, 1, 10), lwd=0, lwd.ticks=1)
log10_axis(side=2, at=c(.01, .1, 1, 10), lwd=0, lwd.ticks=1)
mtext(side = 1, text = "VAF in cfDNA (%)", line = 4, cex = 1.85)
mtext(side = 2, text = "VAF in WBC (%)", line = 4, cex = 1.85)
legend(x=0.009, y=19, pch=21, col="black", pt.bg=cols, pt.cex=1.55, pt.lwd=.5, legend=c("Biopsy matched", "Biopsy subthreshold", "VUSo", "WBC matched"), box.lwd=-1, cex=1.15)
dev.off()
