#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
source('config.R')

if (!dir.exists("../res/figureS5")) {
	dir.create("../res/figureS5")
}

#==================================================
# barplot of recurrent genes
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
 		   
# remove hypermutators
variants = variants %>%
		   filter(!(patient_id %in% hypermutators$patient_id) | (patient_id %in% msi_hypermutators$patient_id))
 		   
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
 						 filter(bio_source %in% c("biopsy_matched", "IMPACT-BAM_matched", "VUSo"), gene %in% gene_list) %>%
 						 dplyr::select(subj_type, gene, percent_patient, bio_source)
   
top_cancer_genes_table = top_cancer_genes_table %>%
 						 mutate(gene = factor(gene, levels = top_cancer_genes_ordered),
 						 subj_type = factor(subj_type, levels = c("Control", "Breast", "Lung", "Prostate")),
 						 bio_source = factor(bio_source, levels = c("VUSo", "biopsy_matched", "IMPACT-BAM_matched")))
  
bar_fill_cols = variant_cols[c("biopsy_matched", "IMPACT-BAM_matched", "VUSo")]

pdf(file="../res/figureS5/recurrent_genes_combined.pdf", width=10, height=8.7)
par(mar = c(6.1, 6, 4.1, 1))
zz = split.screen(figs=matrix(c(0,1, 3/4-.1, 1,
								0,1, 1/2-.05, 3/4+.05,
								0,1, 1/4-.05, 1/2+.05,
								0,1, 0,   1/4+.1,
								0,1,0,1),
				  nrow=5, ncol=4, byrow=TRUE))
cancer_types = c("Control", "Breast", "Lung", "Prostate")
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
 	if (i==4) {
 		axis(side=1, at=zzz, labels=top_cancer_genes_ordered, las=2, font=3, cex=.75, lwd=-1)
 	}
 	title(main=paste0("\n\n",cancer_types[i]), cex.main=1.25)
}
screen(zz[5])
plot(0, 0, type="n", xlim=c(0,10), ylim=c(0,10), xlab="", ylab="", axes=FALSE, frame.plot=FALSE)
mtext(side = 2, text = "% of patients", line = 3, cex = 1.15)
close.screen(all.screens=TRUE)
dev.off()

#==================================================
# number of VUSo per gene versus gene length
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

variants_hyper = variants %>%
 				 filter((patient_id %in% hypermutators$patient_id) | (patient_id %in% msi_hypermutators$patient_id)) %>%
 				 filter(bio_source=="VUSo", is_nonsyn)

unique_genes = na.omit(unique(variants_hyper$gene))
subj_types = c("Breast", "Lung", "Prostate")
n = matrix(NA, nrow=length(unique_genes), ncol=3, dimnames=list(unique_genes, subj_types))
for (i in 1:length(unique_genes)) {
	for (j in 1:length(subj_types)) {
		index = which(variants_hyper$gene==unique_genes[i] & variants_hyper$subj_type==subj_types[j])
		n[i,j] = length(index)
 	}
}
variants_by_gene = apply(n, 1, sum)
all_genes = names(variants_by_gene)
all_genes[all_genes=="FOXL2NB"] = "C3orf72"
names(variants_by_gene)[names(variants_by_gene)=="FOXL2NB"] = "C3orf72"
bed = read.csv(file=common_bed_annotated_w_introns, header=FALSE, sep="\t", stringsAsFactors=FALSE)
target_lengths = vector(mode="numeric", length=length(all_genes))
intron_sizes = vector(mode="numeric", length=length(all_genes))
for (i in 1:length(all_genes)) {
 	index = bed[,4]==all_genes[i]
 	target_lengths[i] = sum(bed[index,3]-bed[index,2])
 	intron_sizes[i] = t(as.matrix(bed[index,3]-bed[index,2])) %*% as.matrix(bed[index,5,drop=FALSE])
}
index = target_lengths == intron_sizes
intron_sizes[index] = 0
target_lengths = target_lengths - intron_sizes
index = order(target_lengths, decreasing=TRUE)
variants_by_gene = variants_by_gene[index]
target_lengths = target_lengths[index]
all_genes = all_genes[index]

tmp.x = data_frame(x = target_lengths,
				   y_0 = variants_by_gene,
				   y_1 = variants_by_gene/length(unique(variants_hyper$patient_id)), 
				   z_0 = all_genes) %>%
		mutate(z_1 = "+")
				   
variants_nohyper = variants %>%
 				 filter(!(patient_id %in% hypermutators$patient_id) | (patient_id %in% msi_hypermutators$patient_id)) %>%
 				 filter(bio_source=="VUSo", is_nonsyn)

unique_genes = na.omit(unique(variants_nohyper$gene))
subj_types = c("Breast", "Lung", "Prostate")
n = matrix(NA, nrow=length(unique_genes), ncol=3, dimnames=list(unique_genes, subj_types))
for (i in 1:length(unique_genes)) {
	for (j in 1:length(subj_types)) {
		index = which(variants_nohyper$gene==unique_genes[i] & variants_nohyper$subj_type==subj_types[j])
		n[i,j] = length(index)
 	}
}
variants_by_gene = apply(n, 1, sum)
all_genes = names(variants_by_gene)
all_genes[all_genes=="FOXL2NB"] = "C3orf72"
names(variants_by_gene)[names(variants_by_gene)=="FOXL2NB"] = "C3orf72"
bed = read.csv(file=common_bed_annotated_w_introns, header=FALSE, sep="\t", stringsAsFactors=FALSE)
target_lengths = vector(mode="numeric", length=length(all_genes))
intron_sizes = vector(mode="numeric", length=length(all_genes))
for (i in 1:length(all_genes)) {
 	index = bed[,4]==all_genes[i]
 	target_lengths[i] = sum(bed[index,3]-bed[index,2])
 	intron_sizes[i] = t(as.matrix(bed[index,3]-bed[index,2])) %*% as.matrix(bed[index,5,drop=FALSE])
}
index = target_lengths == intron_sizes
intron_sizes[index] = 0
target_lengths = target_lengths - intron_sizes
index = order(target_lengths, decreasing=TRUE)
variants_by_gene = variants_by_gene[index]
target_lengths = target_lengths[index]
all_genes = all_genes[index]

tmp.y = data_frame(x = target_lengths,
				   y_0 = variants_by_gene,
				   y_1 = variants_by_gene/length(unique(variants_nohyper$patient_id)), 
				   z_0 = all_genes) %>%
		mutate(z_1 = "-")
				   
tmp.0 = bind_rows(tmp.x, tmp.y)

plot.0 = ggplot(tmp.0, aes(x=x/10e3, y=y_1, fill=z_1)) +
		 geom_point(alpha=1, size=2.5, shape=21) +
		 geom_smooth(method="lm", se=TRUE, level=0.95) +
		 labs(x="\nCoding target length (Kb)\n", y="Number of variants / patient\n") +
		 theme_classic(base_size=15) +
		 scale_x_log10() +
		 ylim(0, 2.5)

pdf(file="../res/figureS5/VUSo_target_lengths_hyper.pdf", width=7, height=7)
par(mar = c(6.1, 7, 4.1, 1))
shapes = 1
plot(1, 1, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, xlim=c(.001,2), ylim=c(0.005,2.5), log="x")
axis(1, at=c(1e-3, 1e-2, 1e-1, 1, 2), labels=c(expression(10^-3), expression(10^-2), expression(10^-1), 1, 2), cex.axis = 1.75, padj = 0.25, lwd=1.85, lwd.ticks=1.75)
axis(2, at=NULL, labels=NULL, cex.axis = 1.75, las = 1, lwd=1.85, lwd.ticks=1.75)
mtext(side = 1, text = "Coding target length (Kb)", line = 4, cex = 1.75)
mtext(side = 2, text = "VUSo / patient", line = 4.5, cex = 1.75)
x = as.numeric(tmp.0$x[tmp.0$z_1=="+"])/10e3
y = as.numeric(tmp.0$y_1[tmp.0$z_1=="+"])
points(x, y, pch=1, col=transparent_rgb("salmon", 205), cex=1.5, lwd=1.5)
fit = lm(y ~ x^2)
y_p = predict(object = fit, nexdata=data.frame(x=seq(.01, 10, l=1000)))
points(x, y_p, type="l", col="goldenrod3", lwd=3, lty=1)
x = as.numeric(tmp.0$x[tmp.0$z_1=="-"])/10e3
y = as.numeric(tmp.0$y_1[tmp.0$z_1=="-"])
points(x, y, pch=1, col=transparent_rgb("steelblue", 205), cex=1.5, lwd=1.5)
fit = lm(y ~ x)
y_p = predict(object = fit, nexdata=data.frame(x=seq(.01, 10, l=1000)))
points(x, y_p, type="l", col="goldenrod3", lwd=3, lty=1)
log10_axis(side=1, at=c(1e-3, 1e-2, 1e-1, 1), lwd=0, lwd.ticks=1)
legend(x=1e-3, y=2.5, pch=shapes, col=c("salmon", "steelblue"), legend=c("Hypermutated", "Non-hypermutated"), box.lwd=-1, pt.cex=1.35, pt.lwd=1.5)
dev.off()

pdf(file="../res/figureS5/VUSo_target_lengths_nohyper.pdf", width=7, height=7)
par(mar = c(6.1, 7, 4.1, 1))
shapes = 1
plot(1, 1, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, xlim=c(.001,2), ylim=c(0.005,.1), log="xy")
axis(1, at=c(1e-3, 1e-2, 1e-1, 1, 2), labels=c(expression(10^-3), expression(10^-2), expression(10^-1), 1, 2), cex.axis = 1.75, padj = 0.25, lwd=1.85, lwd.ticks=1.75)
axis(2, at=c(5e-3, 1e-2, 2e-2, 5e-2, 1e-1), labels=c(expression("5e"^-3), expression("1e"^-2), expression("2e"^-2), expression("5e"^-2), 0.1), cex.axis = 1.75, las = 1, lwd=1.85, lwd.ticks=1.75)
mtext(side = 1, text = "Coding target length (Kb)", line = 4, cex = 1.75)
mtext(side = 2, text = "VUSo / patient", line = 4.5, cex = 1.75)
x = as.numeric(tmp.0$x[tmp.0$z_1=="-"])/10e3
y = as.numeric(tmp.0$y_1[tmp.0$z_1=="-"])
points(x, y, pch=1, col=transparent_rgb("steelblue", 205), cex=1.5, lwd=1.5)
fit = lm(y ~ x)
y_p = predict(object = fit, nexdata=data.frame(x=seq(.01, 10, l=1000)))
points(x, y_p, type="l", col="goldenrod3", lwd=3, lty=1)
log10_axis(side=1, at=c(1e-3, 1e-2, 1e-1, 1), lwd=0, lwd.ticks=1)
log10_axis(side=2, at=c(5e-3, 1e-2, 2e-2, 5e-2, 1e-1), lwd=0, lwd.ticks=1)
dev.off()

