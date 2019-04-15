#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
source('config.R')

if (!dir.exists("../res/figureS9")) {
	dir.create("../res/figureS9")
}

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

pdf(file="../res/figureS9/target_lengths_number_variants.pdf", width=8, height=8)
par(mar = c(6.1, 6, 4.1, 1))
plot(target_lengths, variants_by_gene, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, ylim=c(1,30), xlim=c(0,20000))
points(target_lengths, variants_by_gene, pch=21, col="black", bg="salmon", cex=1.35)
axis(1, at = NULL, labels = NULL, cex.axis = 1.5, las = 1, lwd=1.5, lwd.ticks=1.35)
axis(2, at = c(1,seq(5,30,by=5)), labels=c(1,seq(5,30,by=5)), cex.axis = 1.5, las = 1, lwd=1.5, lwd.ticks=1.35)
mtext(side = 1, text = "Coding target lengths", line = 4, cex = 1.85)
mtext(side = 2, text = "Number of variants", line = 4, cex = 1.85)
index = variants_by_gene>10
text(target_lengths[index], variants_by_gene[index], labels=names(variants_by_gene[index]), pos=3, cex=.95, font=3)
text(2000, 28.5, labels="p = 4.44e-16", pos=3, cex=1.25)
dev.off()
