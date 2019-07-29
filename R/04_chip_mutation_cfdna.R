#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/figure4")) {
	dir.create("../res/figure4")
}

#==================================================
# pie-charts of the number of variants classified
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

variants_control = variants %>%
		   filter(subj_type=="Control") %>%
		   filter(is_nonsyn | bio_source=="biopsy_only") %>%
		   filter(bio_source %in% c("VUSo","biopsy_matched", "biopsy_only", "IMPACT-BAM_matched", "WBC_matched"))

variants_nohyper = variants %>%
		   filter(subj_type!="Control") %>%
		   filter(!patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id)) %>%
		   filter(is_nonsyn | bio_source=="biopsy_only") %>%
		   filter(bio_source %in% c("VUSo","biopsy_matched", "biopsy_only", "IMPACT-BAM_matched", "WBC_matched"))
 
variants_hyper = variants %>%
		   filter(subj_type!="Control") %>%
		   filter(patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id)) %>%
		   filter(is_nonsyn | bio_source=="biopsy_only") %>%
		   filter(bio_source %in% c("VUSo","biopsy_matched", "biopsy_only", "IMPACT-BAM_matched", "WBC_matched"))
 		   
Mx = matrix(NA, nrow=5, ncol=5, dimnames=list(c("biopsy_only", "biopsy_matched", "IMPACT-BAM_matched", "WBC_matched", "VUSo"), c("Control", "Non-Hypermutators", "Breast", "Lung", "Prostate")))
for (i in 1:nrow(Mx)) {
	Mx[i,1] = sum(variants_control$bio_source==rownames(Mx)[i])
	Mx[i,2] = sum(variants_nohyper$bio_source==rownames(Mx)[i])
	Mx[i,3] = sum(variants_nohyper$bio_source==rownames(Mx)[i] & variants_nohyper$subj_type=="Breast")
	Mx[i,4] = sum(variants_nohyper$bio_source==rownames(Mx)[i] & variants_nohyper$subj_type=="Lung")
	Mx[i,5] = sum(variants_nohyper$bio_source==rownames(Mx)[i] & variants_nohyper$subj_type=="Prostate")
}	

My = matrix(NA, nrow=5, ncol=5, dimnames=list(c("biopsy_only", "biopsy_matched", "IMPACT-BAM_matched", "WBC_matched", "VUSo"), c("Control", "Hypermutators", "Breast", "Lung", "Prostate")))
for (i in 1:nrow(My)) {
	My[i,1] = 0
	My[i,2] = sum(variants_hyper$bio_source==rownames(My)[i])
	My[i,3] = sum(variants_hyper$bio_source==rownames(My)[i] & variants_hyper$subj_type=="Breast")
	My[i,4] = sum(variants_hyper$bio_source==rownames(My)[i] & variants_hyper$subj_type=="Lung")
	My[i,5] = sum(variants_hyper$bio_source==rownames(My)[i] & variants_hyper$subj_type=="Prostate")
}
Mx = t(Mx)
My = t(My)

Mx[2,5] = Mx[2,5]-1
Mx[2,2] = Mx[2,2]+1
My[2,5] = My[2,5]-1
My[2,2] = My[2,2]+1
 
pdf(file="../res/figure4/pie_bio_sources.pdf", height=7, width=15)
par(mar = c(6.1, 6.5, 4.1, 1.1))
z = split.screen(figs=matrix(c(0, 0.33, 0, 1, .33, .67, 0, 1, .67, 1, 0, 1), nrow=3, ncol=4, byrow=TRUE))
screen(z[1])
pie(x=Mx[1,4:5], col=variant_cols[colnames(Mx)[4:5]], labels=Mx[1,4:5])
screen(z[2])
pie(x=Mx[2,2:5], col=variant_cols[colnames(Mx)[2:5]], labels=Mx[2,2:5])
screen(z[3])
pie(x=My[2,2:5], col=variant_cols[colnames(My)[2:5]], labels=My[2,2:5])
close.screen(all.screens=TRUE)
dev.off()

#==================================================
# barplot of mutation burden and sources of mutation
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

variants = variants %>%
		   filter(!(patient_id %in% hypermutators$patient_id)) %>%
		   filter(!(patient_id %in% msi_hypermutators$patient_id))

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

pdf(file="../res/figure4/barplot_bio_sources.pdf", width=12, height=9)
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
  	
  	barplot(wbc_matched[,"num"]-.9, col=variant_cols["WBC_matched"], border="black", space=.16, axes=FALSE, ylim=c(-30,30), lwd=.01)
  	abline(h=0, col="white", lwd=4)
  	barplot(vuso[,"num"]+biopsy_matched[,"num"]+biopsy_subthreshold[,"num"], col=variant_cols["VUSo"], border="black", space=.16, add=TRUE, axes=FALSE, lwd=.01)
  	barplot(biopsy_matched[,"num"]+biopsy_subthreshold[,"num"], col=variant_cols["IMPACT-BAM_matched"], border="black", space=.16, add=TRUE, axes=FALSE, lwd=.01)
  	barplot(biopsy_matched[,"num"], col=variant_cols["biopsy_matched"], border="black", space=.16, add=TRUE, axes=FALSE, lwd=.01)
  	if (i==1 | i==3) {
  		axis(2, at = c(-30,-20,-10,0,10,20,30), labels=c(30,20,10,0,10,20,30), cex.axis = 1.4, las = 1)
  	} else {
		axis(2, at = c(-30,-20,-10,0,10,20,30), labels=rep("",7), cex.axis = 1.4, las = 1)
	}
	  title(main=subj[i])
}
screen(zz[5])
plot(0,0, type="n", xlim=c(0,10), ylim=c(-1,1), xlab="", ylab="", axes=FALSE, frame.plot=FALSE)
mtext(side = 2, text = expression("Somatic cfDNA variants / Mb"), line = 1.5, cex = 1.4)
close.screen(all.screens=TRUE)
dev.off()

#==================================================
# scatter plot of mutation burden versus age for
# biopsy matched, vuso and wbc-matched mutations
#==================================================
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
pdf(file="../res/figure4/mutational_burden_cfdna_combined.pdf", width=12, height=12)
par(mar = c(6.1, 6.5, 4.1, 1.1))
zz = split.screen(figs=matrix(c(0,.5,.5,1, .5,1,.5,1, 0,.5,0,.5, .5,1,0,.5), nrow=4, ncol=4, byrow=TRUE))

screen(zz[1])
plot(1, 1, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, ylim=c(1,50), xlim=c(20,90), log="y")
burden = bind_rows(burden_cancer, burden_healthy)
data = data.frame(burden[,c("patient_id", "age", "num_called", "subj_type", "bio_source"),drop=FALSE]) %>%
 	   filter(bio_source=="WBC_matched")
indx = unique(burden[!(burden$patient_id %in% data[,"patient_id"]),"patient_id",drop=TRUE])
if (length(indx)!=0) {
	tmp = matrix(NA, nrow=length(indx), ncol=5)
	tmp[,1] = indx
	tmp[,3] = 0
	tmp[,5] = "WBC_matched"
	zzz = data.frame(clinical[clinical$patient_id %in% indx,c("patient_id","age","subj_type")])
	rownames(zzz) = zzz[,"patient_id"]
	tmp[,2] = zzz[indx,"age"]
	tmp[,4] = zzz[indx,"subj_type"]
	colnames(tmp) = colnames(data)
	tmp = data.frame(tmp)
	data = rbind(data, tmp)
}
p0 = zeroinfl(as.numeric(num_called) ~ as.numeric(age), dist = "poisson", data = data)
text(42, 49, cex = 1.25, labels = paste("(P = ", toupper(signif(summary(p0)$coefficients$count[2,4], 3)), ")", sep = ""), pos = 4)
y0 = fun_zerob(x=as.numeric(data[,"age"]), y=as.numeric(data[,"num_called"]), n=1000, seed=0)
x1 = seq(20, 90, l=100)
y1 = t(apply(y0, 1, quantile, probs=c(.025, .5, .975), na.rm=TRUE))
polygon(c(x1, rev(x1)), c(y1[,"97.5%"], rev(y1[,"2.5%"])), col=transparent_rgb("grey10", 55), border=NA)
points(x1, y1[,"50%"], type="l", col="grey50", lwd=3)
data = data.frame(burden_healthy[,c("age", "num_called", "study", "bio_source"),drop=FALSE]) %>%
  	   filter(bio_source=="WBC_matched")
points(as.numeric(data$age), as.numeric(data$num_called), bg="#1F6784", col="#231F20", pch=21, cex=1.35, lwd=.5)
data = data.frame(burden_cancer[,c("age", "num_called", "study", "bio_source"),drop=FALSE]) %>%
 	   filter(bio_source=="WBC_matched")
points(as.numeric(data$age), as.numeric(data$num_called), bg="#EAA411", col="#231F20", pch=21, cex=1.35, lwd=.5)
axis(1, at = NULL, cex.axis = 1.5, padj = 0.25, lwd=1.5, lwd.ticks=1.35)
axis(1, at = seq(30, 90, by=20), labels=rep("", 4), padj = 0.25, lwd=-1, lwd.ticks=.85, tcl=-.35)
axis(2, at = c(1,2,5,10,20,30,50), labels=c(1,2,5,10,20,30,50), cex.axis = 1.85, las = 1, lwd=1.5, lwd.ticks=1.35)
mtext(side = 2, text = expression("Somatic cfDNA variants / Mb"), line = 5, cex = 1.85)
mtext(side = 1, text = "Age (years)", line = 4, cex = 1.5)
title(main="WBC matched", cex.main=1.5)
legend("topleft", legend=c("Cancer", "Control"), pch=21, col="#231F20", pt.bg=c("#EAA411", "#1F6784"), box.lwd=-1, pt.lwd=0, pt.cex=1.35, cex=1.15)

screen(zz[2])
plot(1, 1, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, ylim=c(1,50), xlim=c(20,90), log="y")
burden = bind_rows(burden_cancer, burden_healthy)
data = data.frame(burden[,c("patient_id", "age", "num_called", "bio_source", "subj_type"),drop=FALSE]) %>%
 	   filter(bio_source=="VUSo")
indx = unique(burden[!(burden$patient_id %in% data[,"patient_id"]),"patient_id",drop=TRUE])
if (length(indx)!=0) {
	tmp = matrix(NA, nrow=length(indx), ncol=5)
	tmp[,1] = indx
	tmp[,3] = 0
	tmp[,4] = "VUSo"
	zzz = data.frame(clinical[clinical$patient_id %in% indx,c("patient_id","age","subj_type")])
	rownames(zzz) = zzz[,"patient_id"]
	tmp[,2] = zzz[indx,"age"]
	tmp[,5] = zzz[indx,"subj_type"]
	colnames(tmp) = colnames(data)
	tmp = data.frame(tmp)
	data = rbind(data, tmp)
}
p0 = zeroinfl(as.numeric(num_called) ~ as.numeric(age), dist = "poisson", data = data)
text(42, 49, cex = 1.25, labels = paste("(P = ", toupper(signif(summary(p0)$coefficients$count[2,4], 3)), ")", sep = ""), pos = 4)
y0 = fun_zerob(x=as.numeric(data[,"age"]), y=as.numeric(data[,"num_called"]), n=1000, seed=0)
x1 = seq(20, 90, l=100)
y1 = t(apply(y0, 1, quantile, probs=c(.025, .5, .975), na.rm=TRUE))
polygon(c(x1, rev(x1)), c(y1[,"97.5%"], rev(y1[,"2.5%"])), col=transparent_rgb("grey10", 55), border=NA)
points(x1, y1[,"50%"], type="l", col="grey50", lwd=3)
data = data.frame(burden_healthy[,c("age", "num_called", "study", "bio_source"),drop=FALSE]) %>%
 	   filter(bio_source=="VUSo")
points(as.numeric(data$age), as.numeric(data$num_called), bg="#1F6784", col="#231F20", pch=21, cex=1.35, lwd=.5)
data = data.frame(burden_cancer[,c("age", "num_called", "study", "bio_source"),drop=FALSE]) %>%
 	   filter(bio_source=="VUSo")
points(as.numeric(data$age), as.numeric(data$num_called), bg="#EAA411", col="#231F20", pch=21, cex=1.35, lwd=.5)
axis(1, at = NULL, cex.axis = 1.5, padj = 0.25, lwd=1.5, lwd.ticks=1.35)
axis(1, at = seq(30, 90, by=20), labels=rep("", 4), padj = 0.25, lwd=-1, lwd.ticks=.85, tcl=-.35)
axis(2, at = c(1,2,5,10,20,30,50), labels=c(1,2,5,10,20,30,50), cex.axis = 1.85, las = 1, lwd=1.5, lwd.ticks=1.35)
mtext(side = 2, text = expression("Somatic cfDNA variants / Mb"), line = 5, cex = 1.85)
mtext(side = 1, text = "Age (years)", line = 4, cex = 1.5)
title(main="VUSo", cex.main=1.5)
legend("topleft", legend=c("Cancer", "Control"), pch=21, col="#231F20", pt.bg=c("#EAA411", "#1F6784"), box.lwd=-1, pt.lwd=0, pt.cex=1.35, cex=1.15)

screen(zz[3])
plot(1, 1, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, ylim=c(1,50), xlim=c(20,90), log="y")
data = data.frame(burden_cancer[,c("patient_id", "age", "num_called", "bio_source", "subj_type"),drop=FALSE]) %>%
 	   filter(bio_source=="biopsy_matched")
indx = unique(burden_cancer[!(burden_cancer$patient_id[burden_cancer$subj_type!="Control"] %in% data[,"patient_id"]),"patient_id",drop=TRUE])
if (length(indx)!=0) {
	tmp = matrix(NA, nrow=length(indx), ncol=5)
	tmp[,1] = indx
	tmp[,3] = 0
	tmp[,4] = "biopsy_matched"
	zzz = data.frame(clinical[clinical$patient_id %in% indx,c("patient_id","age","subj_type")])
	rownames(zzz) = zzz[,"patient_id"]
	tmp[,2] = zzz[indx,"age"]
	tmp[,5] = zzz[indx,"subj_type"]
	colnames(tmp) = colnames(data)
	tmp = data.frame(tmp)
	data = rbind(data, tmp)
}
p0 = zeroinfl(as.numeric(num_called) ~ as.numeric(age), dist = "poisson", data = data)
text(45, 49, cex = 1.25, labels = paste("(P = ", toupper(signif(summary(p0)$coefficients$count[2,4], 3)), ")", sep = ""), pos = 4)
y0 = fun_zerob(x=as.numeric(data[,"age"]), y=as.numeric(data[,"num_called"]), n=1000, seed=0)
x1 = seq(20, 90, l=100)
y1 = t(apply(y0, 1, quantile, probs=c(.025, .5, .975), na.rm=TRUE))
polygon(c(x1, rev(x1)), c(y1[,"97.5%"], rev(y1[,"2.5%"])), col=transparent_rgb("grey10", 55), border=NA)
points(x1, y1[,"50%"], type="l", col="grey50", lwd=3)
points(as.numeric(data$age), as.numeric(data$num_called), bg="#EAA411", col="#231F20", pch=21, cex=1.35, lwd=0.5)
axis(1, at = NULL, cex.axis = 1.5, padj = 0.25, lwd=1.5, lwd.ticks=1.35)
axis(1, at = seq(30, 90, by=20), labels=rep("", 4), padj = 0.25, lwd=-1, lwd.ticks=.85, tcl=-.35)
axis(2, at = c(1,2,5,10,20,30,50), labels=c(1,2,5,10,20,30,50), cex.axis = 1.85, las = 1, lwd=1.5, lwd.ticks=1.35)
mtext(side = 2, text = expression("Somatic cfDNA variants / Mb"), line = 5, cex = 1.85)
mtext(side = 1, text = "Age (years)", line = 4, cex = 1.5)
title(main="Biopsy matched", cex.main=1.5)

screen(zz[4])
plot(1, 1, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, ylim=c(1,50), xlim=c(20,90), log="y")
data = data.frame(burden_cancer[,c("patient_id", "age", "num_called", "bio_source", "subj_type"),drop=FALSE]) %>%
 	   filter(bio_source=="IMPACT-BAM_matched")
indx = unique(burden_cancer[!(burden_cancer$patient_id[burden_cancer$subj_type!="Control"] %in% data[,"patient_id"]),"patient_id",drop=TRUE])
if (length(indx)!=0) {
	tmp = matrix(NA, nrow=length(indx), ncol=5)
	tmp[,1] = indx
	tmp[,3] = 0
	tmp[,4] = "biopsy_matched"
	zzz = data.frame(clinical[clinical$patient_id %in% indx,c("patient_id","age","subj_type")])
	rownames(zzz) = zzz[,"patient_id"]
	tmp[,2] = zzz[indx,"age"]
	tmp[,5] = zzz[indx,"subj_type"]
	colnames(tmp) = colnames(data)
	tmp = data.frame(tmp)
	data = rbind(data, tmp)
}
p0 = zeroinfl(as.numeric(num_called) ~ as.numeric(age), dist = "poisson", data = data)
text(45, 49, cex = 1.25, labels = paste("(P = ", toupper(signif(summary(p0)$coefficients$count[2,4], 3)), ")", sep = ""), pos = 4)
y0 = fun_zerob(x=as.numeric(data[,"age"]), y=as.numeric(data[,"num_called"]), n=1000, seed=0)
x1 = seq(20, 90, l=100)
y1 = t(apply(y0, 1, quantile, probs=c(.025, .5, .975), na.rm=TRUE))
polygon(c(x1, rev(x1)), c(y1[,"97.5%"], rev(y1[,"2.5%"])), col=transparent_rgb("grey10", 55), border=NA)
points(x1, y1[,"50%"], type="l", col="grey50", lwd=3)
points(as.numeric(data$age), as.numeric(data$num_called), bg="#EAA411", col="#231F20", pch=21, cex=1.35, lwd=0.5)
axis(1, at = NULL, cex.axis = 1.5, padj = 0.25, lwd=1.5, lwd.ticks=1.35)
axis(1, at = seq(30, 90, by=20), labels=rep("", 4), padj = 0.25, lwd=-1, lwd.ticks=.85, tcl=-.35)
axis(2, at = c(1,2,5,10,20,30,50), labels=c(1,2,5,10,20,30,50), cex.axis = 1.85, las = 1, lwd=1.5, lwd.ticks=1.35)
mtext(side = 2, text = expression("Somatic cfDNA variants / Mb"), line = 5, cex = 1.85)
mtext(side = 1, text = "Age (years)", line = 4, cex = 1.5)
title(main="Biopsy subthreshold", cex.main=1.5)

close.screen(all.screens=TRUE)
dev.off()

#==================================================
# heat map of wbc-matched variants number of variants
#==================================================
variants_nohyper = variants %>%
 				   filter(!(patient_id %in% hypermutators$patient_id)) %>%
 				   filter(!(patient_id %in% msi_hypermutators$patient_id)) %>%
 				   filter(bio_source=="WBC_matched", is_nonsyn)
 				   
nx = variants_nohyper %>%
	dplyr::group_by(subj_type, gene) %>%
	dplyr::summarize(n = n())

ny = matrix(0, nrow = length(unique(nx %>% .[["gene"]])), ncol=4, dimnames=list(unique(nx %>% .[["gene"]]), c("Control", "Breast", "Lung", "Prostate")))
for (i in 1:ncol(ny)) {
	x = nx %>%
		filter(subj_type==colnames(ny)[i])
	ny[x$gene,i] = x$n
}
index = order(apply(ny, 1, sum), decreasing=TRUE)
ny = ny[index,,drop=FALSE]
pdf(file="../res/figure4/heatmap_wbc_matched_variants_by_var.pdf", width=5, height=21)
corrplot(corr=ny[1:30,,drop=FALSE], method="color", type="full", add=FALSE,
		 col=colorRampPalette(c(rep("#ffffff",2), "#19547b"))(100), bg=NULL, is.corr=FALSE,
		 addgrid.col="white", diag=TRUE, outline=FALSE, mar=c(1, 0, 1, 0), cl.pos="n",
 		 addCoef.col = "black", number.font = 1, number.cex = 1.15, tl.cex = 1.5, tl.col = "black", tl.offset = 0.5)
dev.off()

#==================================================
# heat map of wbc-matched variants gene recurrences
#==================================================
variants_nohyper = variants %>%
 				   filter(!(patient_id %in% hypermutators$patient_id)) %>%
 				   filter(!(patient_id %in% msi_hypermutators$patient_id)) %>%
 				   filter(bio_source=="WBC_matched", is_nonsyn)
subj_num_smry = variants_nohyper %>%
				group_by(subj_type) %>%
				summarise(subj_num = length(unique(patient_id))) %>%
				ungroup() %>%
				mutate(subj_type_num = paste0(subj_type, "(", subj_num, ")"))

gene_recurrences = variants_nohyper %>%
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
colnames(gene_stats_matrix) = unlist(lapply(strsplit(subj_num_smry$subj_type_num[idx_names], "(", fixed=TRUE), function(x) {x[1]}))

pdf(file="../res/figure4/heatmap_wbc_matched_variants_by_subj.pdf", width=5, height=21)
corrplot(corr=gene_stats_matrix[1:30,,drop=FALSE], method="color", type="full", add=FALSE,
		 col=colorRampPalette(c(rep("#ffffff",2), "#19547b"))(100), bg=NULL, is.corr=FALSE,
		 addgrid.col="white", diag=TRUE, outline=FALSE, mar=c(1, 0, 1, 0), cl.pos="n",
 		 addCoef.col = "black", number.font = 1, number.cex = 1.15, tl.cex = 1.5, tl.col = "black", tl.offset = 0.5)
dev.off()

#==================================================
# scatter of vaf in cfdna versus vaf in wbc
#==================================================
variants_nohyper = variants %>%
		   		   filter(!patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id))
cfdna_vs_gdna_vaf_plot = variants_nohyper %>%
						  filter(is_nonsyn) %>%
						  filter(!is.na(afmeancfdna)) %>%
						  filter(bio_source %in% c("WBC_matched","biopsy_matched","VUSo", "IMPACT-BAM_matched"))

pdf(file="../res/figure4/vaf2vaf_scatterplot.pdf", width=8, height=8)
par(mar = c(6.1, 6, 4.1, 1))
plot(1,1, type="n", xlim=c(0.009,15), ylim=c(0.009,15), xlab="", ylab="", axes=FALSE, frame.plot=FALSE, log="xy")
points(x=c(0.01, 15), y=c(0.01,15), type="l", lty=1, lwd=2, col="goldenrod3")
x = cfdna_vs_gdna_vaf_plot$afmeancfdna*100
y = cfdna_vs_gdna_vaf_plot$afmeangdna*100
x[x>15] = NA
y[y>15] = NA
for (i in c("VUSo", "WBC_matched", "biopsy_matched", "IMPACT-BAM_matched")) {
	index = cfdna_vs_gdna_vaf_plot$bio_source==i
	points(x[index], y[index], type="p", col="black", bg=variant_cols[i], pch=21, lwd=.5, cex=1.2)
}
axis(1, at = c(.01, .05, .10, .50, 1.0, 5.0, 10.0, 15.0), labels = c(".01", ".05", ".1", ".5", "1", "5", "10", ""), cex.axis = 1.5, padj = 0.25, lwd=1.5, lwd.ticks=1.35)
axis(1, at = 15, labels = "15", cex.axis = 1.5, padj = 0.25, lwd=-1, las=1)
axis(2, at = c(.01, .05, .10, .50, 1.0, 5.0, 10.0, 15.0), labels = c(".01", ".05", ".1", ".5", "1", "5", "10", ""), cex.axis = 1.5, las = 1, lwd=1.5, lwd.ticks=1.35)
axis(2, at = 15, labels = "15", cex.axis = 1.5, padj = 0.25, lwd=-1, las=1)
mtext(side = 1, text = "VAF in cfDNA (%)", line = 4, cex = 1.85)
mtext(side = 2, text = "VAF in WBC (%)", line = 4, cex = 1.85)
legend(x=0.01, y=19, pch=21, col="black", pt.bg=variant_cols[c("biopsy_matched", "IMPACT-BAM_matched", "VUSo", "WBC_matched")], pt.cex=1.55, pt.lwd=.5, legend=c("Biopsy matched", "Biopsy subthreshold", "VUSo", "WBC matched"), box.lwd=-1, cex=1.15)
dev.off()
