#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/figure3")) {
	dir.create("../res/figure3")
}

all_msi_df = read_tsv(file=url_msi_processed_data)

alias_df = read_tsv(file=url_subject_alias)

type_frame_df = data_frame(type=c("VB","VL","VP"), expand=c("Metastatic Breast","Metastatic Lung","Metastatic Prostate"), subject_type=c("Breast","Lung","Prostate"))

useful_msi_df = all_msi_df %>%
				left_join(alias_df %>% mutate(id=patient_id)) %>%
				mutate(id=alias) %>%
				dplyr::select(patient_id,original_msi,fixed_msi,subject_type,sample)
				
#==================================================
# scatter plot of msi scores
#==================================================
pdf(file="../res/figure3/msi_scatter.pdf", width=7, height=7)
par(mar = c(6.1, 6, 4.1, 1))
zz = split.screen(figs=matrix(c(0,1,0,1, 0,1,0,1, 0,1,0,1), nrow=3, ncol=4, byrow=TRUE))
x = useful_msi_df %>%
	filter(sample=="tumor")
x = (100*x$fixed_msi)
y = useful_msi_df %>%
	filter(sample=="cfdna")
y = (100*y$fixed_msi)
z = cohort_cols[as.vector(useful_msi_df$subject_type)]
z2 = rep(21, length(useful_msi_df$subject_type))

screen(zz[1])
plot(x[x<5], y[x<5], xlim=c(0,4.15), ylim=c(0,4), pch=z2[x<5], bg=unlist(lapply(z[x<5], transparent_rgb, 255)), col="#231F20", xlab="", ylab="", axes=FALSE, frame.plot=FALSE, cex=1.55, lwd=1)
axis(1, at = c(0,1,2), labels=c(0,1,2), cex.axis = 1.5, padj = 0.25)
axis(2, at = c(0,1,2), labels=c(0,1,2), cex.axis = 1.5, las = 1)
axis(1, at = c(2,2.9), labels=c("",""), lwd.ticks=-1)
axis(2, at = c(2,2.9), labels=c("",""), lwd.ticks=-1)
axis.break(axis=1, breakpos=2.13,pos=NA,bgcol="white",breakcol="black",style="slash",brw=0.02)
axis.break(axis=2, breakpos=2.3,pos=NA,bgcol="white",breakcol="black",style="slash",brw=0.02)
mtext(side = 1, text = "Tumor MSI score", line = 4, cex = 1.5)
mtext(side = 2, text = "cfDNA MSI score", line = 4, cex = 1.5)
legend(x=.03, y=4.1, legend=c("Breast", "Lung", "Prostate"), pt.bg=cohort_cols[c("Breast", "Lung", "Prostate")], col="#231F20", pch=21, pt.cex=1.55, box.lwd=-1, pt.lwd=1)
screen(zz[2])
plot(x[x>2 & x<10], y[x>2 & x<10], xlim=c(-10,25), ylim=c(0,20), pch=z2[x>5 & x<10], bg=unlist(lapply(z[x>5 & x<10], transparent_rgb, 255)), col="#231F20", xlab="", ylab="", axes=FALSE, frame.plot=FALSE, cex=1.55, lwd=1)
axis(1, at = c(9,15,20,25), labels=c(9,15,20,25), cex.axis = 1.5, padj = 0.25)
screen(zz[3])
plot(x[x>10], y[x>10], xlim=c(-10,25), ylim=c(-10,20), pch=z2[x>10], bg=unlist(lapply(z[x>10], transparent_rgb, 255)), col="#231F20", xlab="", ylab="", axes=FALSE, frame.plot=FALSE, cex=1.55, lwd=1)
axis(2, at = c(9,15,20), labels=c(9,15,20), cex.axis = 1.5, las = 1)
close.screen(all.screens=TRUE)
dev.off()

#==================================================
# scatter plot of tmb estimated using cfdna versus
# msk-impact
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
		   filter(is_nonsyn) %>%
		   filter(bio_source %in% c("VUSo","biopsy_matched", "biopsy_only"))
		   
summ_per_patient_cfdna = variants %>%
				   filter(subj_type!="Control") %>%
				   filter(bio_source %in% c("VUSo", "biopsy_matched")) %>%
 				   group_by(patient_id) %>%
 				   summarize(num_called = n()) %>%
 				   mutate(TMB = num_called/total_bed_Mb) %>%
 				   ungroup()
 
summ_per_patient_impact = variants %>%
				   filter(subj_type!="Control") %>%
				   filter(bio_source %in% c("biopsy_matched", "biopsy_only")) %>%
 				   group_by(patient_id) %>%
 				   summarize(num_called = n()) %>%
 				   mutate(TMB = num_called/total_bed_Mb) %>%
 				   ungroup()
 				   
pdf(file="../res/figure3/tmb_scatter.pdf", width=7, height=7)
par(mar = c(6.1, 6, 4.1, 1))
zz = split.screen(figs=matrix(c(0,1,0,1, 0,1,0,1, 0,1,0,1), nrow=3, ncol=4, byrow=TRUE))

valid_patient_ids = tracker_grail %>%
					filter(study=="TechVal") %>%
					filter(patient_id %in% tracker_impact$patient_id) %>%
					.[["patient_id"]]
x = y = rep(0, length(valid_patient_ids))
names(x) = names(y) = valid_patient_ids
x[summ_per_patient_impact[,1,drop=TRUE]] = summ_per_patient_impact[,3,drop=TRUE]
y[summ_per_patient_cfdna[,1,drop=TRUE]] = summ_per_patient_cfdna[,3,drop=TRUE]
subj_type = rep("Breast", length(valid_patient_ids))
subj_type[grepl("VL", valid_patient_ids)] = "Lung"
subj_type[grepl("VP", valid_patient_ids)] = "Prostate"
z1 = cohort_cols[subj_type]
z2 = rep(21, length(subj_type))
z2[y>30] = 24
screen(zz[1])
plot(x[y<150], y[y<150], xlim=c(0,40), ylim=c(0,150), pch=z2[y<150], bg=unlist(lapply(z1[y<150], transparent_rgb, 255)), col="#231F20", xlab="", ylab="", axes=FALSE, frame.plot=FALSE, cex=1.55, lwd=1)
points(x=c(-10,40), y=rep(22.6996,2), type="l", lty=3, col="goldenrod3", lwd=1.5)
points(x=rep(13.8,2), y=c(-10,150), type="l", lty=3, col="goldenrod3", lwd=1.5)
axis(1, at = c(0,10,20,30,40), labels=c(0,10,20,30,40), cex.axis = 1.5, padj = 0.25)
axis(2, at = c(0,25,50,75,100,110), labels=c(0,25,50,75,100,110), cex.axis = 1.5, las = 1)
axis(2, at = c(100,150), labels=c("",""), lwd.ticks=-1)
axis.break(axis=2, breakpos=120, pos=NA, bgcol="white",breakcol="black",style="slash",brw=0.02)
mtext(side = 1, text = "Tumor mutation burden", line = 4, cex = 1.5)
mtext(side = 2, text = "cfDNA mutation burden", line = 4, cex = 1.5)
legend(x=30, y=155, legend=c("Breast", "Lung", "Prostate"), col="#231F20", pt.bg=cohort_cols[c("Breast", "Lung", "Prostate")], pch=21, pt.cex=1.55, box.lwd=-1, pt.lwd=1)
legend(x=30, y=120, legend="Non-hyper-\nmutated", pch=21, col="#231F20", pt.bg="#231F20", pt.cex=1.55, box.lwd=-1, pt.lwd=1)
legend(x=30, y=105, legend="Hyper-\nmutated", pch=24, col="#231F20", pt.bg="#231F20", pt.cex=1.55, box.lwd=-1, pt.lwd=1)
screen(zz[2])
plot(x[y>=150], y[y>=150], xlim=c(0,50), ylim=c(150,600), pch=24, col="#231F20", bg=unlist(lapply(z1[y>=150], transparent_rgb, 255)), xlab="", ylab="", axes=FALSE, frame.plot=FALSE, cex=1.55, lwd=1)
axis(2, at = c(550,600), labels=c(550,600), cex.axis = 1.5, las = 1)
close.screen(all.screens=TRUE)
dev.off()

#==================================================
# bar plot of mutational signatures
#==================================================
variants_cfdna = variants %>%
 				 filter(patient_id %in% c(msi_hypermutators$patient_id, hypermutators$patient_id)) %>%
 				 filter(bio_source %in% c("VUSo", "biopsy_matched"))
base96 = mut.to.sigs.input(mut.ref = data.frame(variants_cfdna),
  						   sample.id = "patient_id",
  						   chr = "chrom",
  						   pos = "position",
  						   ref = "ref_orig",
  						   alt = "alt_orig")
sig2extract = c(1:6,10,13,15,20,26)
sigs = foreach (i=1:nrow(base96)) %dopar% {
  	patient_id = rownames(base96)[i]
  	tmp = whichSignatures(base96, sample.id=patient_id, signatures.ref=signatures.cosmic[sig2extract,,drop=FALSE], contexts.needed=TRUE)
  	x = as.vector(tmp$weights)
  	x = unlist(c(x, 1-sum(x)))
 	return(invisible(x))
}
sigs = do.call(cbind, sigs)
colnames(sigs) = rownames(base96)
rownames(sigs) = c(paste0("Signature ", sig2extract), "Other")
fsigs = matrix(0, nrow=7, ncol=ncol(sigs))
rownames(fsigs) = c("Aging", "APOBEC", "HRD", "MMR", "Smoking", "POLE", "Other")
colnames(fsigs) = colnames(sigs)
fsigs["Aging",] = sigs["Signature 1",] + sigs["Signature 5",]
fsigs["APOBEC",] = sigs["Signature 2",] + sigs["Signature 13",]
fsigs["HRD",] = sigs["Signature 3",]
fsigs["MMR",] = sigs["Signature 6",] + sigs["Signature 15",] + sigs["Signature 20",] + sigs["Signature 26",]
fsigs["Smoking",] = sigs["Signature 4",]
fsigs["POLE",] = sigs["Signature 10",]
fsigs["Other",] = sigs["Other",]
fsigs = fsigs*100
pdf(file="../res/figure3/mutational_signatures.pdf", width=15, height=7)
par(mar = c(6.1, 6, 4.1, 1))
zz = split.screen(figs=matrix(c(0,.7,0,.99, 0, .76, 0, 1), nrow=2, ncol=4, byrow=TRUE))
col = brewer.pal(n=nrow(fsigs), name="Paired")
screen(zz[1])
z = barplot(fsigs, col=col, names.arg=rep("", ncol(fsigs)), beside=FALSE, axes=FALSE, ylim=c(0,100), main="", cex.main=1.65, las=2)
mtext(side = 2, text = "% of signature", line = 4, cex = 1.55)
axis(2, at = NULL, labels=NULL, cex.axis = 1.55, las = 1, line=0)
screen(zz[2])
plot(0,0, type="n", xlim=c(0,10), ylim=c(0,10), xlab="", ylab="", axes=FALSE, frame.plot=FALSE)
legend("topright", pch=15, col=col, legend=rownames(fsigs), box.lwd=-1, pt.cex=1.75)
close.screen(all.screens=TRUE)
dev.off()

summ_per_patient_cfdna = variants_cfdna %>%
				   filter(bio_source %in% c("VUSo", "biopsy_matched")) %>%
 				   group_by(patient_id) %>%
 				   summarize(num_called = n()) %>%
 				   mutate(TMB = num_called/total_bed_Mb) %>%
 				   ungroup()
sig2extract = c(1:6,10,13,15,20,26)
rho = unlist(foreach (i=1:nrow(base96)) %dopar% {
  	patient_id = rownames(base96)[i]
  	tmp = whichSignatures(base96, sample.id=patient_id, signatures.ref=signatures.cosmic[sig2extract,,drop=FALSE], contexts.needed=TRUE)
  	x = cor(as.vector(tmp$tumor), as.vector(tmp$product))
 	return(invisible(x))
})
names(rho) = rownames(base96)
pdf(file="../res/figure3/mutational_signatures_pearson.pdf", width=15, height=3)
par(mar = c(6.1, 6, 4.1, 1))
zz = split.screen(figs=matrix(c(0,.7,0,.99), nrow=1, ncol=4, byrow=TRUE))
screen(zz[1])
plot(z, rho, type="o", pch=21, col="black", bg="white", cex=1.5, lwd=1.5, ylim=c(0.1,1.05), xlab="", ylab="", main="", frame.plot=FALSE, axes=FALSE)
axis(2, at = c(.2, .6, 1), labels=c(".2", ".6", "1"), cex.axis = 1.5, las = 1, line=2.7)
close.screen(all.screens=TRUE)
dev.off()

#==================================================
# venn diagram of number of mutations
#==================================================
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
		   
snv_vars = read_tsv(snv_file$scored, col_types = cols(.default = col_character())) %>%
		   type_convert()  %>%
		   mutate(level_2a = as.character(level_2a)) %>%
		   mutate(level_r1 = as.character(level_r1))
		   
indel_vars = read_tsv(indel_file$scored, col_types = cols(.default = col_character())) %>%
			 type_convert()  %>%
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

variants_nohyper = variants %>%
		   	filter(subj_type!="Control") %>%
		   	filter(!patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id)) %>%
		   	filter(is_nonsyn | bio_source=="biopsy_only") %>%
		   	filter(bio_source %in% c("VUSo","biopsy_matched", "biopsy_only", "IMPACT-BAM_matched"))

variants_hyper = variants %>%
		   	filter(subj_type!="Control") %>%
		   	filter(patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id)) %>%
		   	filter(is_nonsyn | bio_source=="biopsy_only") %>%
		   	filter(bio_source %in% c("VUSo","biopsy_matched", "biopsy_only", "IMPACT-BAM_matched"))

pander(sum(variants_nohyper$bio_source=="VUSo"))
pander(sum(variants_nohyper$subj_type=="Breast" & variants_nohyper$bio_source=="VUSo"))
pander(sum(variants_nohyper$subj_type=="Lung" & variants_nohyper$bio_source=="VUSo"))
pander(sum(variants_nohyper$subj_type=="Prostate" & variants_nohyper$bio_source=="VUSo"))

pander(sum(variants_hyper$bio_source=="VUSo"))
pander(sum(variants_hyper$subj_type=="Breast" & variants_hyper$bio_source=="VUSo"))
pander(sum(variants_hyper$subj_type=="Lung" & variants_hyper$bio_source=="VUSo"))
pander(sum(variants_hyper$subj_type=="Prostate" & variants_hyper$bio_source=="VUSo"))

pdf(file="../res/figure3/venn_diagram_somatic_vuso.pdf", width=14, height=9)
par(mar = c(6.1, 6, 4.1, 1))
fit = euler(c("VUSo_nonh" = sum(variants_nohyper$bio_source=="VUSo"),
			  "Biopsy_only_nonh" = sum(variants_nohyper$bio_source=="biopsy_only"),
			  "VUSo_nonh&Biopsy_only_nonh" = sum(variants_nohyper$bio_source=="biopsy_matched" | variants_nohyper$bio_source=="IMPACT-BAM_matched"),
			  "VUSo_hyp" = sum(variants_hyper$bio_source=="VUSo"),
			  "Biopsy_only_hyp" = sum(variants_hyper$bio_source=="biopsy_only"),
		      "VUSo_hyp&Biopsy_only_hyp" = sum(variants_hyper$bio_source=="biopsy_matched" | variants_hyper$bio_source=="IMPACT-BAM_matched")))
plot(fit, fills=rep(unlist(lapply(c(variant_cols["VUSo"], variant_cols["biopsy_only"]), transparent_rgb, 155)), 2))
dev.off()

#==================================================
# recist and psa
#==================================================
PSA = read_tsv(file=url_psa, col_types = cols(.default = col_character())) %>%
	  type_convert()
RECIST = read_tsv(file=url_recist, col_types = cols(.default = col_character())) %>%
	     type_convert()

pdf(file="../res/figure3/PSA.pdf", width=10, height=7)
par(mar = c(6.1, 9, 4.1, 1))
plot(PSA$Days, PSA$PSA, xlim=c(-150,750), ylim=c(0,14), type="l", col="grey10", xlab="", ylab="", axes=FALSE, frame.plot=FALSE, lwd=2.5)
index = 1:37
points(PSA$Days, PSA$PSA, type="p", pch=21, bg="grey10", col="grey10", lwd=1.5, cex=0.95)
axis(1, at = c(-150,0,150,300,450,600,750), labels=c(-150,0,150,300,450,600,750), cex.axis = 1.5, padj = 0.25, lwd=1.5)
axis(2, at = NULL, labels=NULL, cex.axis = 1.5, las = 1, lwd=1.5)
mtext(side = 1, text = "Time (days)", line = 4, cex = 1.5)
mtext(side = 2, text = "PSA (ng/ml)", line = 4, cex = 1.5)
dev.off()

pdf(file="../res/figure3/RECIST.pdf", width=10, height=7)
par(mar = c(6.1, 9, 4.1, 1))
plot(RECIST$Days, RECIST$RECIST1.1, xlim=c(0,750), ylim=c(-100,0), type="l", col="grey10", xlab="", ylab="", axes=FALSE, frame.plot=FALSE, lwd=2.5)
index = 1:13
points(RECIST$Days, RECIST$RECIST1.1, type="p", pch=21, bg="grey10", col="grey10", lwd=1.5, cex=0.95)
axis(1, at = c(0,150,300,450,600,750), labels=c(0,150,300,450,600,750), cex.axis = 1.5, padj = 0.25, lwd=1.5)
axis(2, at = NULL, labels=NULL, cex.axis = 1.5, las = 1, lwd=1.5)
mtext(side = 1, text = "Time (days)", line = 4, cex = 1.5)
mtext(side = 2, text = "Percent change in\ntumor size\n(RECIST v1.1)", line = 4, cex = 1.5)
dev.off()
