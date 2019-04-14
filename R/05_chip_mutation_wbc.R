#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/figure5")) {
	dir.create("../res/figure5")
}

#==================================================
# Extract mutation calls from wbc in all genes
#==================================================
if (FLAG) {
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
	 		   af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq)) %>%
	 		   mutate(ID_x = make_id_x(patient_id, chrom, position, ref, alt))
	 		   
	tracker_grail = read_csv(file=patient_tracker, col_types = cols(.default = col_character()))  %>%
 					type_convert()
	tracker_impact = read_csv(impact_tracker, col_types = cols(.default = col_character()))  %>%
					type_convert()
 	valid_patient_ids = tracker_grail %>%
 	  				    filter(patient_id %in% tracker_impact$patient_id | tissue == "Healthy") %>%
 	  				    filter(!(tissue %in% c("Breast", "Lung", "Prostate") & study=="Merlin")) %>%
 	  				    .[["patient_id"]]
 	clinical = read_tsv(clinical_file, col_types = cols(.default = col_character())) %>%
 			   type_convert() %>%
 			   rename(localisation = tissue)
 			   
 	valid_patient_ids = intersect(valid_patient_ids, clinical$patient_id)
 			   
 	som_vars_breast = read_tsv(file=somatic_vars_breast) %>%
 					  type_convert()
 	som_vars_lung = read_tsv(file=somatic_vars_lung) %>%
 					type_convert()
 	som_vars_prostate = read_tsv(file=somatic_vars_prostate) %>%
 						type_convert()
 	som_vars = bind_rows(som_vars_breast, som_vars_lung, som_vars_prostate) %>%
 			   mutate(ID_x = make_id_x(CASE, CHROM, POS, REF, ALT)) %>%
 			   mutate(ID_y = make_id_y(CASE, vcf2maf_HGVSp_Short))
 	load(cosmic_file)
 	cosmic_db = cosmic_db %>%
 				filter(SNP==0) %>%
 				filter(as.numeric(n)>=2) %>%
 				mutate(id = paste0("chr", CHROM, ":", POS, "_", REF, ">", ALT))
 	load(gnomad_file)
  	gnomad_db = gnomad_db %>%
  				filter(AF>=0.01) %>%
  				mutate(id = paste0("chr", CHROM, ":", POS, "_", REF, ">", ALT))
  	load(hotspot_file)
  	cancer_hotspot = cancer_hotspot %>%
  					 mutate(id = paste0("chr", Chromosome, ":", Start_Position, "_", Reference_Allele, ">", Tumor_Seq_Allele2))
	 		   
  	all_vars = read_tsv(wbc_variants$scored, col_types = cols(.default = col_character()))  %>%
 			   type_convert()
 	all_annotations = read_tsv(wbc_variants$annotations, comment="#", col_types = cols(.default = col_character()))  %>%
 					type_convert()
 	all_vars = bind_cols(all_vars, all_annotations) %>%
 			   left_join(clinical, by=c("patient_id","study")) %>%
 	           mutate(subj_type = ifelse(subj_type=="Healthy" | is.na(subj_type), "Control", subj_type)) %>%
 			   mutate(is_patient_valid = (patient_id %in% valid_patient_ids)) %>%
 			   mutate(ID_x = make_id_x(patient_id, chrom, pos, ref, alt)) %>%
 			   mutate(ID_y = make_id_y(patient_id, HGVSp_Short)) %>%
 			   mutate(loc_lng = str_c(chrom, ":", pos, "_", ref, ">", alt)) %>%
 			   mutate(loc_srt = str_c(chrom, ":", pos)) %>%
 			   mutate(c_panel = in_common_bed(chrom, pos)) %>%
 			   mutate(is_hypermutator = (patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id))) %>%
 			   mutate(is_lowdepth = (dpnobaq < 500)) %>%
 		   	   mutate(is_lowqual = (qualnobaq < 60)) %>%
 			   mutate(is_tumor_matched = (ID_x %in% som_vars$ID_x) | (ID_y %in% som_vars$ID_y)) %>%
 			   mutate(is_cfdna_matched = (ID_x %in% variants$ID_x[variants$bio_source %in% c("biopsy_matched", "biopsy_only", "germline", "IMPACT-BAM_matched")])) %>%
 			   mutate(in_cosmic = loc_lng %in% cosmic_db$id) %>%
  		   	   mutate(in_exac = !is.na(ExAC_AF) & ExAC_AF>=0.01) %>%
  		   	   mutate(in_gnomad = loc_lng %in% gnomad_db$id) %>%
  		   	   mutate(is_hotspot = loc_lng %in% cancer_hotspot$id)
  		   	   
  	save(all_vars, clinical, file=all_vars_and_clinical)
} else {
 	load(all_vars_and_clinical)
}

#==================================================
# Save variants
#==================================================
save_vars = all_vars
 
#==================================================
# Default filters
#==================================================
all_vars = all_vars %>%
 		   filter(is_patient_valid) %>%
 		   filter(c_panel) %>%
 		   filter(!is_hypermutator) %>%
 		   filter(!is_lowdepth) %>%
 		   filter(!is_lowqual) %>%
 		   filter(!is_tumor_matched) %>%
 		   filter(!is_cfdna_matched)
 
#==================================================
# < 5% recurrence | is_hotspot | frame-shifting
# in CH related gene
#==================================================
n_samples = all_vars %>%
 			distinct(patient_id) %>%
  		    count()
recurrence = all_vars %>%
   			 group_by(loc_lng) %>%
   			 count() %>%
   			 ungroup() %>%
   			 rename(n_recurrence=n) %>%
   			 mutate(f_recurrence=n_recurrence/n_samples$n)
all_vars = all_vars %>%
   		   left_join(recurrence)
all_vars = all_vars %>%
 		   filter(f_recurrence < 0.05 | is_hotspot | (Variant_Classification %in% c("Frame_Shift_Ins", "Frame_Shift_Del", "Nonsense_Mutation", "Nonstop_Mutation") & (SYMBOL %in% chip_genes)))
  		   
#==================================================
# Co-occurring indels filter
#==================================================
recurrence = all_vars %>%
   			 group_by(patient_id, loc_srt) %>%
   			 count() %>%
   			 ungroup() %>%
   			 rename(n_indel=n)
all_vars = all_vars %>%
		   left_join(recurrence)
all_vars = all_vars %>%
		   filter(!(n_indel > 1 & indel))
		   
#==================================================
# Variant classe filter
#==================================================
all_vars = all_vars %>%
 		   filter(Variant_Classification!="3'Flank") %>%
  		   filter(Variant_Classification!="3'UTR") %>%
  		   filter(Variant_Classification!="5'Flank") %>%
  		   filter(Variant_Classification!="5'UTR") %>%
  		   filter(Variant_Classification!="In_Frame_Del") %>%
  		   filter(Variant_Classification!="In_Frame_Ins") %>%
  		   filter(Variant_Classification!="Intron") %>%
  		   filter(Variant_Classification!="RNA") %>%
  		   filter(Variant_Classification!="Silent") %>%
  		   filter(Variant_Classification!="IGR") %>%
  		   filter(Variant_Classification!="Translation_Start_Site")
  		   
#==================================================
# HLA-A
#==================================================
all_vars = all_vars %>%
 		   filter(SYMBOL!="HLA-A")

#==================================================
# Germline filter
#==================================================
all_vars = all_vars %>%
 		   filter((adnobaq/dpnobaq)<=.3 | (Variant_Classification %in% c("Frame_Shift_Ins", "Frame_Shift_Del", "Nonsense_Mutation", "Nonstop_Mutation") & SYMBOL %in% chip_genes))
		   
#==================================================
# ExAC filter
#==================================================
all_vars = all_vars %>%
 		   filter(!in_exac)
 		   
#==================================================
# gnomAD filter
#==================================================
all_vars = all_vars %>%
 		   filter(!in_gnomad)
 		   
#==================================================
# Scatter plot mutation burden versus age for
# all individuals and all genes
#==================================================
burden_healthy = all_vars %>%
  		   		 filter(subj_type=="Control") %>%
  		   		 group_by(patient_id) %>%
  		   		 summarize(num_called = n()) %>%
  		   		 ungroup() %>%
  		   		 left_join(clinical)
patient_ids = save_vars %>%
			  filter(is_patient_valid) %>%
			  filter(subj_type=="Control") %>%
			  select(patient_id) %>%
			  distinct() %>%
			  .[["patient_id"]]
index = (clinical$patient_id %in% patient_ids) & !(clinical$patient_id %in% burden_healthy$patient_id)
tmp = clinical[index,,drop=FALSE]
tmp = cbind(tmp, num_called=rep(0, sum(index)))
burden_healthy = rbind(burden_healthy, tmp[,colnames(burden_healthy)])
burden_cancer = all_vars %>%
  		   		filter(subj_type!="Control") %>%
  		   		group_by(patient_id) %>%
  	 	   		summarize(num_called = n()) %>%
  	 	   		ungroup() %>%
  	 	   		left_join(clinical)
patient_ids = save_vars %>%
			  filter(is_patient_valid) %>%
			  filter(subj_type!="Control") %>%
			  select(patient_id) %>%
			  distinct() %>%
			  .[["patient_id"]]			  
index = (clinical$patient_id %in% patient_ids) & !(clinical$patient_id %in% burden_cancer$patient_id)
tmp = clinical[index,,drop=FALSE]
tmp = cbind(tmp, num_called=rep(0, sum(index)))
burden_cancer = rbind(burden_cancer, tmp[,colnames(burden_cancer)])
 
pdf(file="../res/figure5/mutational_burden_all_genes.pdf")
par(mar = c(6.1, 6.5, 4.1, 1.1))
plot(0, 0, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, ylim=c(.5,50), xlim=c(20,100), log="y")
data = data.frame(burden_cancer[,c("age", "num_called", "study"),drop=FALSE])
p1 = glm(num_called ~ age, family="poisson", data=data)
#text(20, 50, cex = .95, labels = bquote(AIC == .(summary(p1)$aic, 3)), pos = 4)
#text(20, 40, cex = .95, labels = paste("(P = ", toupper(signif(summary(p1)$coefficients[2,4], 3)), ")", sep = ""), pos = 4)
#text(20, 30, cex = .95, labels = paste("N = ", sum(burden_cancer$num_called),sep = ""), pos = 4)
p1 = ggplot_build(qplot(age, num_called, data=burden_cancer) + stat_smooth(method = "glm", method.args = list(family = "poisson"), level = .99))$data[[2]]
points(burden_cancer$age, ifelse(burden_cancer$num_called==0, 0.5, burden_cancer$num_called), bg=transparentRgb("#EAA411", 255), pch=21, cex=1.05, lwd=0)
data = data.frame(burden_healthy[,c("age", "num_called", "study"),drop=FALSE])
p2 = glm(num_called ~ age, family="poisson", data=data)
#text(20, 20, cex = .95, labels = bquote(AIC == .(summary(p2)$aic, 3)), pos = 4)
#text(20, 15, cex = .95, labels = paste("(P = ", toupper(signif(summary(p2)$coefficients[2,4], 3)), ")", sep = ""), pos = 4)
#text(20, 10, cex = .95, labels = paste("N = ", sum(burden_healthy$num_called),sep = ""), pos = 4)
p2 = ggplot_build(qplot(age, num_called, data=burden_healthy) + stat_smooth(method = "glm", method.args = list(family = "poisson"), level = .99))$data[[2]]
points(burden_healthy$age, ifelse(burden_healthy$num_called==0, 0.5, burden_healthy$num_called), bg=transparentRgb("#1F6784", 255), pch=21, cex=1.05, lwd=0)
points(p1[,"x"], p1[,"y"], type="l", col="#EAA411", lwd=2)
points(p2[,"x"], p2[,"y"], type="l", col="#1F6784", lwd=2)
axis(1, at = NULL, cex.axis = 1.5, padj = 0.25, lwd=1.5, lwd.ticks=1.35)
axis(1, at = seq(30, 90, by=20), labels=rep("", 4), padj = 0.25, lwd=-1, lwd.ticks=.85, tcl=-.35)
axis(2, at = c(.5, 1, 2, 5, 10, 20, 50), labels = c("0", "1", "2", "5", "10", "20", "50"), cex.axis = 1.5, las = 1, lwd=1.5, lwd.ticks=1.35)
mtext(side = 1, text = "Age (years)", line = 4, cex = 1.5)
mtext(side = 2, text = "Number of variants", line = 4.5, cex = 1.5)
legend(x=20, y=55, pch=21, legend=c("Cancer","Control"), box.lwd=-1, pt.cex=1.05, cex=.95, pt.bg=unlist(lapply(c("#EAA411", "#1F6784"), transparentRgb, 255)), pt.lwd=0)
dev.off()
 
#==================================================
# Scatter plot mutation burden versus age for
# all individuals and CH genes only
#==================================================
burden_healthy = all_vars %>%
				 filter(SYMBOL %in% chip_genes) %>%
  		   		 filter(subj_type=="Control") %>%
  		   		 group_by(patient_id) %>%
   		   		 summarize(num_called = n()) %>%
  		   		 ungroup() %>%
  		   		 left_join(clinical)
patient_ids = save_vars %>%
			  filter(is_patient_valid) %>%
			  filter(subj_type=="Control") %>%
			  select(patient_id) %>%
			  distinct() %>%
			  .[["patient_id"]]
index = (clinical$patient_id %in% patient_ids) & !(clinical$patient_id %in% burden_healthy$patient_id)
tmp = clinical[index,,drop=FALSE]
tmp = cbind(tmp, num_called=rep(0, sum(index)))
burden_healthy = rbind(burden_healthy, tmp[,colnames(burden_healthy)])
burden_cancer = all_vars %>%
				filter(SYMBOL %in% chip_genes) %>%
  		   		filter(subj_type!="Control") %>%
  		   		group_by(patient_id) %>%
  	 	   		summarize(num_called = n()) %>%
  	 	   		ungroup() %>%
  	 	   		left_join(clinical)
patient_ids = save_vars %>%
			  filter(is_patient_valid) %>%
			  filter(subj_type!="Control") %>%
			  select(patient_id) %>%
			  distinct() %>%
			  .[["patient_id"]]			  
index = (clinical$patient_id %in% patient_ids) & !(clinical$patient_id %in% burden_cancer$patient_id)
tmp = clinical[index,,drop=FALSE]
tmp = cbind(tmp, num_called=rep(0, sum(index)))
burden_cancer = rbind(burden_cancer, tmp[,colnames(burden_cancer)])
 
pdf(file="../res/figure5/mutational_burden_ch_genes.pdf")
par(mar = c(6.1, 6.5, 4.1, 1.1))
plot(0, 0, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, ylim=c(.5,50), xlim=c(20,100), log="y")
data = data.frame(burden_cancer[,c("age", "num_called", "study"),drop=FALSE])
p1 = glm(num_called ~ age, family="poisson", data=data)
#text(20, 50, cex = .95, labels = bquote(AIC == .(summary(p1)$aic, 3)), pos = 4)
#text(20, 40, cex = .95, labels = paste("(P = ", toupper(signif(summary(p1)$coefficients[2,4], 3)), ")", sep = ""), pos = 4)
#text(20, 30, cex = .95, labels = paste("N = ", sum(burden_cancer$num_called),sep = ""), pos = 4)
p1 = ggplot_build(qplot(age, num_called, data=burden_cancer) + stat_smooth(method = "glm", method.args = list(family = "poisson"), level = .99))$data[[2]]
points(burden_cancer$age, ifelse(burden_cancer$num_called==0, 0.5, burden_cancer$num_called), bg=transparentRgb("#EAA411", 255), pch=21, cex=1.05, lwd=0)
data = data.frame(burden_healthy[,c("age", "num_called", "study"),drop=FALSE])
p2 = glm(num_called ~ age, family="poisson", data=data)
#text(20, 20, cex = .95, labels = bquote(AIC == .(summary(p2)$aic, 3)), pos = 4)
#text(20, 15, cex = .95, labels = paste("(P = ", toupper(signif(summary(p2)$coefficients[2,4], 3)), ")", sep = ""), pos = 4)
#text(20, 10, cex = .95, labels = paste("N = ", sum(burden_healthy$num_called),sep = ""), pos = 4)
p2 = ggplot_build(qplot(age, num_called, data=burden_healthy) + stat_smooth(method = "glm", method.args = list(family = "poisson"), level = .99))$data[[2]]
points(burden_healthy$age, ifelse(burden_healthy$num_called==0, 0.5, burden_healthy$num_called), bg=transparentRgb("#1F6784", 255), pch=21, cex=1.05, lwd=0)
points(p1[,"x"], p1[,"y"], type="l", col="#EAA411", lwd=2)
points(p2[,"x"], p2[,"y"], type="l", col="#1F6784", lwd=2)
axis(1, at = NULL, cex.axis = 1.5, padj = 0.25, lwd=1.5, lwd.ticks=1.35)
axis(1, at = seq(30, 90, by=20), labels=rep("", 4), padj = 0.25, lwd=-1, lwd.ticks=.85, tcl=-.35)
axis(2, at = c(.5, 1, 2, 5, 10, 20, 50), labels = c("0", "1", "2", "5", "10", "20", "50"), cex.axis = 1.5, las = 1, lwd=1.5, lwd.ticks=1.35)
mtext(side = 1, text = "Age (years)", line = 4, cex = 1.5)
mtext(side = 2, text = "Number of variants", line = 4.5, cex = 1.5)
legend(x=20, y=55, pch=21, legend=c("Cancer","Control"), box.lwd=-1, pt.cex=1.05, cex=.95, pt.bg=unlist(lapply(c("#EAA411", "#1F6784"), transparentRgb, 255)), pt.lwd=0)
dev.off()
  
#==================================================
# Barplot of mutation burden for all individuals
#==================================================
snvs = read_tsv(somatic_snvs_grail$scored, col_types = cols(.default = col_character()))  %>%
	   mutate(ID=paste0(patient_id, ":", chrom, ":", position, ":", ref_orig, ">", alt_orig))
indels = read_tsv(somatic_indels_grail$scored, col_types = cols(.default = col_character()))  %>%
 		 mutate(ID=paste0(patient_id, ":", chrom, ":", position, ":", ref_orig, ">", alt_orig))
feature_names = intersect(colnames(indels), colnames(snvs))
som_vars_grail = bind_rows(indels[,feature_names], snvs[,feature_names])
all_vars = all_vars %>%
		   mutate(is_cfdna_matched = ID_x %in% som_vars_grail$ID)

burden_healthy = all_vars %>%
  		   		 filter(subj_type=="Control") %>%
  		   		 group_by(patient_id) %>%
  	 	   		 summarize(num_called = n()) %>%
  	 	   		 ungroup() %>%
  	 	   		 left_join(clinical)
patient_ids = save_vars %>%
			  filter(is_patient_valid) %>%
			  filter(subj_type=="Control") %>%
			  select(patient_id) %>%
			  distinct() %>%
			  .[["patient_id"]]
index = (clinical$patient_id %in% patient_ids) & !(clinical$patient_id %in% burden_healthy$patient_id)  	 	   		 
tmp = clinical[index,,drop=FALSE]
tmp = cbind(tmp, num_called=rep(0, sum(index)))
burden_healthy = rbind(burden_healthy, tmp[,colnames(burden_healthy)])
burden_healthy_incfdna = all_vars %>%
  		   		 filter(subj_type=="Control") %>%
  		   		 filter(is_cfdna_matched) %>%
  		   		 group_by(patient_id) %>%
  	 	   		 summarize(num_called = n()) %>%
  	 	   		 ungroup() %>%
  	 	   		 left_join(clinical)
index = (burden_healthy$patient_id %in% burden_healthy_incfdna$patient_id)
tmp = burden_healthy[!index,,drop=FALSE]
tmp$num_called = 0
burden_healthy_incfdna = rbind(burden_healthy_incfdna, tmp)
burden_healthy_ch = all_vars %>%
				    filter(SYMBOL %in% chip_genes) %>%
  		   		    filter(subj_type=="Control") %>%
  		   		    group_by(patient_id) %>%
  	 	   		    summarize(num_called = n()) %>%
  	 	   		    ungroup() %>%
  	 	   		    left_join(clinical)
index = (clinical$patient_id %in% patient_ids) & !(clinical$patient_id %in% burden_healthy_ch$patient_id)  
tmp = clinical[index,,drop=FALSE]
tmp = cbind(tmp, num_called=rep(0, sum(index)))
burden_healthy_ch = rbind(burden_healthy_ch, tmp[,colnames(burden_healthy_ch)])
burden_cancer = all_vars %>%
  		   		filter(subj_type!="Control") %>%
  		   		group_by(patient_id) %>%
  	 	   		summarize(num_called = n()) %>%
  	 	   		ungroup() %>%
  	 	   		left_join(clinical)
patient_ids = save_vars %>%
			  filter(is_patient_valid) %>%
			  filter(subj_type!="Control") %>%
			  select(patient_id) %>%
			  distinct() %>%
			  .[["patient_id"]]			  
index = (clinical$patient_id %in% patient_ids) & !(clinical$patient_id %in% burden_cancer$patient_id)
tmp = clinical[index,,drop=FALSE]
tmp = cbind(tmp, num_called=rep(0, sum(index)))
burden_cancer = rbind(burden_cancer, tmp[,colnames(burden_cancer)])
burden_cancer_incfdna = all_vars %>%
  		   		 filter(subj_type!="Control") %>%
  		   		 filter(is_cfdna_matched) %>%
  		   		 group_by(patient_id) %>%
  	 	   		 summarize(num_called = n()) %>%
  	 	   		 ungroup() %>%
  	 	   		 left_join(clinical)
index = (burden_cancer$patient_id %in% burden_cancer_incfdna$patient_id)
tmp = burden_cancer[!index,,drop=FALSE]
tmp$num_called = 0
burden_cancer_incfdna = rbind(burden_cancer_incfdna, tmp)
burden_cancer_ch = all_vars %>%
				   filter(SYMBOL %in% chip_genes) %>%
  		   		   filter(subj_type!="Control") %>%
  		   		   group_by(patient_id) %>%
  	 	   		   summarize(num_called = n()) %>%
  	 	   		   ungroup() %>%
  	 	   		   left_join(clinical)
index = (clinical$patient_id %in% patient_ids) & !(clinical$patient_id %in% burden_cancer_ch$patient_id)
tmp = clinical[index,,drop=FALSE]
tmp = cbind(tmp, num_called=rep(0, sum(index)))
burden_cancer_ch = rbind(burden_cancer_ch, tmp[,colnames(burden_cancer_ch)])

pdf(file="../res/figure5/mutational_burden_barplot_combined.pdf", width=21, height=9)
par(mar = c(6.1, 6, 4.1, 1))
zz = split.screen(figs=matrix(c(0,1,0.36,1, 0,1,0,0.595), nrow=2, ncol=4, byrow=TRUE))
screen(zz[1])
data_.cancer_all = data.frame(burden_cancer[,c("age", "num_called", "study"),drop=FALSE], row.names=burden_cancer$patient_id)
data_.healthy_all = data.frame(burden_healthy[,c("age", "num_called", "study"),drop=FALSE], row.names=burden_healthy$patient_id)
data_.cancer_cfdna = data.frame(burden_cancer_incfdna[,c("age", "num_called", "study"),drop=FALSE], row.names=burden_cancer_incfdna$patient_id)
data_.healthy_cfdna = data.frame(burden_healthy_incfdna[,c("age", "num_called", "study"),drop=FALSE], row.names=burden_healthy_incfdna$patient_id)
data_.cancer_ch = data.frame(burden_cancer_ch[,c("age", "num_called", "study"),drop=FALSE], row.names=burden_cancer_ch$patient_id)
data_.healthy_ch = data.frame(burden_healthy_ch[,c("age", "num_called", "study"),drop=FALSE], row.names=burden_healthy_ch$patient_id)
data_.cancer_ch = cbind(data_.cancer_ch, matrix(NA, nrow=nrow(data_.cancer_ch), ncol=length(chip_genes), dimnames=list(rownames(data_.cancer_ch), chip_genes)))
data_.healthy_ch = cbind(data_.healthy_ch, matrix(NA, nrow=nrow(data_.healthy_ch), ncol=length(chip_genes), dimnames=list(rownames(data_.healthy_ch), chip_genes)))
for (i in 1:nrow(data_.cancer_ch)) {
	for (j in 1:length(chip_genes)) {
 		index = all_vars$patient_id==rownames(data_.cancer_ch)[i] & all_vars$SYMBOL==chip_genes[j]
 		if (sum(index)!=0) {
 			data_.cancer_ch[i,chip_genes[j]] = max(100*all_vars[index,"adnobaq"]/all_vars[index,"dpnobaq"])
 		}
 	}
}
for (i in 1:nrow(data_.healthy_ch)) {
	for (j in 1:length(chip_genes)) {
 		index = all_vars$patient_id==rownames(data_.healthy_ch)[i] & all_vars$SYMBOL==chip_genes[j]
 		if (sum(index)!=0) {
 			data_.healthy_ch[i,chip_genes[j]] = max(100*all_vars[index,"adnobaq"]/all_vars[index,"dpnobaq"])
 		}
 	}
}

data_.cancer_cfdna = data_.cancer_cfdna[order(data_.cancer_cfdna$num_called),,drop=FALSE]
data_.healthy_cfdna = data_.healthy_cfdna[order(data_.healthy_cfdna$num_called),,drop=FALSE]
data_.cancer_all = data_.cancer_all[rownames(data_.cancer_cfdna),,drop=FALSE]
data_.healthy_all = data_.healthy_all[rownames(data_.healthy_cfdna),,drop=FALSE]
data_.cancer_ch = data_.cancer_ch[rownames(data_.cancer_cfdna),,drop=FALSE]
data_.healthy_ch = data_.healthy_ch[rownames(data_.healthy_cfdna),,drop=FALSE]

data_.cancer_all = data_.cancer_all[order(data_.cancer_all$num_called),,drop=FALSE]
data_.healthy_all = data_.healthy_all[order(data_.healthy_all$num_called),,drop=FALSE]
data_.cancer_cfdna = data_.cancer_cfdna[rownames(data_.cancer_all),,drop=FALSE]
data_.healthy_cfdna = data_.healthy_cfdna[rownames(data_.healthy_all),,drop=FALSE]
data_.cancer_ch = data_.cancer_ch[rownames(data_.cancer_all),,drop=FALSE]
data_.healthy_ch = data_.healthy_ch[rownames(data_.healthy_all),,drop=FALSE]

age_cat = list(start = c(20, 40, 50, 60, 70, 80),
			   end	 = c(40, 50, 60, 70, 80, 100))

index_.cancer = list()
for (i in 1:6) {
	indx = which(data_.cancer_all$age>=age_cat$start[i] & data_.cancer_all$age<age_cat$end[i])
	index_.cancer[[i]] = indx
}
index_.healthy = list()
for (i in 1:6) {
	indx = which(data_.healthy_all$age>=age_cat$start[i] & data_.healthy_all$age<age_cat$end[i])
	index_.healthy[[i]] = indx
}
data_all = data_cfdna = data_ch = list()
for (i in 1:6) {
	data_all[[i]] = rbind(data_.healthy_all[index_.healthy[[i]],], data_.cancer_all[index_.cancer[[i]],])
	data_cfdna[[i]] = rbind(data_.healthy_cfdna[index_.healthy[[i]],], data_.cancer_cfdna[index_.cancer[[i]],])
	data_ch[[i]] = rbind(data_.healthy_ch[index_.healthy[[i]],], data_.cancer_ch[index_.cancer[[i]],])
}
data_all = do.call(rbind, data_all)
data_cfdna = do.call(rbind, data_cfdna)
data_cfdna = data_cfdna[rownames(data_all),,drop=FALSE]
data_ch = do.call(rbind, data_ch)
data_ch = data_ch[rownames(data_all),,drop=FALSE]
cols = ifelse(grepl("MSK", rownames(data_all)), "#EAA411", "#1F6784")
zzz = barplot(data_all[,"num_called"], col=unlist(lapply(cols, transparentRgb, 105)), border=cols, space=.15, axes=FALSE, ylim=c(-3,41), xlim=c(0,220))
barplot(data_cfdna[,"num_called"], col=unlist(lapply(cols, transparentRgb, 255)), border=NA, space=.15, lwd=0.01, add=TRUE, axes=FALSE)
axis(2, at = NULL, cex.axis = 1.25, las = 1)
mtext(side = 2, text = "Number of variants", line = 3.75, cex = 1.35)
legend(x=0, y=42, legend=c("Cancer: WBC only", "Cancer: WBC and cfDNA", "Control: WBC only", "Control: WBC and cfDNA"), col=c(transparentRgb("#EAA411"), "#EAA411", transparentRgb("#1F6784"), "#1F6784"), pch=15, box.lwd=-1, cex=1.05, pt.cex=1.65)
legend(x=65, y=42, legend=c("DNMT3A", "TP53", "TET2", "ASXL1", "PPM1D", "Other CH"),
				      col=c("#01985C", "#F7DC02", "#1E4665", "#FF7175", "#8CDB5E", "#3D98D3"),
				      pch=19, box.lwd=-1, cex=1.05, pt.cex=1.45, text.font=3)
ix = NULL
for (i in 1:6) {
 	indx = which(data_all$age>=age_cat$start[i] & data_all$age<age_cat$end[i])
 	ix = c(ix, length(indx))
}
ix = cumsum(ix)
colss = unlist(lapply(c("#F7F7F7", "#D9D9D9", "#BDBDBD", "#969696", "#636363", "#252525"), transparentRgb, 225))
colss2 = c(rep("black", 4), rep("white", 2))
for (i in 2:length(ix)) {
 	rect(xleft=zzz[ix[i-1],1]+0.7, xright=zzz[ix[i],1]+.18, ybottom=-3, ytop=-.5, col=colss[i], border="white", lwd=1)
 	if (i!=length(ix)) {
		text(x=mean(c(zzz[ix[i-1],1]+0.7, zzz[ix[i],1])+.18), y=-1.61, labels=paste0(age_cat$start[i], " - ", age_cat$end[i]-1), cex=.75, col=colss2[i])
	} else {
		text(x=mean(c(zzz[ix[i-1],1]+0.7, zzz[ix[i],1])+.18), y=-1.61, labels=paste0(age_cat$start[i], "+"), cex=.75, col=colss2[i])
	}
}
rect(xleft=.5, xright=zzz[ix[1],1]+.18, ybottom=-3, ytop=-.5, col=colss[1], border="white", lwd=1)
text(x=mean(c(.5, zzz[ix[1],1]+.18)), y=-1.61, labels=paste0(age_cat$start[1], " - ", age_cat$end[1]-1), cex=.65, col=colss2[1])
screen(zz[2])
max_vaf = vector(mode="numeric", length=nrow(data_all))
for (i in 1:nrow(data_all)) {
 	if (sum(all_vars$patient_id==rownames(data_all)[i])!=0) {
 		max_vaf[i] = max(100*(all_vars$adnobaq/all_vars$dpnobaq)[all_vars$patient_id==rownames(data_all)[i]])
	} else {
		max_vaf[i] = 0.05
	}
}
plot(1,1, type="n", xlab="", ylab="", axes=FALSE, frame.plot=FALSE, xlim=c(0,220), ylim=c(100,0.05), log="y")
 for (i in 1:nrow(zzz)) {
  	points(rep(zzz[i,1], 2), c(0.05, (max_vaf)[i]), type="l", col=cols[i], lwd=1.35)
}
to_plot = c("DNMT3A", "TP53", "TET2", "ASXL1", "PPM1D")
to_plot_also = c("JAK2", "RUNX1", "SF3B1", "SRSF2", "IDH1", "IDH2", "U2AF1", "CBL", "ATM", "CHEK2")
cols = c("#01985C", "#F7DC02", "#1E4665", "#FF7175", "#8CDB5E", "#3D98D3")
for (i in 1:length(to_plot)) {
 	points(zzz[,1], data_ch[,to_plot[i]], type="p", pch=19, col=cols[i], cex=.85)
}
for (i in 1:length(to_plot_also)) {
	points(zzz[,1], data_ch[,to_plot_also[i]], type="p", pch=19, col=cols[6], cex=.75)
}
axis(2, at = c(0.05,0.1,0.5,1,5,10,50), labels = c("0.05","0.1","0.5","1.0","5.0","10.0","50.0"), cex.axis = 1.25, las = 1)
mtext(side = 2, text = "VAF (%)", line = 3.75, cex = 1.35)
close.screen(all.screens=TRUE)
dev.off()

#==================================================
# Stacked barplot of dominant CH mutation
#==================================================
pdf(file="../res/figure5/stacked_barplot_max_vaf_ch_only_tx_combined.pdf", width=7, height=5)
par(mar = c(6.1, 6, 4.1, 1))
zz = split.screen(figs=matrix(c(0,.9,.4,1, 0,.9,0,.6, 0,1,0,1), nrow=3, ncol=4, byrow=TRUE))
screen(zz[1])
index_control = grepl("W", rownames(data_ch), fixed=TRUE)
index_late = grepl("V", rownames(data_ch), fixed=TRUE)
tmp = data_ch[index_control,-c(1:3),drop=FALSE]
tmp  = cbind(tmp[,to_plot,drop=FALSE],
			 "OTHER CH"=apply(tmp[,to_plot_also], 1, max, na.rm=TRUE))
tmp[is.infinite(as.matrix(tmp))] = NA
inx = apply(tmp, 1, function(x) sum(is.na(x)))==ncol(tmp)
f1 = vector(mode="numeric", length=length(inx))
names(f1) = rownames(tmp)
f1[inx] = "NO CH"
f1[!inx] = colnames(tmp)[apply(tmp[!inx,,drop=FALSE], 1, which.max)]
f1 = table(f1)
tmp = data_ch[index_late,-c(1:3),drop=FALSE]
tmp  = cbind(tmp[,to_plot,drop=FALSE],
 			 "OTHER CH"=apply(tmp[,to_plot_also], 1, max, na.rm=TRUE))
tmp[is.infinite(as.matrix(tmp))] = NA
inx = apply(tmp, 1, function(x) sum(is.na(x)))==ncol(tmp)
f3 = vector(mode="numeric", length=length(inx))
names(f3) = rownames(tmp)
f3[inx] = "NO CH"
f3[!inx] = colnames(tmp)[apply(tmp[!inx,,drop=FALSE], 1, which.max)]
f3 = table(f3)
arg_names = c(to_plot, "OTHER CH", "NO CH")
f = n = matrix(0, nrow=length(arg_names), ncol=2, dimnames=list(arg_names, c("Cancer","Control")))
f[names(f1),"Control"] = 100*f1/sum(f1)
f[names(f3),"Cancer"] = 100*f3/sum(f3)
n[names(f1),"Control"] = f1
n[names(f3),"Cancer"] = f3
pc1 = chisq.test(x=f, simulate.p.value=TRUE, B=10000)
barplot(f, horiz=TRUE, col=c(cols,"grey80"), axes=FALSE, xlab="", las=1)
screen(zz[2])
tx = read.csv(file="../res/etc/prior_tx_techval_0818.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE, row.names=1)
index = !(rownames(data_ch) %in% rownames(tx))
tx2 = matrix(0,nrow=sum(index), ncol=4)
rownames(tx2) = rownames(data_ch)[index]
colnames(tx2) = colnames(tx)[1:4]
tx = rbind(tx[,1:4,drop=FALSE], tx2)
tx = tx[rownames(data_ch),,drop=FALSE]
index_control = !(as.logical(tx$prior_rt) | as.logical(tx$prior_chemo))
index_trc = !index_control
tmp = data_ch[index_control,-c(1:3),drop=FALSE]
tmp  = cbind(tmp[,to_plot,drop=FALSE],
			 "OTHER CH"=apply(tmp[,to_plot_also], 1, max, na.rm=TRUE))
tmp[is.infinite(as.matrix(tmp))] = NA
inx = apply(tmp, 1, function(x) sum(is.na(x)))==ncol(tmp)
f1 = vector(mode="numeric", length=length(inx))
names(f1) = rownames(tmp)
f1[inx] = "NO CH"
f1[!inx] = colnames(tmp)[apply(tmp[!inx,,drop=FALSE], 1, which.max)]
f1 = table(f1)
tmp = data_ch[index_trc,-c(1:3),drop=FALSE]
tmp  = cbind(tmp[,to_plot,drop=FALSE],
			 "OTHER CH"=apply(tmp[,to_plot_also], 1, max, na.rm=TRUE))
tmp[is.infinite(as.matrix(tmp))] = NA
inx = apply(tmp, 1, function(x) sum(is.na(x)))==ncol(tmp)
f2 = vector(mode="numeric", length=length(inx))
names(f2) = rownames(tmp)
f2[inx] = "NO CH"
f2[!inx] = colnames(tmp)[apply(tmp[!inx,,drop=FALSE], 1, which.max)]
f2 = table(f2)
arg_names = c(to_plot, "OTHER CH", "NO CH")
f = n = matrix(0, nrow=length(arg_names), ncol=2, dimnames=list(arg_names, c("RT or CT", "No prior\nRT or CT")))
f[names(f1),"No prior\nRT or CT"] = 100*f1/sum(f1)
f[names(f2),"RT or CT"] = 100*f2/sum(f2)
n[names(f1),"No prior\nRT or CT"] = f1
n[names(f2),"RT or CT"] = f2
pc2 = chisq.test(x=f, simulate.p.value=TRUE, B=10000)
barplot(f, horiz=TRUE, col=c(cols,"grey80"), axes=FALSE, xlab="", las=1)
axis(1, at = NULL, cex.axis = 1.35, padj = 0.25, lwd=1.25, lwd.ticks=1.25, line=.5)
mtext(side = 1, text = "Fraction of patients (%)", line = 4, cex = 1.45)
screen(zz[3])
plot(0, 0, type="n", xlab="", ylab="", xlim=c(0,100), ylim=c(0,100), axes=FALSE, frame.plot=FALSE)
points(c(92,93), c(95,95), type="l")
points(c(92,93), c(77,77), type="l")
points(c(93,93), c(77,95), type="l")
text(x=94, y=86, labels="*", cex=.5)
points(c(92,93), c(23,23), type="l")
points(c(92,93), c(5,5), type="l")
points(c(93,93), c(5,23), type="l")
text(x=94, y=14, labels="**", cex=.5)
close.screen(all.screens=TRUE)
dev.off()

#==================================================
# Bubble plot of ratio indels to snvs
#==================================================
top_ch = c("DNMT3A", "TP53", "TET2", "ASXL1", "PPM1D")
other_ch = c("JAK2", "RUNX1", "SF3B1", "SRSF2", "IDH1", "IDH2", "U2AF1", "CBL", "ATM", "CHEK2")

index = grepl("W", all_vars$patient_id)
tmp_vars = all_vars[index,,drop=FALSE]
index = tmp_vars[,"SYMBOL",drop=TRUE] %in% c(top_ch, other_ch)
tmp_vars[!index,"SYMBOL"] = "OTHER"
index = tmp_vars[,"SYMBOL",drop=TRUE] %in% other_ch
tmp_vars[index,"SYMBOL"] = "OTHER CH"
tmp_vars = tmp_vars %>%
		   mutate(is_truncating = nchar(ref)!=1 | nchar(alt)!=1 | Variant_Classification=="Nonsense_Mutation" | Variant_Classification=="Nonstop_Mutation")
variants_control = table(tmp_vars$is_truncating, tmp_vars$SYMBOL)/sum(grepl("W", unique(all_vars$patient_id)))

index = grepl("V", all_vars$patient_id)
tmp_vars = all_vars[index,,drop=FALSE]
index = tmp_vars[,"SYMBOL",drop=TRUE] %in% c(top_ch, other_ch)
tmp_vars[!index,"SYMBOL"] = "OTHER"
index = tmp_vars[,"SYMBOL",drop=TRUE] %in% other_ch
tmp_vars[index,"SYMBOL"] = "OTHER CH"
tmp_vars = tmp_vars %>%
		   mutate(is_truncating = nchar(ref)!=1 | nchar(alt)!=1 | Variant_Classification=="Nonsense_Mutation" | Variant_Classification=="Nonstop_Mutation")
variants_cancer = table(tmp_vars$is_truncating, tmp_vars$SYMBOL)/sum(grepl("V", unique(all_vars$patient_id)))

tx = read.csv(file="../res/etc/prior_tx_techval_0818.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
index = !(all_vars$patient_id %in% tx[tx$prior_rt==1 | tx$prior_chemo==1,1])
tmp_vars = all_vars[index,,drop=FALSE]
index = tmp_vars[,"SYMBOL",drop=TRUE] %in% c(top_ch, other_ch)
tmp_vars[!index,"SYMBOL"] = "OTHER"
index = tmp_vars[,"SYMBOL",drop=TRUE] %in% other_ch
tmp_vars[index,"SYMBOL"] = "OTHER CH"
tmp_vars = tmp_vars %>%
		   mutate(is_truncating = nchar(ref)!=1 | nchar(alt)!=1 | Variant_Classification=="Nonsense_Mutation" | Variant_Classification=="Nonstop_Mutation")
variants_tn = table(tmp_vars$is_truncating, tmp_vars$SYMBOL)/length(unique(all_vars$patient_id[!(all_vars$patient_id %in% tx[tx$prior_rt==1 | tx$prior_chemo==1,1])]))

index_rct = all_vars$patient_id %in% tx[tx$prior_rt==1 | tx$prior_chemo==1,1]
tmp_vars = all_vars[index_rct,,drop=FALSE]
index = tmp_vars[,"SYMBOL",drop=TRUE] %in% c(top_ch, other_ch)
tmp_vars[!index,"SYMBOL"] = "OTHER"
index = tmp_vars[,"SYMBOL",drop=TRUE] %in% other_ch
tmp_vars[index,"SYMBOL"] = "OTHER CH"
tmp_vars = tmp_vars %>%
		   mutate(is_truncating = nchar(ref)!=1 | nchar(alt)!=1 | Variant_Classification=="Nonsense_Mutation" | Variant_Classification=="Nonstop_Mutation")
variants_rctx = table(tmp_vars$is_truncating, tmp_vars$SYMBOL)/length(unique(all_vars$patient_id[all_vars$patient_id %in% tx[tx$prior_rt==1 | tx$prior_chemo==1,1]]))

pdf(file="../res/figure5/corrplot_snvs_to_indels_ch_only.pdf", width=10, height=10)
m = rbind(apply(variants_control, 2, sum),
		  apply(variants_cancer, 2, sum),
		  apply(variants_tn, 2, sum),
		  apply(variants_rctx, 2, sum))
m2 = rbind(100*variants_control[2,]/apply(variants_control, 2, sum),
		   100*variants_cancer[2,]/apply(variants_cancer, 2, sum),
		   100*variants_tn[2,]/apply(variants_tn, 2, sum),
		   100*variants_rctx[2,]/apply(variants_rctx, 2, sum))
m = m[,c("DNMT3A", "TP53", "TET2", "ASXL1", "PPM1D", "OTHER CH"),drop=FALSE]
m2 = m2[,c("DNMT3A", "TP53", "TET2", "ASXL1", "PPM1D", "OTHER CH"),drop=FALSE]
colnames(m) = colnames(m2) = c("DNMT3A", "TP53", "TET2", "ASXL1", "PPM1D", "Other CH")
rownames(m) = rownames(m2) = c("Control", "Cancer", "Treatment naïve", "RT or CT")
corrplot2(corr=m, corr2=m2, method="circle", is.corr=FALSE, addgrid.col=NA,
		  col=colorRampPalette(c("#ffffff", rep("#c33764", 4), "#1d2671"))(100),
		  cl.lim=c(-0.01,1.6), tl.cex=1.25, tl.col="black", cl.pos = "n")
dev.off()