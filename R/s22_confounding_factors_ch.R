#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/rebuttal")) {
	dir.create("../res/rebuttal")
}

#==================================================
# Age by cancer status and treatment history
#==================================================
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

index_control = grepl("W", valid_patient_ids, fixed=TRUE)
index_cancer = grepl("V", valid_patient_ids, fixed=TRUE)

tmp = clinical %>%
	  dplyr::select(patient_id, age) %>%
	  mutate(cat = "") %>%
	  mutate(cat = ifelse(patient_id %in% valid_patient_ids[index_control], "Control", cat)) %>%
	  mutate(cat = ifelse(patient_id %in% valid_patient_ids[index_cancer], "Cancer", cat)) %>%
	  filter(cat != "") %>%
	  mutate(facet = "Age by cancer status")
	  
p = signif(wilcox.test(tmp$age[tmp$cat=="Cancer"], tmp$age[tmp$cat=="Control"], alternative="two.sided", correct=FALSE)$p.value, 3)
		
plot.0 = ggplot(tmp, aes(x=cat, y=age, fill = cat)) +
		 geom_boxplot(stat = "boxplot") +
		 scale_fill_manual(values = rev(c("#FDAE61", "salmon"))) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="", y="Age (years)\n") +
 		 coord_cartesian(ylim = c(20,100)) +
		 guides(fill=FALSE) +
		 facet_wrap(~facet) +
		 theme_bw(base_size=15) +
  		 annotate(geom="text", x=1.5, y=100, label=paste0("P = ", p))
	  	
pdf(file="../res/rebuttal/Age_by_Cancer.pdf", width=5, height=5)
print(plot.0)
dev.off()

tx = read.csv(file=url_prior_tx, header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
	 rename(patient_id = ID)
clinical = 	left_join(clinical, tx, by="patient_id") %>%
			mutate(prior_rt = ifelse(subj_type=="Healthy", 0, prior_rt)) %>%
			mutate(prior_chemo = ifelse(subj_type=="Healthy", 0, prior_chemo)) %>%
			filter(patient_id %in% valid_patient_ids)

tmp = clinical %>%
	  dplyr::select(patient_id, age) %>%
	  mutate(cat = "No RT/CT") %>%
	  mutate(cat = ifelse(clinical$prior_rt==1 | clinical$prior_chemo==1, "RT/CT", cat)) %>%
	  mutate(cat = factor(cat, levels=c("RT/CT","No RT/CT"))) %>%
	  mutate(facet = "Age by treatment history")
	  
p = signif(wilcox.test(tmp$age[tmp$cat=="No RT/CT"], tmp$age[tmp$cat=="RT/CT"], alternative="two.sided", correct=FALSE)$p.value, 3)
		
plot.0 = ggplot(tmp, aes(x=cat, y=age, fill = cat)) +
		 geom_boxplot(stat = "boxplot") +
		 scale_fill_manual(values = rev(c("#FDAE61", "salmon"))) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="", y="Age (years)\n") +
 		 coord_cartesian(ylim = c(20,100)) +
		 guides(fill=FALSE) +
		 facet_wrap(~facet) +
		 theme_bw(base_size=15) +
  		 annotate(geom="text", x=1.5, y=100, label=paste0("P = ", p))
	  	
pdf(file="../res/rebuttal/Age_by_Treatment.pdf", width=5, height=5)
print(plot.0)
dev.off()

#==================================================
# Smoking history by treatment history
#==================================================
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
clinical = clinical %>%
		   filter(patient_id %in% valid_patient_ids)
smoking_history_prostate = read_csv(file=url_smoking_history$prostate, col_types = cols(.default = col_character()))  %>%
						   type_convert() %>%
						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
smoking_history_breast = read_csv(file=url_smoking_history$breast, col_types = cols(.default = col_character()))  %>%
						 type_convert() %>%
						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
		   
tx = read.csv(file=url_prior_tx, header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
	 rename(patient_id = ID)
clinical = 	left_join(clinical, tx, by="patient_id") %>%
			mutate(prior_rt = ifelse(subj_type=="Healthy", 0, prior_rt)) %>%
			mutate(prior_chemo = ifelse(subj_type=="Healthy", 0, prior_chemo))

tmp = clinical %>%
	  dplyr::select(patient_id, smoking_history=have_you_ever_smoked_) %>%
	  mutate(cat = "No RT/CT") %>%
	  mutate(cat = ifelse(clinical$prior_rt==1 | clinical$prior_chemo==1, "RT/CT", cat)) %>%
	  mutate(cat = factor(cat, levels=c("RT/CT","No RT/CT"))) %>%
	  filter(!is.na(smoking_history)) %>%
	  group_by(cat, smoking_history) %>%
	  count_() %>%
	  mutate(facet = "Smoking by treatment history")
tmp[1:2,"n"] = 100*tmp[1:2,"n"]/sum(tmp[1:2,"n"])
tmp[3:4,"n"] = 100*tmp[3:4,"n"]/sum(tmp[3:4,"n"])
tmp = as.data.frame(tmp) %>%
	  mutate(cat = paste0(cat, "\n", smoking_history))
	  
plot.0 = ggplot(tmp, aes(x=cat, y=n, fill = smoking_history)) +
		 geom_bar(stat="identity", color="black") +
		 scale_fill_manual(values = rev(c("#FDAE61", "salmon"))) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="", y=" Frequency (%)\n") +
		 guides(fill=FALSE) +
		 coord_cartesian(ylim = c(0,100)) +
		 facet_wrap(~facet) +
		 theme_bw(base_size=15)
	  	
pdf(file="../res/rebuttal/Smoking_by_Treatment.pdf", width=5, height=5)
print(plot.0)
dev.off()

tmp = clinical %>%
	  dplyr::select(patient_id, smoking_history=have_you_ever_smoked_) %>%
	  mutate(cat = "No RT/CT") %>%
	  mutate(cat = ifelse(clinical$prior_rt==1 | clinical$prior_chemo==1, "RT/CT", cat)) %>%
	  mutate(cat = factor(cat, levels=c("RT/CT","No RT/CT"))) %>%
	  filter(!is.na(smoking_history))
p = fisher.test(table(tmp$smoking_history, tmp$cat))

#==================================================
# Smoking history as confounding factor
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

#==================================================
# Scatter plot of mutation burden versus age for
# biopsy matched, VUSo and WBC matched mutations
# adjusted for smoking history
#==================================================
tracker_grail = read_csv(file=patient_tracker, col_types = cols(.default = col_character()))  %>%
				type_convert()
tracker_impact = read_csv(impact_tracker, col_types = cols(.default = col_character()))  %>%
				 type_convert()
valid_patient_ids = tracker_grail %>%
  				    filter(patient_id %in% tracker_impact$patient_id | tissue == "Healthy") %>%
  				    filter(!(tissue %in% c("Breast", "Lung", "Prostate") & study=="Merlin")) %>%
  				    .[["patient_id"]]
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
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
data = data %>% type_convert()
p0 = zeroinfl(num_called ~ age, dist = "poisson", data = data)

smoking_history_prostate = read_csv(file=url_smoking_history$prostate, col_types = cols(.default = col_character()))  %>%
						   type_convert() %>%
						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
smoking_history_breast = read_csv(file=url_smoking_history$breast, col_types = cols(.default = col_character()))  %>%
						 type_convert() %>%
						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
data = left_join(data, clinical %>% dplyr::select(patient_id, smoking_history = have_you_ever_smoked_), by="patient_id") %>%
	   type_convert()
p1 = zeroinfl(num_called ~ age + smoking_history, dist = "poisson", data = data)
tmp = data %>%
	  mutate(subj_type = ifelse(subj_type=="Control", "Control", "Cancer")) %>%
	  mutate(subj_type = as.factor(subj_type)) %>%
	  mutate(bio_source = "WBC-matched") %>%
	  mutate(num_called = ifelse(num_called==0, 0.5, num_called))
	   
plot.0 = ggplot(tmp , aes(x = age, y = num_called, fill = subj_type)) +
		 geom_point(alpha=1, size=3.5, shape = 21, color = "#231F20") +
		 scale_fill_manual(values = c("Control"="#FDAE61", "Cancer"="salmon")) +
		 geom_smooth(formula = y ~ x, method="glm", aes(x = age, y = num_called), data=tmp %>% filter(num_called!=0.5), inherit.aes = FALSE, color="grey50", fill="grey50", fullrange=TRUE) +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.85), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\n Age (years)\n", y="Somatic cfDNA variants / Mb\n") +
 		 scale_y_log10(
 		 	breaks = function(x) { c(.5, 1, 2, 5, 10, 20, 30, 50) },
 		 	labels = function(x) { c("0", "1", "2", "5", "10", "20", "30", "50") }
 		 ) +
		 coord_cartesian(xlim = c(20, 90), ylim = c(.5, 50)) +
		 annotation_logticks(side="l") +
		 guides(fill=guide_legend(title=c("Cancer status"))) +
		 annotate(geom="text", x=55, y=50, label=paste0("p = ", signif(summary(p0)$coefficients$count[2,4], 3)), size=4.5) +
		 annotate(geom="text", x=56.4, y=37, label=paste0("p* = ", signif(summary(p1)$coefficients$count[2,4], 3)), size=4.5)
	  	
pdf(file="../res/rebuttal/WBC_matched_vs_Age_Combined.pdf", width=5.5, height=6)
print(plot.0)
dev.off()

burden = bind_rows(burden_cancer, burden_healthy)
data = data.frame(burden[,c("patient_id", "age", "num_called", "subj_type", "bio_source"),drop=FALSE]) %>%
 	   filter(bio_source=="VUSo")
indx = unique(burden[!(burden$patient_id %in% data[,"patient_id"]),"patient_id",drop=TRUE])
if (length(indx)!=0) {
	tmp = matrix(NA, nrow=length(indx), ncol=5)
	tmp[,1] = indx
	tmp[,3] = 0
	tmp[,5] = "VUSo"
	zzz = data.frame(clinical[clinical$patient_id %in% indx,c("patient_id","age","subj_type")])
	rownames(zzz) = zzz[,"patient_id"]
	tmp[,2] = zzz[indx,"age"]
	tmp[,4] = zzz[indx,"subj_type"]
	colnames(tmp) = colnames(data)
	tmp = data.frame(tmp)
	data = rbind(data, tmp)
}
data = data %>% type_convert()
p0 = zeroinfl(num_called ~ age, dist = "poisson", data = data)

smoking_history_prostate = read_csv(file=url_smoking_history$prostate, col_types = cols(.default = col_character()))  %>%
						   type_convert() %>%
						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
smoking_history_breast = read_csv(file=url_smoking_history$breast, col_types = cols(.default = col_character()))  %>%
						 type_convert() %>%
						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
data = left_join(data, clinical %>% dplyr::select(patient_id, smoking_history = have_you_ever_smoked_), by="patient_id") %>%
	   type_convert()
p1 = zeroinfl(num_called ~ age + smoking_history, dist = "poisson", data = data)

tmp = data %>%
	  mutate(subj_type = ifelse(subj_type=="Control", "Control", "Cancer")) %>%
	  mutate(subj_type = as.factor(subj_type)) %>%
	  mutate(num_called = ifelse(num_called==0, 0.5, num_called))
	   
plot.0 = ggplot(tmp , aes(x = age, y = num_called, fill = subj_type)) +
		 geom_point(alpha=1, size=3.5, shape = 21, color = "#231F20") +
		 scale_fill_manual(values = c("Control"="#FDAE61", "Cancer"="salmon")) +
		 geom_smooth(formula = y ~ x, method="glm", aes(x = age, y = num_called), data=tmp %>% filter(num_called!=0.5), inherit.aes = FALSE, color="grey50", fill="grey50", fullrange=TRUE) +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.85), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\n Age (years)\n", y="Somatic cfDNA variants / Mb\n") +
 		 scale_y_log10(
 		 	breaks = function(x) { c(.5, 1, 2, 5, 10, 20, 30, 50) },
 		 	labels = function(x) { c("0", "1", "2", "5", "10", "20", "30", "50") }
 		 ) +
		 coord_cartesian(xlim = c(20, 90), ylim = c(.5, 50)) +
		 annotation_logticks(side="l") +
		 guides(fill=guide_legend(title=c("Cancer status"))) +
		 annotate(geom="text", x=55, y=50, label=paste0("p = ", signif(summary(p0)$coefficients$count[2,4], 3)), size=4.5) +
		 annotate(geom="text", x=54, y=37, label=paste0("p* = ", signif(summary(p1)$coefficients$count[2,4], 3)), size=4.5)
	  	
pdf(file="../res/rebuttal/VUSo_vs_Age_Combined.pdf", width=5.5, height=6)
print(plot.0)
dev.off()

burden = bind_rows(burden_cancer, burden_healthy)
data = data.frame(burden[,c("patient_id", "age", "num_called", "subj_type", "bio_source"),drop=FALSE]) %>%
 	   filter(bio_source=="biopsy_matched")
indx = unique(burden[!(burden$patient_id %in% data[,"patient_id"]),"patient_id",drop=TRUE])
if (length(indx)!=0) {
	tmp = matrix(NA, nrow=length(indx), ncol=5)
	tmp[,1] = indx
	tmp[,3] = 0
	tmp[,5] = "biopsy_matched"
	zzz = data.frame(clinical[clinical$patient_id %in% indx,c("patient_id","age","subj_type")])
	rownames(zzz) = zzz[,"patient_id"]
	tmp[,2] = zzz[indx,"age"]
	tmp[,4] = zzz[indx,"subj_type"]
	colnames(tmp) = colnames(data)
	tmp = data.frame(tmp)
	data = rbind(data, tmp)
}
data = data %>% type_convert()
p0 = zeroinfl(num_called ~ age, dist = "poisson", data = data)
signif(summary(p0)$coefficients$count[2,4], 3)

smoking_history_prostate = read_csv(file=url_smoking_history$prostate, col_types = cols(.default = col_character()))  %>%
						   type_convert() %>%
						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
smoking_history_breast = read_csv(file=url_smoking_history$breast, col_types = cols(.default = col_character()))  %>%
						 type_convert() %>%
						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
data = left_join(data, clinical %>% dplyr::select(patient_id, smoking_history = have_you_ever_smoked_), by="patient_id") %>%
	   type_convert()
p1 = zeroinfl(num_called ~ age + smoking_history, dist = "poisson", data = data)
signif(summary(p0)$coefficients$count[2,4], 3)
tmp = data %>%
	  mutate(subj_type = ifelse(subj_type=="Control", "Control", "Cancer")) %>%
	  mutate(subj_type = as.factor(subj_type)) %>%
	  mutate(bio_source = "Biopsy-matched") %>%
	  mutate(num_called = ifelse(num_called==0, 0.5, num_called)) %>%
	  filter(subj_type=="Cancer")
	   
plot.0 = ggplot(tmp , aes(x = age, y = num_called, fill = subj_type)) +
		 geom_point(alpha=1, size=3.5, shape = 21, color = "#231F20") +
		 scale_fill_manual(values = c("Control"="#FDAE61", "Cancer"="salmon")) +
		 geom_smooth(formula = y ~ x, method="glm", aes(x = age, y = num_called), data=tmp %>% filter(num_called!=0.5), inherit.aes = FALSE, color="grey50", fill="grey50", fullrange=TRUE) +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.85), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\n Age (years)\n", y="Somatic cfDNA variants / Mb\n") +
 		 scale_y_log10(
 		 	breaks = function(x) { c(.5, 1, 2, 5, 10, 20, 30, 50) },
 		 	labels = function(x) { c("0", "1", "2", "5", "10", "20", "30", "50") }
 		 ) +
		 coord_cartesian(xlim = c(20, 90), ylim = c(.5, 50)) +
		 annotation_logticks(side="l") +
		 theme(legend.position="none") +
		 annotate(geom="text", x=55, y=50, label=paste0("p = ", signif(summary(p0)$coefficients$count[2,4], 3)), size=4.5) +
		 annotate(geom="text", x=56.4, y=37, label=paste0("p* = ", signif(summary(p1)$coefficients$count[2,4], 3)), size=4.5)
	  	
pdf(file="../res/rebuttal/Biopsy_matched_vs_Age_Combined.pdf", width=5.5, height=6)
print(plot.0)
dev.off()

burden = bind_rows(burden_cancer, burden_healthy)
data = data.frame(burden[,c("patient_id", "age", "num_called", "subj_type", "bio_source"),drop=FALSE]) %>%
 	   filter(bio_source=="IMPACT-BAM_matched")
indx = unique(burden[!(burden$patient_id %in% data[,"patient_id"]),"patient_id",drop=TRUE])
if (length(indx)!=0) {
	tmp = matrix(NA, nrow=length(indx), ncol=5)
	tmp[,1] = indx
	tmp[,3] = 0
	tmp[,5] = "IMPACT-BAM_matched"
	zzz = data.frame(clinical[clinical$patient_id %in% indx,c("patient_id","age","subj_type")])
	rownames(zzz) = zzz[,"patient_id"]
	tmp[,2] = zzz[indx,"age"]
	tmp[,4] = zzz[indx,"subj_type"]
	colnames(tmp) = colnames(data)
	tmp = data.frame(tmp)
	data = rbind(data, tmp)
}
data = data %>% type_convert()
p0 = zeroinfl(num_called ~ age, dist = "poisson", data = data)
signif(summary(p0)$coefficients$count[2,4], 3)

smoking_history_prostate = read_csv(file=url_smoking_history$prostate, col_types = cols(.default = col_character()))  %>%
						   type_convert() %>%
						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
smoking_history_breast = read_csv(file=url_smoking_history$breast, col_types = cols(.default = col_character()))  %>%
						 type_convert() %>%
						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
data = left_join(data, clinical %>% dplyr::select(patient_id, smoking_history = have_you_ever_smoked_), by="patient_id") %>%
	   type_convert()
p1 = zeroinfl(num_called ~ age + smoking_history, dist = "poisson", data = data)
signif(summary(p0)$coefficients$count[2,4], 3)
tmp = data %>%
	  mutate(subj_type = ifelse(subj_type=="Control", "Control", "Cancer")) %>%
	  mutate(subj_type = as.factor(subj_type)) %>%
	  mutate(bio_source = "Biopsy-subthreshold") %>%
	  mutate(num_called = ifelse(num_called==0, 0.5, num_called)) %>%
	  filter(subj_type=="Cancer")
	   
plot.0 = ggplot(tmp , aes(x = age, y = num_called, fill = subj_type)) +
		 geom_point(alpha=1, size=3.5, shape = 21, color = "#231F20") +
		 scale_fill_manual(values = c("Control"="#FDAE61", "Cancer"="salmon")) +
		 geom_smooth(formula = y ~ x, method="glm", aes(x = age, y = num_called), data=tmp %>% filter(num_called!=0.5), inherit.aes = FALSE, color="grey50", fill="grey50", fullrange=TRUE) +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.85), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\n Age (years)\n", y="Somatic cfDNA variants / Mb\n") +
 		 scale_y_log10(
 		 	breaks = function(x) { c(.5, 1, 2, 5, 10, 20, 30, 50) },
 		 	labels = function(x) { c("0", "1", "2", "5", "10", "20", "30", "50") }
 		 ) +
		 coord_cartesian(xlim = c(20, 90), ylim = c(.5, 50)) +
		 annotation_logticks(side="l") +
		 theme(legend.position="none") +
		 annotate(geom="text", x=55, y=50, label=paste0("p = ", signif(summary(p0)$coefficients$count[2,4], 3)), size=4.5) +
		 annotate(geom="text", x=55.5, y=37, label=paste0("p* = ", signif(summary(p1)$coefficients$count[2,4], 3)), size=4.5)
	  	
pdf(file="../res/rebuttal/Biopsy_subthreshold_vs_Age_Combined.pdf", width=5.5, height=6)
print(plot.0)
dev.off()
