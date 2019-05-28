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
	  select(patient_id, age) %>%
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

tx = read.csv(file="../res/etc/prior_tx_techval_0818.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
	 rename(patient_id = ID)
clinical = 	left_join(clinical, tx, by="patient_id") %>%
			mutate(prior_rt = ifelse(subj_type=="Healthy", 0, prior_rt)) %>%
			mutate(prior_chemo = ifelse(subj_type=="Healthy", 0, prior_chemo)) %>%
			filter(patient_id %in% valid_patient_ids)

tmp = clinical %>%
	  select(patient_id, age) %>%
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
smoking_history_prostate = read_csv(file="../res/etc/MSK_VP_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						   type_convert() %>%
						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
smoking_history_breast = read_csv(file="../res/etc/MSK_VB_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						 type_convert() %>%
						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
		   
tx = read.csv(file="../res/etc/prior_tx_techval_0818.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
	 rename(patient_id = ID)
clinical = 	left_join(clinical, tx, by="patient_id") %>%
			mutate(prior_rt = ifelse(subj_type=="Healthy", 0, prior_rt)) %>%
			mutate(prior_chemo = ifelse(subj_type=="Healthy", 0, prior_chemo))

tmp = clinical %>%
	  select(patient_id, smoking_history=have_you_ever_smoked_) %>%
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
	  select(patient_id, smoking_history=have_you_ever_smoked_) %>%
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
data = data %>%
	   type_convert()
p0 = zeroinfl(num_called ~ age, dist = "poisson", data = data)


smoking_history_prostate = read_csv(file="../res/etc/MSK_VP_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						   type_convert() %>%
						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
smoking_history_breast = read_csv(file="../res/etc/MSK_VB_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						 type_convert() %>%
						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
data = left_join(data, clinical %>% select(patient_id, smoking_history = have_you_ever_smoked_), by="patient_id") %>%
	   type_convert()
p1 = zeroinfl(num_called ~ age + smoking_history, dist = "poisson", data = data)
tmp = data %>%
	  mutate(subj_type = ifelse(subj_type=="Control", "Control", "Cancer")) %>%
	  mutate(subj_type = as.factor(subj_type)) %>%
	  filter(num_called !=0) %>%
	  mutate(bio_source = "WBC-matched")
	   
plot.0 = ggplot(tmp , aes(x = age, y = num_called, fill = subj_type)) +
		 geom_point(alpha=.85, size=3.5, shape = 21, color = "#231F20") +
		 scale_fill_manual(values = c("Control"="#FDAE61", "Cancer"="salmon")) +
		 geom_smooth(formula = y ~ x, method="glm", aes(x = age, y = num_called), data=tmp, inherit.aes = FALSE, color="grey50", fill="grey50") +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.85), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\n Age (years)\n", y="Somatic cfDNA variants / Mb\n") +
 		 scale_y_log10(
 		 	breaks = function(x) { c(1, 2, 5, 10, 20, 30, 50) },
 		 	labels = function(x) { c("1", "2", "5", "10", "20", "30", "50") }
 		 ) +
		 coord_cartesian(xlim = c(20, 90), ylim = c(1, 50)) +
		 annotation_logticks(side="l") +
		 guides(fill=guide_legend(title=c("Cancer status"))) +
		 annotate(geom="text", x=55, y=50, label=paste0("p = ", signif(summary(p0)$coefficients$count[2,4], 3)), size=3.5) +
		 annotate(geom="text", x=56.4, y=40, label=paste0("p* = ", signif(summary(p1)$coefficients$count[2,4], 3)), size=3.5)
	  	
pdf(file="../res/rebuttal/WBC_matched_vs_Age_Combined.pdf", width=5, height=6)
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
	tmp[,5] = "WBC_matched"
	zzz = data.frame(clinical[clinical$patient_id %in% indx,c("patient_id","age","subj_type")])
	rownames(zzz) = zzz[,"patient_id"]
	tmp[,2] = zzz[indx,"age"]
	tmp[,4] = zzz[indx,"subj_type"]
	colnames(tmp) = colnames(data)
	tmp = data.frame(tmp)
	data = rbind(data, tmp)
}
data = data %>%
	   type_convert()
p0 = zeroinfl(num_called ~ age, dist = "poisson", data = data)

smoking_history_prostate = read_csv(file="../res/etc/MSK_VP_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						   type_convert() %>%
						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
smoking_history_breast = read_csv(file="../res/etc/MSK_VB_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						 type_convert() %>%
						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
data = left_join(data, clinical %>% select(patient_id, smoking_history = have_you_ever_smoked_), by="patient_id") %>%
	   type_convert()
p1 = zeroinfl(num_called ~ age + smoking_history, dist = "poisson", data = data)

tmp = data %>%
	  mutate(subj_type = ifelse(subj_type=="Control", "Control", "Cancer")) %>%
	  mutate(subj_type = as.factor(subj_type)) %>%
	  filter(num_called !=0)
	   
plot.0 = ggplot(tmp , aes(x = age, y = num_called, fill = subj_type)) +
		 geom_point(alpha=.85, size=3.5, shape = 21, color = "#231F20") +
		 scale_fill_manual(values = c("Control"="#FDAE61", "Cancer"="salmon")) +
		 geom_smooth(formula = y ~ x, method="glm", aes(x = age, y = num_called), data=tmp, inherit.aes = FALSE, color="grey50", fill="grey50") +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.85), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\n Age (years)\n", y="Somatic cfDNA variants / Mb\n") +
 		 scale_y_log10(
 		 	breaks = function(x) { c(1, 2, 5, 10, 20, 30, 50) },
 		 	labels = function(x) { c("1", "2", "5", "10", "20", "30", "50") }
 		 ) +
		 coord_cartesian(xlim = c(20, 90), ylim = c(1, 50)) +
		 annotation_logticks(side="l") +
		 guides(fill=guide_legend(title=c("Cancer status"))) +
		 annotate(geom="text", x=55, y=50, label=paste0("p = ", signif(summary(p0)$coefficients$count[2,4], 3)), size=3.5) +
		 annotate(geom="text", x=54, y=40, label=paste0("p* = ", signif(summary(p1)$coefficients$count[2,4], 3)), size=3.5)
	  	
pdf(file="../res/rebuttal/VUSo_vs_Age_Combined.pdf", width=5, height=6)
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
	tmp[,5] = "WBC_matched"
	zzz = data.frame(clinical[clinical$patient_id %in% indx,c("patient_id","age","subj_type")])
	rownames(zzz) = zzz[,"patient_id"]
	tmp[,2] = zzz[indx,"age"]
	tmp[,4] = zzz[indx,"subj_type"]
	colnames(tmp) = colnames(data)
	tmp = data.frame(tmp)
	data = rbind(data, tmp)
}
data = data %>%
	   type_convert()
p0 = zeroinfl(num_called ~ age, dist = "poisson", data = data)
signif(summary(p0)$coefficients$count[2,4], 3)

smoking_history_prostate = read_csv(file="../res/etc/MSK_VP_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						   type_convert() %>%
						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
smoking_history_breast = read_csv(file="../res/etc/MSK_VB_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						 type_convert() %>%
						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
data = left_join(data, clinical %>% select(patient_id, smoking_history = have_you_ever_smoked_), by="patient_id") %>%
	   type_convert()
p1 = zeroinfl(num_called ~ age + smoking_history, dist = "poisson", data = data)
signif(summary(p0)$coefficients$count[2,4], 3)
tmp = data %>%
	  mutate(subj_type = ifelse(subj_type=="Control", "Control", "Cancer")) %>%
	  mutate(subj_type = as.factor(subj_type)) %>%
	  filter(num_called !=0) %>%
	  mutate(bio_source = "Biopsy-matched")
	   
plot.0 = ggplot(tmp , aes(x = age, y = num_called, fill = subj_type)) +
		 geom_point(alpha=.85, size=3.5, shape = 21, color = "#231F20") +
		 scale_fill_manual(values = c("Control"="#FDAE61", "Cancer"="salmon")) +
		 geom_smooth(formula = y ~ x, method="glm", aes(x = age, y = num_called), data=tmp, inherit.aes = FALSE, color="grey50", fill="grey50") +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.85), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\n Age (years)\n", y="Somatic cfDNA variants / Mb\n") +
 		 scale_y_log10(
 		 	breaks = function(x) { c(1, 2, 5, 10, 20, 30, 50) },
 		 	labels = function(x) { c("1", "2", "5", "10", "20", "30", "50") }
 		 ) +
		 coord_cartesian(xlim = c(20, 90), ylim = c(1, 50)) +
		 annotation_logticks(side="l") +
		 theme(legend.position="none") +
		 annotate(geom="text", x=55, y=50, label=paste0("p = ", signif(summary(p0)$coefficients$count[2,4], 3)), size=3.5) +
		 annotate(geom="text", x=56.4, y=40, label=paste0("p* = ", signif(summary(p1)$coefficients$count[2,4], 3)), size=3.5)
	  	
pdf(file="../res/rebuttal/Biopsy_matched_vs_Age_Combined.pdf", width=5, height=6)
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
	tmp[,5] = "WBC_matched"
	zzz = data.frame(clinical[clinical$patient_id %in% indx,c("patient_id","age","subj_type")])
	rownames(zzz) = zzz[,"patient_id"]
	tmp[,2] = zzz[indx,"age"]
	tmp[,4] = zzz[indx,"subj_type"]
	colnames(tmp) = colnames(data)
	tmp = data.frame(tmp)
	data = rbind(data, tmp)
}
data = data %>%
	   type_convert()
p0 = zeroinfl(num_called ~ age, dist = "poisson", data = data)
signif(summary(p0)$coefficients$count[2,4], 3)

smoking_history_prostate = read_csv(file="../res/etc/MSK_VP_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						   type_convert() %>%
						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
smoking_history_breast = read_csv(file="../res/etc/MSK_VB_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						 type_convert() %>%
						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
data = left_join(data, clinical %>% select(patient_id, smoking_history = have_you_ever_smoked_), by="patient_id") %>%
	   type_convert()
p1 = zeroinfl(num_called ~ age + smoking_history, dist = "poisson", data = data)
signif(summary(p0)$coefficients$count[2,4], 3)
tmp = data %>%
	  mutate(subj_type = ifelse(subj_type=="Control", "Control", "Cancer")) %>%
	  mutate(subj_type = as.factor(subj_type)) %>%
	  filter(num_called !=0) %>%
	  mutate(bio_source = "Biopsy-subthreshold")
	   
plot.0 = ggplot(tmp , aes(x = age, y = num_called, fill = subj_type)) +
		 geom_point(alpha=.85, size=3.5, shape = 21, color = "#231F20") +
		 scale_fill_manual(values = c("Control"="#FDAE61", "Cancer"="salmon")) +
		 geom_smooth(formula = y ~ x, method="glm", aes(x = age, y = num_called), data=tmp, inherit.aes = FALSE, color="grey50", fill="grey50") +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.85), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\n Age (years)\n", y="Somatic cfDNA variants / Mb\n") +
 		 scale_y_log10(
 		 	breaks = function(x) { c(1, 2, 5, 10, 20, 30, 50) },
 		 	labels = function(x) { c("1", "2", "5", "10", "20", "30", "50") }
 		 ) +
		 coord_cartesian(xlim = c(20, 90), ylim = c(1, 50)) +
		 annotation_logticks(side="l") +
		 theme(legend.position="none") +
		 annotate(geom="text", x=55, y=50, label=paste0("p = ", signif(summary(p0)$coefficients$count[2,4], 3)), size=3.5) +
		 annotate(geom="text", x=55.5, y=40, label=paste0("p* = ", signif(summary(p1)$coefficients$count[2,4], 3)), size=3.5)
	  	
pdf(file="../res/rebuttal/Biopsy_subthreshold_vs_Age_Combined.pdf", width=5, height=6)
print(plot.0)
dev.off()

#==================================================
# Scatter plot of mutation burden versus age for
# biopsy matched, VUSo and WBC matched mutations
# for smokers
#==================================================
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

smoking_history_prostate = read_csv(file="../res/etc/MSK_VP_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						   type_convert() %>%
						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
smoking_history_breast = read_csv(file="../res/etc/MSK_VB_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						 type_convert() %>%
						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
data = left_join(data, clinical %>% select(patient_id, smoking_history = have_you_ever_smoked_), by="patient_id") %>%
	   type_convert() %>%
	   filter(smoking_history=="Yes")
p0 = zeroinfl(num_called ~ age, dist = "poisson", data = data)
tmp = data %>%
	  mutate(subj_type = ifelse(subj_type=="Control", "Control", "Cancer")) %>%
	  mutate(subj_type = as.factor(subj_type)) %>%
	  filter(num_called !=0) %>%
	  mutate(bio_source = "WBC-matched")
	   
plot.0 = ggplot(tmp , aes(x = age, y = num_called, fill = subj_type)) +
		 geom_point(alpha=.85, size=3.5, shape = 21, color = "#231F20") +
		 scale_fill_manual(values = c("Control"="#FDAE61", "Cancer"="salmon")) +
		 geom_smooth(formula = y ~ x, method="glm", aes(x = age, y = num_called), data=tmp, inherit.aes = FALSE, color="grey50", fill="grey50") +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.85), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\n Age (years)\n", y="Somatic cfDNA variants / Mb\n") +
 		 scale_y_log10(
 		 	breaks = function(x) { c(1, 2, 5, 10, 20, 30, 50) },
 		 	labels = function(x) { c("1", "2", "5", "10", "20", "30", "50") }
 		 ) +
		 coord_cartesian(xlim = c(20, 90), ylim = c(1, 50)) +
		 annotation_logticks(side="l") +
		 theme(legend.position="none") +
		 annotate(geom="text", x=55, y=50, label=paste0("p = ", signif(summary(p0)$coefficients$count[2,4], 3)), size=3.5)
	  	
pdf(file="../res/rebuttal/WBC_matched_vs_Age_Smoker.pdf", width=5, height=6)
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

smoking_history_prostate = read_csv(file="../res/etc/MSK_VP_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						   type_convert() %>%
						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
smoking_history_breast = read_csv(file="../res/etc/MSK_VB_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						 type_convert() %>%
						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
data = left_join(data, clinical %>% select(patient_id, smoking_history = have_you_ever_smoked_), by="patient_id") %>%
	   type_convert() %>%
	   filter(smoking_history=="Yes")
p0 = zeroinfl(num_called ~ age, dist = "poisson", data = data)
tmp = data %>%
	  mutate(subj_type = ifelse(subj_type=="Control", "Control", "Cancer")) %>%
	  mutate(subj_type = as.factor(subj_type)) %>%
	  filter(num_called !=0) %>%
	  mutate(bio_source = "VUSo")
	   
plot.0 = ggplot(tmp , aes(x = age, y = num_called, fill = subj_type)) +
		 geom_point(alpha=.85, size=3.5, shape = 21, color = "#231F20") +
		 scale_fill_manual(values = c("Control"="#FDAE61", "Cancer"="salmon")) +
		 geom_smooth(formula = y ~ x, method="glm", aes(x = age, y = num_called), data=tmp, inherit.aes = FALSE, color="grey50", fill="grey50") +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.85), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\n Age (years)\n", y="Somatic cfDNA variants / Mb\n") +
 		 scale_y_log10(
 		 	breaks = function(x) { c(1, 2, 5, 10, 20, 30, 50) },
 		 	labels = function(x) { c("1", "2", "5", "10", "20", "30", "50") }
 		 ) +
		 coord_cartesian(xlim = c(20, 90), ylim = c(1, 50)) +
		 annotation_logticks(side="l") +
		 theme(legend.position="none") +
		 annotate(geom="text", x=55, y=50, label=paste0("p = ", signif(summary(p0)$coefficients$count[2,4], 3)), size=3.5)
	  	
pdf(file="../res/rebuttal/VUSo_vs_Age_Smoker.pdf", width=5, height=6)
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

smoking_history_prostate = read_csv(file="../res/etc/MSK_VP_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						   type_convert() %>%
						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
smoking_history_breast = read_csv(file="../res/etc/MSK_VB_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						 type_convert() %>%
						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
data = left_join(data, clinical %>% select(patient_id, smoking_history = have_you_ever_smoked_), by="patient_id") %>%
	   type_convert() %>%
	   filter(smoking_history=="Yes")
p0 = zeroinfl(num_called ~ age, dist = "poisson", data = data)
tmp = data %>%
	  mutate(subj_type = ifelse(subj_type=="Control", "Control", "Cancer")) %>%
	  mutate(subj_type = as.factor(subj_type)) %>%
	  filter(num_called !=0) %>%
	  mutate(bio_source = "Biopsy-matched")
	   
plot.0 = ggplot(tmp , aes(x = age, y = num_called, fill = subj_type)) +
		 geom_point(alpha=.85, size=3.5, shape = 21, color = "#231F20") +
		 scale_fill_manual(values = c("Control"="#FDAE61", "Cancer"="salmon")) +
		 geom_smooth(formula = y ~ x, method="glm", aes(x = age, y = num_called), data=tmp, inherit.aes = FALSE, color="grey50", fill="grey50") +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.85), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\n Age (years)\n", y="Somatic cfDNA variants / Mb\n") +
 		 scale_y_log10(
 		 	breaks = function(x) { c(1, 2, 5, 10, 20, 30, 50) },
 		 	labels = function(x) { c("1", "2", "5", "10", "20", "30", "50") }
 		 ) +
		 coord_cartesian(xlim = c(20, 90), ylim = c(1, 50)) +
		 annotation_logticks(side="l") +
		 theme(legend.position="none") +
		 annotate(geom="text", x=55, y=50, label=paste0("p = ", signif(summary(p0)$coefficients$count[2,4], 3)), size=3.5)
	  	
pdf(file="../res/rebuttal/Biopsy_matched_vs_Age_Smoker.pdf", width=5, height=6)
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

smoking_history_prostate = read_csv(file="../res/etc/MSK_VP_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						   type_convert() %>%
						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
smoking_history_breast = read_csv(file="../res/etc/MSK_VB_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						 type_convert() %>%
						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
data = left_join(data, clinical %>% select(patient_id, smoking_history = have_you_ever_smoked_), by="patient_id") %>%
	   type_convert() %>%
	   filter(smoking_history=="Yes")
p0 = zeroinfl(num_called ~ age, dist = "poisson", data = data)
tmp = data %>%
	  mutate(subj_type = ifelse(subj_type=="Control", "Control", "Cancer")) %>%
	  mutate(subj_type = as.factor(subj_type)) %>%
	  filter(num_called !=0) %>%
	  mutate(bio_source = "Biopsy-subrthreshold")
	   
plot.0 = ggplot(tmp , aes(x = age, y = num_called, fill = subj_type)) +
		 geom_point(alpha=.85, size=3.5, shape = 21, color = "#231F20") +
		 scale_fill_manual(values = c("Control"="#FDAE61", "Cancer"="salmon")) +
		 geom_smooth(formula = y ~ x, method="glm", aes(x = age, y = num_called), data=tmp, inherit.aes = FALSE, color="grey50", fill="grey50") +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.85), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\n Age (years)\n", y="Somatic cfDNA variants / Mb\n") +
 		 scale_y_log10(
 		 	breaks = function(x) { c(1, 2, 5, 10, 20, 30, 50) },
 		 	labels = function(x) { c("1", "2", "5", "10", "20", "30", "50") }
 		 ) +
		 coord_cartesian(xlim = c(20, 90), ylim = c(1, 50)) +
		 annotation_logticks(side="l") +
		 theme(legend.position="none") +
		 annotate(geom="text", x=55, y=50, label=paste0("p = ", signif(summary(p0)$coefficients$count[2,4], 3)), size=3.5)
	  	
pdf(file="../res/rebuttal/Biopsy_subthreshold_vs_Age_Smoker.pdf", width=5, height=6)
print(plot.0)
dev.off()

#==================================================
# Scatter plot of mutation burden versus age for
# biopsy matched, VUSo and WBC matched mutations
# for never smokers
#==================================================
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

smoking_history_prostate = read_csv(file="../res/etc/MSK_VP_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						   type_convert() %>%
						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
smoking_history_breast = read_csv(file="../res/etc/MSK_VB_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						 type_convert() %>%
						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
data = left_join(data, clinical %>% select(patient_id, smoking_history = have_you_ever_smoked_), by="patient_id") %>%
	   type_convert() %>%
	   filter(smoking_history=="No")
p0 = zeroinfl(num_called ~ age, dist = "poisson", data = data)
tmp = data %>%
	  mutate(subj_type = ifelse(subj_type=="Control", "Control", "Cancer")) %>%
	  mutate(subj_type = as.factor(subj_type)) %>%
	  filter(num_called !=0) %>%
	  mutate(bio_source = "WBC-matched")
	   
plot.0 = ggplot(tmp , aes(x = age, y = num_called, fill = subj_type)) +
		 geom_point(alpha=.85, size=3.5, shape = 21, color = "#231F20") +
		 scale_fill_manual(values = c("Control"="#FDAE61", "Cancer"="salmon")) +
		 geom_smooth(formula = y ~ x, method="glm", aes(x = age, y = num_called), data=tmp, inherit.aes = FALSE, color="grey50", fill="grey50") +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.85), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\n Age (years)\n", y="Somatic cfDNA variants / Mb\n") +
 		 scale_y_log10(
 		 	breaks = function(x) { c(1, 2, 5, 10, 20, 30, 50) },
 		 	labels = function(x) { c("1", "2", "5", "10", "20", "30", "50") }
 		 ) +
		 coord_cartesian(xlim = c(20, 90), ylim = c(1, 50)) +
		 annotation_logticks(side="l") +
		 theme(legend.position="none") +
		 annotate(geom="text", x=55, y=50, label=paste0("p = ", signif(summary(p0)$coefficients$count[2,4], 3)), size=3.5)
	  	
pdf(file="../res/rebuttal/WBC_matched_vs_Age_Non-Smoker.pdf", width=5, height=6)
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

smoking_history_prostate = read_csv(file="../res/etc/MSK_VP_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						   type_convert() %>%
						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
smoking_history_breast = read_csv(file="../res/etc/MSK_VB_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						 type_convert() %>%
						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
data = left_join(data, clinical %>% select(patient_id, smoking_history = have_you_ever_smoked_), by="patient_id") %>%
	   type_convert() %>%
	   filter(smoking_history=="No")
p0 = zeroinfl(num_called ~ age, dist = "poisson", data = data)
tmp = data %>%
	  mutate(subj_type = ifelse(subj_type=="Control", "Control", "Cancer")) %>%
	  mutate(subj_type = as.factor(subj_type)) %>%
	  filter(num_called !=0) %>%
	  mutate(bio_source = "VUSo")
	   
plot.0 = ggplot(tmp , aes(x = age, y = num_called, fill = subj_type)) +
		 geom_point(alpha=.85, size=3.5, shape = 21, color = "#231F20") +
		 scale_fill_manual(values = c("Control"="#FDAE61", "Cancer"="salmon")) +
		 geom_smooth(formula = y ~ x, method="glm", aes(x = age, y = num_called), data=tmp, inherit.aes = FALSE, color="grey50", fill="grey50") +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.85), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\n Age (years)\n", y="Somatic cfDNA variants / Mb\n") +
 		 scale_y_log10(
 		 	breaks = function(x) { c(1, 2, 5, 10, 20, 30, 50) },
 		 	labels = function(x) { c("1", "2", "5", "10", "20", "30", "50") }
 		 ) +
		 coord_cartesian(xlim = c(20, 90), ylim = c(1, 50)) +
		 annotation_logticks(side="l") +
		 theme(legend.position="none") +
		 annotate(geom="text", x=55, y=50, label=paste0("p = ", signif(summary(p0)$coefficients$count[2,4], 3)), size=3.5)
	  	
pdf(file="../res/rebuttal/VUSo_vs_Age_Non-Smoker.pdf", width=5, height=6)
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

smoking_history_prostate = read_csv(file="../res/etc/MSK_VP_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						   type_convert() %>%
						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
smoking_history_breast = read_csv(file="../res/etc/MSK_VB_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						 type_convert() %>%
						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
data = left_join(data, clinical %>% select(patient_id, smoking_history = have_you_ever_smoked_), by="patient_id") %>%
	   type_convert() %>%
	   filter(smoking_history=="No")
p0 = zeroinfl(num_called ~ age, dist = "poisson", data = data)
tmp = data %>%
	  mutate(subj_type = ifelse(subj_type=="Control", "Control", "Cancer")) %>%
	  mutate(subj_type = as.factor(subj_type)) %>%
	  filter(num_called !=0) %>%
	  mutate(bio_source = "Biopsy-matched")
	   
plot.0 = ggplot(tmp , aes(x = age, y = num_called, fill = subj_type)) +
		 geom_point(alpha=.85, size=3.5, shape = 21, color = "#231F20") +
		 scale_fill_manual(values = c("Control"="#FDAE61", "Cancer"="salmon")) +
		 geom_smooth(formula = y ~ x, method="glm", aes(x = age, y = num_called), data=tmp, inherit.aes = FALSE, color="grey50", fill="grey50") +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.85), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\n Age (years)\n", y="Somatic cfDNA variants / Mb\n") +
 		 scale_y_log10(
 		 	breaks = function(x) { c(1, 2, 5, 10, 20, 30, 50) },
 		 	labels = function(x) { c("1", "2", "5", "10", "20", "30", "50") }
 		 ) +
		 coord_cartesian(xlim = c(20, 90), ylim = c(1, 50)) +
		 annotation_logticks(side="l") +
		 theme(legend.position="none") +
		 annotate(geom="text", x=55, y=50, label=paste0("p = ", signif(summary(p0)$coefficients$count[2,4], 3)), size=3.5)
	  	
pdf(file="../res/rebuttal/Biopsy_matched_vs_Age_Non-Smoker.pdf", width=5, height=6)
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

smoking_history_prostate = read_csv(file="../res/etc/MSK_VP_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						   type_convert() %>%
						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
smoking_history_breast = read_csv(file="../res/etc/MSK_VB_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
						 type_convert() %>%
						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
data = left_join(data, clinical %>% select(patient_id, smoking_history = have_you_ever_smoked_), by="patient_id") %>%
	   type_convert() %>%
	   filter(smoking_history=="No")
p0 = zeroinfl(num_called ~ age, dist = "poisson", data = data)
tmp = data %>%
	  mutate(subj_type = ifelse(subj_type=="Control", "Control", "Cancer")) %>%
	  mutate(subj_type = as.factor(subj_type)) %>%
	  filter(num_called !=0) %>%
	  mutate(bio_source = "Biopsy-subrthreshold")
	   
plot.0 = ggplot(tmp , aes(x = age, y = num_called, fill = subj_type)) +
		 geom_point(alpha=.85, size=3.5, shape = 21, color = "#231F20") +
		 scale_fill_manual(values = c("Control"="#FDAE61", "Cancer"="salmon")) +
		 geom_smooth(formula = y ~ x, method="glm", aes(x = age, y = num_called), data=tmp, inherit.aes = FALSE, color="grey50", fill="grey50") +
		 facet_wrap(~bio_source) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.85), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\n Age (years)\n", y="Somatic cfDNA variants / Mb\n") +
 		 scale_y_log10(
 		 	breaks = function(x) { c(1, 2, 5, 10, 20, 30, 50) },
 		 	labels = function(x) { c("1", "2", "5", "10", "20", "30", "50") }
 		 ) +
		 coord_cartesian(xlim = c(20, 90), ylim = c(1, 50)) +
		 annotation_logticks(side="l") +
		 theme(legend.position="none") +
		 annotate(geom="text", x=55, y=50, label=paste0("p = ", signif(summary(p0)$coefficients$count[2,4], 3)), size=3.5)
	  	
pdf(file="../res/rebuttal/Biopsy_subthreshold_vs_Age_Non-Smoker.pdf", width=5, height=6)
print(plot.0)
dev.off()

# load(all_vars_and_clinical)
# save_vars = all_vars
# all_vars = all_vars %>%
#  		   filter(is_patient_valid) %>%
#  		   filter(c_panel) %>%
#  		   filter(!is_hypermutator) %>%
#  		   filter(!is_lowdepth) %>%
#  		   filter(!is_lowqual) %>%
#  		   filter(!is_tumor_matched) %>%
#  		   filter(!is_cfdna_matched)
#  
# n_samples = all_vars %>%
#  			distinct(patient_id) %>%
#   		    count()
# recurrence = all_vars %>%
#    			 group_by(loc_lng) %>%
#    			 count() %>%
#    			 ungroup() %>%
#    			 rename(n_recurrence=n) %>%
#    			 mutate(f_recurrence=n_recurrence/n_samples$n)
# all_vars = all_vars %>%
#    		   left_join(recurrence)
# all_vars = all_vars %>%
#  		   filter(f_recurrence < 0.05 | is_hotspot | (Variant_Classification %in% c("Frame_Shift_Ins", "Frame_Shift_Del", "Nonsense_Mutation", "Nonstop_Mutation") & (SYMBOL %in% chip_genes)))
#   		   
# recurrence = all_vars %>%
#    			 group_by(patient_id, loc_srt) %>%
#    			 count() %>%
#    			 ungroup() %>%
#    			 rename(n_indel=n)
# all_vars = all_vars %>%
# 		   left_join(recurrence)
# all_vars = all_vars %>%
# 		   filter(!(n_indel > 1 & indel))
# 		   
# all_vars = all_vars %>%
#  		   filter(Variant_Classification!="3'Flank") %>%
#   		   filter(Variant_Classification!="3'UTR") %>%
#   		   filter(Variant_Classification!="5'Flank") %>%
#   		   filter(Variant_Classification!="5'UTR") %>%
#   		   filter(Variant_Classification!="In_Frame_Del") %>%
#   		   filter(Variant_Classification!="In_Frame_Ins") %>%
#   		   filter(Variant_Classification!="Intron") %>%
#   		   filter(Variant_Classification!="RNA") %>%
#   		   filter(Variant_Classification!="Silent") %>%
#   		   filter(Variant_Classification!="IGR") %>%
#   		   filter(Variant_Classification!="Translation_Start_Site")
#   		   
# all_vars = all_vars %>%
#  		   filter(SYMBOL!="HLA-A")
# 
# all_vars = all_vars %>%
#  		   filter((adnobaq/dpnobaq)<=.3 | (Variant_Classification %in% c("Frame_Shift_Ins", "Frame_Shift_Del", "Nonsense_Mutation", "Nonstop_Mutation") & SYMBOL %in% chip_genes))
# 		   
# all_vars = all_vars %>%
#  		   filter(!in_exac)
#  		   
# all_vars = all_vars %>%
#  		   filter(!in_gnomad)
#  		   
# snvs = read_tsv(somatic_snvs_grail$scored, col_types = cols(.default = col_character()))  %>%
# 	   mutate(ID=paste0(patient_id, ":", chrom, ":", position, ":", ref_orig, ">", alt_orig))
# indels = read_tsv(somatic_indels_grail$scored, col_types = cols(.default = col_character()))  %>%
#  		 mutate(ID=paste0(patient_id, ":", chrom, ":", position, ":", ref_orig, ">", alt_orig))
# feature_names = intersect(colnames(indels), colnames(snvs))
# som_vars_grail = bind_rows(indels[,feature_names], snvs[,feature_names])
# all_vars = all_vars %>%
# 		   mutate(is_cfdna_matched = ID_x %in% som_vars_grail$ID)
# 
# burden_healthy = all_vars %>%
#   		   		 filter(subj_type=="Control") %>%
#   		   		 group_by(patient_id) %>%
#   	 	   		 summarize(num_called = n()) %>%
#   	 	   		 ungroup() %>%
#   	 	   		 left_join(clinical)
# patient_ids = save_vars %>%
# 			  filter(is_patient_valid) %>%
# 			  filter(subj_type=="Control") %>%
# 			  select(patient_id) %>%
# 			  distinct() %>%
# 			  .[["patient_id"]]
# index = (clinical$patient_id %in% patient_ids) & !(clinical$patient_id %in% burden_healthy$patient_id)  	 	   		 
# tmp = clinical[index,,drop=FALSE]
# tmp = cbind(tmp, num_called=rep(0, sum(index)))
# burden_healthy = rbind(burden_healthy, tmp[,colnames(burden_healthy)])
# burden_healthy_incfdna = all_vars %>%
#   		   		 filter(subj_type=="Control") %>%
#   		   		 filter(is_cfdna_matched) %>%
#   		   		 group_by(patient_id) %>%
#   	 	   		 summarize(num_called = n()) %>%
#   	 	   		 ungroup() %>%
#   	 	   		 left_join(clinical)
# index = (burden_healthy$patient_id %in% burden_healthy_incfdna$patient_id)
# tmp = burden_healthy[!index,,drop=FALSE]
# tmp$num_called = 0
# burden_healthy_incfdna = rbind(burden_healthy_incfdna, tmp)
# burden_healthy_ch = all_vars %>%
# 				    filter(SYMBOL %in% chip_genes) %>%
#   		   		    filter(subj_type=="Control") %>%
#   		   		    group_by(patient_id) %>%
#   	 	   		    summarize(num_called = n()) %>%
#   	 	   		    ungroup() %>%
#   	 	   		    left_join(clinical)
# index = (clinical$patient_id %in% patient_ids) & !(clinical$patient_id %in% burden_healthy_ch$patient_id)  
# tmp = clinical[index,,drop=FALSE]
# tmp = cbind(tmp, num_called=rep(0, sum(index)))
# burden_healthy_ch = rbind(burden_healthy_ch, tmp[,colnames(burden_healthy_ch)])
# burden_cancer = all_vars %>%
#   		   		filter(subj_type!="Control") %>%
#   		   		group_by(patient_id) %>%
#   	 	   		summarize(num_called = n()) %>%
#   	 	   		ungroup() %>%
#   	 	   		left_join(clinical)
# patient_ids = save_vars %>%
# 			  filter(is_patient_valid) %>%
# 			  filter(subj_type!="Control") %>%
# 			  select(patient_id) %>%
# 			  distinct() %>%
# 			  .[["patient_id"]]			  
# index = (clinical$patient_id %in% patient_ids) & !(clinical$patient_id %in% burden_cancer$patient_id)
# tmp = clinical[index,,drop=FALSE]
# tmp = cbind(tmp, num_called=rep(0, sum(index)))
# burden_cancer = rbind(burden_cancer, tmp[,colnames(burden_cancer)])
# burden_cancer_incfdna = all_vars %>%
#   		   		 filter(subj_type!="Control") %>%
#   		   		 filter(is_cfdna_matched) %>%
#   		   		 group_by(patient_id) %>%
#   	 	   		 summarize(num_called = n()) %>%
#   	 	   		 ungroup() %>%
#   	 	   		 left_join(clinical)
# index = (burden_cancer$patient_id %in% burden_cancer_incfdna$patient_id)
# tmp = burden_cancer[!index,,drop=FALSE]
# tmp$num_called = 0
# burden_cancer_incfdna = rbind(burden_cancer_incfdna, tmp)
# burden_cancer_ch = all_vars %>%
# 				   filter(SYMBOL %in% chip_genes) %>%
#   		   		   filter(subj_type!="Control") %>%
#   		   		   group_by(patient_id) %>%
#   	 	   		   summarize(num_called = n()) %>%
#   	 	   		   ungroup() %>%
#   	 	   		   left_join(clinical)
# index = (clinical$patient_id %in% patient_ids) & !(clinical$patient_id %in% burden_cancer_ch$patient_id)
# tmp = clinical[index,,drop=FALSE]
# tmp = cbind(tmp, num_called=rep(0, sum(index)))
# burden_cancer_ch = rbind(burden_cancer_ch, tmp[,colnames(burden_cancer_ch)])
# 
# data_.cancer_all = data.frame(burden_cancer[,c("age", "num_called", "study"),drop=FALSE], row.names=burden_cancer$patient_id)
# data_.healthy_all = data.frame(burden_healthy[,c("age", "num_called", "study"),drop=FALSE], row.names=burden_healthy$patient_id)
# data_.cancer_cfdna = data.frame(burden_cancer_incfdna[,c("age", "num_called", "study"),drop=FALSE], row.names=burden_cancer_incfdna$patient_id)
# data_.healthy_cfdna = data.frame(burden_healthy_incfdna[,c("age", "num_called", "study"),drop=FALSE], row.names=burden_healthy_incfdna$patient_id)
# data_.cancer_ch = data.frame(burden_cancer_ch[,c("age", "num_called", "study"),drop=FALSE], row.names=burden_cancer_ch$patient_id)
# data_.healthy_ch = data.frame(burden_healthy_ch[,c("age", "num_called", "study"),drop=FALSE], row.names=burden_healthy_ch$patient_id)
# data_.cancer_ch = cbind(data_.cancer_ch, matrix(NA, nrow=nrow(data_.cancer_ch), ncol=length(chip_genes), dimnames=list(rownames(data_.cancer_ch), chip_genes)))
# data_.healthy_ch = cbind(data_.healthy_ch, matrix(NA, nrow=nrow(data_.healthy_ch), ncol=length(chip_genes), dimnames=list(rownames(data_.healthy_ch), chip_genes)))
# for (i in 1:nrow(data_.cancer_ch)) {
# 	for (j in 1:length(chip_genes)) {
#  		index = all_vars$patient_id==rownames(data_.cancer_ch)[i] & all_vars$SYMBOL==chip_genes[j]
#  		if (sum(index)!=0) {
#  			data_.cancer_ch[i,chip_genes[j]] = max(100*all_vars[index,"adnobaq"]/all_vars[index,"dpnobaq"])
#  		}
#  	}
# }
# for (i in 1:nrow(data_.healthy_ch)) {
# 	for (j in 1:length(chip_genes)) {
#  		index = all_vars$patient_id==rownames(data_.healthy_ch)[i] & all_vars$SYMBOL==chip_genes[j]
#  		if (sum(index)!=0) {
#  			data_.healthy_ch[i,chip_genes[j]] = max(100*all_vars[index,"adnobaq"]/all_vars[index,"dpnobaq"])
#  		}
#  	}
# }
# 
# data_.cancer_cfdna = data_.cancer_cfdna[order(data_.cancer_cfdna$num_called),,drop=FALSE]
# data_.healthy_cfdna = data_.healthy_cfdna[order(data_.healthy_cfdna$num_called),,drop=FALSE]
# data_.cancer_all = data_.cancer_all[rownames(data_.cancer_cfdna),,drop=FALSE]
# data_.healthy_all = data_.healthy_all[rownames(data_.healthy_cfdna),,drop=FALSE]
# data_.cancer_ch = data_.cancer_ch[rownames(data_.cancer_cfdna),,drop=FALSE]
# data_.healthy_ch = data_.healthy_ch[rownames(data_.healthy_cfdna),,drop=FALSE]
# 
# data_.cancer_all = data_.cancer_all[order(data_.cancer_all$num_called),,drop=FALSE]
# data_.healthy_all = data_.healthy_all[order(data_.healthy_all$num_called),,drop=FALSE]
# data_.cancer_cfdna = data_.cancer_cfdna[rownames(data_.cancer_all),,drop=FALSE]
# data_.healthy_cfdna = data_.healthy_cfdna[rownames(data_.healthy_all),,drop=FALSE]
# data_.cancer_ch = data_.cancer_ch[rownames(data_.cancer_all),,drop=FALSE]
# data_.healthy_ch = data_.healthy_ch[rownames(data_.healthy_all),,drop=FALSE]
# 
# age_cat = list(start = c(20, 40, 50, 60, 70, 80),
# 			   end	 = c(40, 50, 60, 70, 80, 100))
# 
# index_.cancer = list()
# for (i in 1:6) {
# 	indx = which(data_.cancer_all$age>=age_cat$start[i] & data_.cancer_all$age<age_cat$end[i])
# 	index_.cancer[[i]] = indx
# }
# index_.healthy = list()
# for (i in 1:6) {
# 	indx = which(data_.healthy_all$age>=age_cat$start[i] & data_.healthy_all$age<age_cat$end[i])
# 	index_.healthy[[i]] = indx
# }
# data_all = data_cfdna = data_ch = list()
# for (i in 1:6) {
# 	data_all[[i]] = rbind(data_.healthy_all[index_.healthy[[i]],], data_.cancer_all[index_.cancer[[i]],])
# 	data_cfdna[[i]] = rbind(data_.healthy_cfdna[index_.healthy[[i]],], data_.cancer_cfdna[index_.cancer[[i]],])
# 	data_ch[[i]] = rbind(data_.healthy_ch[index_.healthy[[i]],], data_.cancer_ch[index_.cancer[[i]],])
# }
# data_all = do.call(rbind, data_all)
# data_cfdna = do.call(rbind, data_cfdna)
# data_cfdna = data_cfdna[rownames(data_all),,drop=FALSE]
# data_ch = do.call(rbind, data_ch)
# data_ch = data_ch[rownames(data_all),,drop=FALSE]
# ix = NULL
# for (i in 1:6) {
#  	indx = which(data_all$age>=age_cat$start[i] & data_all$age<age_cat$end[i])
#  	ix = c(ix, length(indx))
# }
# ix = cumsum(ix)
# colss = unlist(lapply(c("#F7F7F7", "#D9D9D9", "#BDBDBD", "#969696", "#636363", "#252525"), transparentRgb, 225))
# colss2 = c(rep("black", 4), rep("white", 2))
# max_vaf = vector(mode="numeric", length=nrow(data_all))
# for (i in 1:nrow(data_all)) {
#  	if (sum(all_vars$patient_id==rownames(data_all)[i])!=0) {
#  		max_vaf[i] = max(100*(all_vars$adnobaq/all_vars$dpnobaq)[all_vars$patient_id==rownames(data_all)[i]])
# 	} else {
# 		max_vaf[i] = 0.05
# 	}
# }
# to_plot = c("DNMT3A", "TP53", "TET2", "ASXL1", "PPM1D")
# to_plot_also = c("JAK2", "RUNX1", "SF3B1", "SRSF2", "IDH1", "IDH2", "U2AF1", "CBL", "ATM", "CHEK2")
# cols = c("#01985C", "#F7DC02", "#1E4665", "#FF7175", "#8CDB5E", "#3D98D3")
# 
# index_control = grepl("W", rownames(data_ch), fixed=TRUE)
# index_late = grepl("V", rownames(data_ch), fixed=TRUE)
# tmp = data_ch[index_control,-c(1:3),drop=FALSE]
# tmp  = cbind(tmp[,to_plot,drop=FALSE],
# 			 "OTHER_CH"=apply(tmp[,to_plot_also], 1, max, na.rm=TRUE))
# tmp[is.infinite(as.matrix(tmp))] = NA
# tmp[is.na(tmp)] = 0
# tmp = cbind(tmp, "NO_CH"=rep(0, nrow(tmp)))
# 
# for (i in 1:nrow(tmp)) {
# 	if (sum(tmp[i,])==0) {
# 		tmp[i,"NO_CH"] = 99
# 	} else {
# 		index = which.max(tmp[i,])
# 		tmp[i,index] = 99
# 	}
# }
# tmp[tmp!=99] = 0
# tmp[tmp==99] = 1
# tmp = data.frame(tmp, "subj_type" = "Control")
# 
# tmp2 = data_ch[index_late,-c(1:3),drop=FALSE]
# tmp2 = cbind(tmp2[,to_plot,drop=FALSE],
#  			 "OTHER_CH"=apply(tmp2[,to_plot_also], 1, max, na.rm=TRUE))
# tmp2[is.infinite(as.matrix(tmp2))] = NA
# tmp2[is.na(tmp2)] = 0
# tmp2 = cbind(tmp2, "NO_CH"=rep(0, nrow(tmp2)))
# for (i in 1:nrow(tmp2)) {
# 	if (sum(tmp2[i,])==0) {
# 		tmp2[i,"NO_CH"] = 99
# 	} else {
# 		index = which.max(tmp2[i,])
# 		tmp2[i,index] = 99
# 	}
# }
# tmp2[tmp2!=99] = 0
# tmp2[tmp2==99] = 1
# tmp2 = data.frame(tmp2, "subj_type" = "Cancer")
# 
# data = rbind(tmp, tmp2)
# data = cbind(data, "patient_id" = rownames(data))
# 
# tx = read.csv(file="../res/etc/prior_tx_techval_0818.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
# 	 rename(patient_id = ID)
# data = full_join(data, tx, by="patient_id")
# smoking_history_prostate = read_csv(file="../res/etc/MSK_VP_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
# 						   type_convert() %>%
# 						   rename(patient_id = `Patient ID`, have_you_ever_smoked_ = TobaccoHx) %>%
# 						   mutate(have_you_ever_smoked_ = ifelse(have_you_ever_smoked_=="Y", "Yes", "No"))
# smoking_history_breast = read_csv(file="../res/etc/MSK_VB_SMOKING.csv", col_types = cols(.default = col_character()))  %>%
# 						 type_convert() %>%
# 						 rename(patient_id = `MSK Sample ID`, have_you_ever_smoked_ = TobaccoHx)
# smoking_history_updated = bind_rows(smoking_history_prostate, smoking_history_breast)
# clinical = left_join(clinical, smoking_history_updated, by="patient_id") %>%
# 		   mutate(have_you_ever_smoked_ = ifelse(is.na(have_you_ever_smoked_.x), have_you_ever_smoked_.y, have_you_ever_smoked_.x))
# data = left_join(data, clinical %>% select(patient_id, smoking_history = have_you_ever_smoked_), by="patient_id") %>%
# 	   type_convert() %>%
# 	   mutate(prior_rt = ifelse(subj_type=="Control", 0, prior_rt)) %>%
# 	   mutate(prior_chemo = ifelse(subj_type=="Control", 0, prior_chemo)) %>%
# 	   mutate(prior_ctrt = prior_rt==1 | prior_chemo==1)
# 	   
# mylogit <- glm(prior_ctrt ~ factor(PPM1D) + smoking_history, data = data, family=binomial("logit"))
