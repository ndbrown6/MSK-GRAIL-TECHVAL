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
# Permutation based p-value
#==================================================
all_vars = read_tsv(file=wbc_scored_annotated_and_clinical$scored, col_types = cols(.default = col_character())) %>%
		   type_convert()
clinical = read_tsv(file=wbc_scored_annotated_and_clinical$clinical, col_types = cols(.default = col_character())) %>%
 		   type_convert()
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
# Variant class filter
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
			  dplyr::select(patient_id) %>%
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
			  dplyr::select(patient_id) %>%
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

max_vaf = vector(mode="numeric", length=nrow(data_all))
for (i in 1:nrow(data_all)) {
 	if (sum(all_vars$patient_id==rownames(data_all)[i])!=0) {
 		max_vaf[i] = max(100*(all_vars$adnobaq/all_vars$dpnobaq)[all_vars$patient_id==rownames(data_all)[i]])
	} else {
		max_vaf[i] = 0.05
	}
}
to_plot = c("DNMT3A", "TP53", "TET2", "ASXL1", "PPM1D")
to_plot_also = c("JAK2", "RUNX1", "SF3B1", "SRSF2", "IDH1", "IDH2", "U2AF1", "CBL", "ATM", "CHEK2")

index_control = grepl("W", rownames(data_ch), fixed=TRUE)
index_late = grepl("V", rownames(data_ch), fixed=TRUE)
tmp = data_ch[index_control,-c(1:3),drop=FALSE]
tmp  = cbind(tmp[,to_plot,drop=FALSE],
			 "OTHER_CH"=apply(tmp[,to_plot_also], 1, max, na.rm=TRUE))
tmp[is.infinite(as.matrix(tmp))] = NA
tmp[is.na(as.matrix(tmp))] = 0
tmp = cbind(tmp, "NO_CH"=rep(0, nrow(tmp)))
for (i in 1:nrow(tmp)) {
	if (sum(tmp[i,])==0) {
		tmp[i,"NO_CH"] = 1
	} else {
		index = which.max(tmp[i,])
		tmp[i,] = 0
		tmp[i,index] = 1
	}
}
x_ = tmp %>%
	 mutate(subj_type = 0) %>%
	 mutate(patient_id = rownames(tmp))
	 
tmp = data_ch[index_late,-c(1:3),drop=FALSE]
tmp  = cbind(tmp[,to_plot,drop=FALSE],
 			 "OTHER_CH"=apply(tmp[,to_plot_also], 1, max, na.rm=TRUE))
tmp[is.infinite(as.matrix(tmp))] = NA
tmp[is.na(as.matrix(tmp))] = 0
tmp = cbind(tmp, "NO_CH"=rep(0, nrow(tmp)))
for (i in 1:nrow(tmp)) {
	if (sum(tmp[i,])==0) {
		tmp[i,"NO_CH"] = 1
	} else {
		index = which.max(tmp[i,])
		tmp[i,] = 0
		tmp[i,index] = 1
	}
}
y_ = tmp %>%
	 mutate(subj_type = 1) %>%
	 mutate(patient_id = rownames(tmp))
	 
data = bind_rows(x_, y_)
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
data = left_join(data, clinical %>% dplyr::select(patient_id, smoking_history = have_you_ever_smoked_, age), by="patient_id") %>%
	   type_convert() %>%
	   mutate(smoking_history = as.factor(smoking_history))

tx = read.csv(file=url_prior_tx, header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
	 rename(patient_id = ID) %>%
	 mutate(prior_rt_ct = ifelse(prior_rt==1 | prior_chemo==1, 1, 0))

data = left_join(data, tx %>% dplyr::select(patient_id, prior_rt_ct)) %>%
	   mutate(prior_rt_ct = ifelse(is.na(prior_rt_ct) & subj_type==0, 0, prior_rt_ct))

fisher.test(data$PPM1D[data$smoking_history=="Yes"], data$prior_rt_ct[data$smoking_history=="Yes"])

