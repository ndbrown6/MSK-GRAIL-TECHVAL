#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/rebuttal")) {
	dir.create("../res/rebuttal")
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
		   		  
patient_ids = unique(variants$patient_id)

tmp.0 = variants %>%
		filter(bio_source %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
		filter(dpnobaq < 10000) %>%
		group_by(patient_id) %>%
		count(patient_id) %>%
		full_join(data.frame(patient_id = patient_ids), by="patient_id") %>%
		mutate(n_low = ifelse(is.na(n), 0, n)) %>%
		dplyr::select(-n)
	
tmp.1 = variants %>%
		filter(bio_source %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
		filter(dpnobaq >= 10000) %>%
		group_by(patient_id) %>%
		count(patient_id) %>%
		full_join(data.frame(patient_id = patient_ids), by="patient_id") %>%
		mutate(n_high = ifelse(is.na(n), 0, n)) %>%
		dplyr::select(-n)
		
tmp = left_join(tmp.0, tmp.1, by="patient_id") %>%
	  mutate(Hypermutated = ifelse(patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id), "Yes", "No")) %>%
	  mutate(n_high = ifelse(n_high == 0, 0.5, n_high)) %>%
	  mutate(n_low = ifelse(n_low == 0, 0.5, n_low)) %>%
	  mutate(facets = "Somatic mutations per sample")

cols = c("Yes"="#D7191C", "No"="#2B83BA")
plot.0 = ggplot(tmp, aes(x = n_high, y = n_low)) +
		 geom_point(alpha=.85, size=3.5, pch = 21, colour = "black", aes(fill=Hypermutated)) +
		 scale_fill_manual(values = cols) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.8, 0.3), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="Number of mutations < 10,000X\n", x="\nNumber of mutations > 10,000X\n") +
		 scale_x_log10(
		 	breaks = function(x) { c(0.5, 1, 10, 100, 200) },
		 	labels = function(x) { c("0", "1", "10", "100", "200") }
		 ) + 
		 scale_y_log10(
		 	breaks = function(x) { c(0.5, 1, 10, 100, 700) },
		 	labels = function(x) { c("0", "1", "10", "100", "700") }
		 ) +
	 	coord_cartesian(xlim = c(0.5, 200), ylim = c(0.5,700)) +
		annotation_logticks() +
		guides(fill=guide_legend(title=c("Hypermutated"))) +
		facet_wrap(~facets)

pdf(file="../res/rebuttal/Number_Mutations_by_Patient_DP.pdf", width=5, height=6)
print(plot.0)
dev.off()

tmp.0 = variants %>%
		filter(dpnobaq>=10000) %>%
		filter(bio_source %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
		mutate(bio_source = ifelse(bio_source %in% c("biopsy_matched", "IMPACT-BAM_matched"), "Tumor-matched", bio_source)) %>%
		mutate(bio_source = ifelse(bio_source == "WBC_matched", "WBC-matched", bio_source)) %>%
		group_by(patient_id, bio_source) %>%
		count(patient_id, bio_source) %>%
		rename(`Variant category` = bio_source)
		
tmp.1 = variants %>%
		filter(dpnobaq>=10000) %>%
		filter(bio_source %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
		group_by(patient_id) %>%
		count(patient_id) %>%
		rename(N = n)
		
tmp.0 = tmp.0 %>%
		left_join(tmp.1)
		
patient_ids = unique(tmp.0$patient_id)

tmp.1 = NULL
for (i in c("Tumor-matched", "WBC-matched", "VUSo")) {
	tmp.1 = bind_rows(tmp.1, subset(tmp.0, `Variant category`==i) %>%
							 full_join(data.frame(patient_id = patient_ids), by="patient_id")) %>%
			mutate(n = ifelse(is.na(n), 0, n))
}
tmp.1 = tmp.1 %>%
		 filter(!is.na(`Variant category`))

cols = c("Tumor-matched"="#D7191C",
		 "WBC-matched"="#ABDDA4",
		 "VUSo"="#2B83BA")

plot.0 = ggplot(tmp.1, aes(x=reorder(patient_id, -N), y=n)) +
  		 geom_bar(stat="identity", aes(fill=`Variant category`)) +
  		 theme(axis.text.x=element_text(angle=90)) +
  		 scale_fill_manual(values = cols) +
  		 coord_cartesian(ylim = c(-.1, 120)) +
  		 labs(x="", y="Number of mutations > 10,000X\n") +
  		 guides(fill=guide_legend(title=c("Variant category")))

pdf(file="../res/rebuttal/Number_High_Depth_Mutations_by_Patient.pdf", width=10, height=6)
print(plot.0)
dev.off()

tmp.0 = variants %>%
		filter(bio_source %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
		filter(dpnobaq >= 10000) %>%
		group_by(gene) %>%
		count(gene)
		
bed = read_tsv(file=common_bed_annotated_w_introns, col_names = c("chrom", "start", "end", "gene", "intron"))
target_lengths = vector(mode="numeric", length=nrow(tmp.0))
intron_sizes = vector(mode="numeric", length=length(tmp.0))
for (i in 1:nrow(tmp.0)) {
 	index = bed$gene==tmp.0$gene[i]
 	target_lengths[i] = sum(bed[index,3]-bed[index,2])
 	intron_sizes[i] = t(as.matrix(bed[index,3]-bed[index,2])) %*% as.matrix(bed[index,5,drop=FALSE])
}
index = target_lengths == intron_sizes
intron_sizes[index] = 0
target_lengths = target_lengths - intron_sizes

tmp.0 = bind_cols(tmp.0, data.frame(size=target_lengths, facets="Number of mutations per gene"))

plot.0 = ggplot(tmp.0, aes(x = n, y = size, label=gene)) +
		 geom_point(alpha=.85, size=3.5, pch = 21, colour = "black", fill = "#2B83BA") +
	   	 geom_text_repel(data = subset(tmp.0, n >= 4), size=3, force=1, fontface=3)+
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.3, 0.8), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\nNumber of mutations > 10,000X\n", y="Coding target size (bp)\n") +
		 scale_y_log10(
		 	breaks = function(x) { c(50, 100, 1000, 10000, 20000) },
		 	labels = function(x) { c("50", "100", "1000", "10000", "20000") }
		 ) + 
		 scale_x_continuous(
		 	breaks = function(x) { c(0, 2, 4, 6, 8, 10) },
		 	labels = function(x) { c("0", "2", "4", "6", "8", "10") }
		 ) +
	 	 coord_cartesian(xlim = c(0, 10), ylim = c(50,20000)) +
		 annotation_logticks(sides = "l") +
		 facet_wrap(~facets)
		 
pdf(file="../res/rebuttal/High_Depth_Mutations_by_Target_Size.pdf", width=5, height=6)
print(plot.0)
dev.off()

tmp = variants %>%
	  filter(bio_source %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
	  mutate(facets = "Collapsed variant level coverage")
	  
cols = c("biopsy_matched"="#D7191C",
		 "IMPACT-BAM_matched"="#FDAE61",
		 "WBC_matched"="#ABDDA4",
		 "VUSo"="#2B83BA")

plot.0 = ggplot(tmp, aes(x = dpnobaq, y = dpgdna)) +
		 geom_point(alpha=1, size=3.5, pch = 21, colour = "black", aes(fill=bio_source)) +
		 scale_fill_manual(values = cols) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.8, 0.8), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="Variant depth in WBC\n", x="\nVariant depth in cfDNA\n") +
 	 	 coord_cartesian(xlim = c(1, 25000), ylim = c(1,25000)) +
		 guides(fill=guide_legend(title=c("Variant category"))) +
		 facet_wrap(~facets)

pdf(file="../res/rebuttal/High_Depth_cfDNA_vs_gDNA.pdf", width=6, height=6)
print(plot.0)
dev.off()

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

qc_metrics_cfdna = read.csv(file="../modified_v11/QC_metrics/TechVal_Merlin_QC_metrics.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
			 	   dplyr::select(sample_id, patient_id, sample_type, tissue, volume_of_blood_mL, volume_of_DNA_source_mL, DNA_extraction_yield_ng, DNA_input_concentration_ng_uL, Library_preparation_input_ng, raw.MEAN_BAIT_COVERAGE, collapsed.MEAN_BAIT_COVERAGE, collapsed_fragment.MEAN_BAIT_COVERAGE, readErrorRate, readSubstErrorRate, Study) %>%
			 	   filter(sample_type=="cfDNA")
tracker_grail_cfdna = read.csv(file=patient_tracker, header=TRUE, sep=",", stringsAsFactors=FALSE) %>%
					  dplyr::select(patient_id, cfdna_sample_id) %>%
					  rename(msk_id = patient_id, sample_id = cfdna_sample_id)
qc_metrics_cfdna = left_join(qc_metrics_cfdna, tracker_grail_cfdna, by="sample_id")

qc_metrics_wbc = read.csv(file="../modified_v11/QC_metrics/TechVal_Merlin_QC_metrics.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
			 	 dplyr::select(sample_id, patient_id, sample_type, tissue, volume_of_blood_mL, volume_of_DNA_source_mL, DNA_extraction_yield_ng, DNA_input_concentration_ng_uL, Library_preparation_input_ng, raw.MEAN_BAIT_COVERAGE, collapsed.MEAN_BAIT_COVERAGE, collapsed_fragment.MEAN_BAIT_COVERAGE, readErrorRate, readSubstErrorRate, Study) %>%
			 	 filter(sample_type=="gDNA")
tracker_grail_wbc = read.csv(file=patient_tracker, header=TRUE, sep=",", stringsAsFactors=FALSE) %>%
					dplyr::select(patient_id, gdna_sample_id) %>%
					rename(msk_id = patient_id, sample_id = gdna_sample_id)
qc_metrics_wbc = left_join(qc_metrics_wbc, tracker_grail_wbc, by="sample_id")

qc_metrics = rbind(qc_metrics_cfdna, qc_metrics_wbc) %>%
			 filter(msk_id %in% valid_patient_ids) %>%
			 dplyr::select(-sample_id, -msk_id) %>%
			 rename(Patient_ID = patient_id,
			 		Sample_Type = sample_type,
			 		Tissue = tissue,
			 		Volume_of_blood_mL = volume_of_blood_mL,
			 		Volume_of_DNA_source_mL = volume_of_DNA_source_mL,
			 		Uncollapsed_Mean_Coverage = raw.MEAN_BAIT_COVERAGE,
				   	Collapsed_Mean_Coverage = collapsed.MEAN_BAIT_COVERAGE,
				   	Collapsed_Fragment_Mean_Coverage = collapsed_fragment.MEAN_BAIT_COVERAGE,
				   	Indel_and_Substitution_Error_Rate = readErrorRate,
				   	Substitution_Error_Rate = readSubstErrorRate,
				   	Assay_Version = Study) %>%
			mutate(Assay_Version = ifelse(Assay_Version=="TechVal", "V1", "V2")) %>%
			arrange(Patient_ID, Sample_Type) %>%
			mutate(Library_preparation_input_ng = ifelse(Library_preparation_input_ng>75, 75, Library_preparation_input_ng))
	  
tmp = left_join(tmp,
				qc_metrics %>%
				filter(Sample_Type=="cfDNA") %>%
				rename(patient_id = Patient_ID),
				by = "patient_id")
				
plot.0 = ggplot(tmp, aes(y = dpnobaq, x = Collapsed_Mean_Coverage)) +
		 geom_point(alpha=1, size=2, pch = 21, colour = "black", aes(fill=bio_source)) +
		 scale_fill_manual(values = cols) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.25, 0.8), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="Depth in cfDNA\n", x="\nMean collapsed coverage\n") +
		 scale_y_log10(
 		 ) + 
		 annotation_logticks() +
		 guides(fill=FALSE) +
		 facet_wrap(~bio_source)

pdf(file="../res/rebuttal/High_Depth_Mutations_by_Mean_Coverage.pdf", width=6, height=6)
print(plot.0)
dev.off()

tmp = tmp %>%
	  mutate(bio_source = ifelse(bio_source %in% c("biopsy_matched","IMPACT-BAM_matched"), "tumor_matched", bio_source)) %>%
	  mutate(cat_1 = case_when(
	 	   		Library_preparation_input_ng<30 ~ "<30",
		   		Library_preparation_input_ng>=30 & Library_preparation_input_ng<45 ~ "30-44",
		   		Library_preparation_input_ng>=45 & Library_preparation_input_ng<60 ~ "45-59",
		   		Library_preparation_input_ng>=60 & Library_preparation_input_ng<75 ~ "60-74",
		   		Library_preparation_input_ng>=75 ~ "75",
		   		)) %>%
	  mutate(cat_2 = factor(paste0(bio_source, ":", cat_1), levels=paste0(rep(c("tumor_matched", "VUSo", "WBC_matched"), each=5), ":", c("<30", "30-44", "45-59", "60-74", "75")))) %>%
	  mutate(dpnobaq = ifelse(dpnobaq>59000, NA, dpnobaq)) %>%
  	  mutate(dpnobaq = ifelse(dpnobaq<100, NA, dpnobaq))

plot.0 = ggplot(tmp, aes(y = dpnobaq, x = cat_2)) +
		 geom_violin(trim=FALSE) +
		 stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
		 facet_wrap(~facets) + 
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="\nInput DNA for library\npreparation (ng)", y="Depth in cfDNA\n") +
		 scale_y_log10(
		 	breaks = function(x) { c(100, 1000, 10000, 50000) },
		 	labels = function(x) { c(100, 1000, 10000, 50000) }
 		 ) +
 		 coord_cartesian(ylim = c(100,60000)) +
		 guides(fill=FALSE)

pdf(file="../res/rebuttal/High_Depth_Mutations_by_input_DNA.pdf", width=10, height=6)
print(plot.0)
dev.off()

cols = c("Tumor matched"="#D7191C",
		 "WBC matched"="#ABDDA4",
		 "VUSo"="#2B83BA")

tmp = variants %>%
	  filter(bio_source %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
	  mutate(bio_source = ifelse(bio_source %in% c("biopsy_matched", "IMPACT-BAM_matched"), "tumor_matched", bio_source)) %>%
	  mutate(bio_source = case_when(
	  	bio_source == "tumor_matched" ~ "Tumor matched",
	  	bio_source == "WBC_matched" ~ "WBC matched",
	  	bio_source == "VUSo" ~ "VUSo"))
		
plot.0 = ggplot(tmp, aes(y = dpnobaq, x = af_nobaq)) +
		 geom_point(alpha=1, size=3, pch = 21, colour = "black", aes(fill=bio_source)) +
		 scale_fill_manual(values = cols) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.25, 0.8), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="Collapsed variant depth in cfDNA\n", x="\nVAF (%)\n") +
		 scale_x_log10(
		 	breaks = function(x) { c(.1, 1, 10, 100) },
		 	labels = function(x) { c(".1", "1", "10", "100") }
 		 ) + 
		 scale_y_log10(
 		 ) + 
		 annotation_logticks(side="bl") +
		 guides(fill=FALSE) +
		 facet_wrap(~bio_source)

pdf(file="../res/rebuttal/High_Depth_Mutations_by_VAF.pdf", width=10, height=6)
print(plot.0)
dev.off()

tmp = variants %>%
	  filter(bio_source %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
	  mutate(bio_source = ifelse(bio_source %in% c("biopsy_matched", "IMPACT-BAM_matched"), "tumor_matched", bio_source)) %>%
	  filter(patient_id == "MSK-VB-0023") %>%
	  mutate(chrom = gsub("chr", "", chrom, fixed=TRUE)) %>%
	  mutate(chrom = ifelse(chrom=="X", "23", chrom)) %>%
	  mutate(chrom = as.numeric(chrom))
	  
cols = c("tumor_matched"="#D7191C",
		 "WBC_matched"="#ABDDA4",
		 "VUSo"="#2B83BA")

pdf(file="../res/rebuttal/MSK-VB-0023_Manhattan-Plot_High_DP.pdf", width=10, height=2.25)
par(mar=c(5, 5, 4, 2)+.1)
data(CytoBand)
end = NULL
for (i in 1:23) {
	end = c(end, max(CytoBand[CytoBand[,1]==i,"End"]))
}
end = cumsum(end)
start = c(1, end[1:22]+1)
CytoBand = cbind(start, end)
index = NULL
for (i in 1:23) {
	index = c(index, seq(from = CytoBand[i, "start"], to=CytoBand[i, "end"], length=sum(tmp$chrom==i)))
}
indx = tmp$dpnobaq>=10000
plot(0, 0, type="n", axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-.07,.1), xlim=c(1,max(index)))
points(index[indx], rep(0, sum(indx)), type="p", pch="|", cex=2, col=cols[tmp$bio_source[indx]])
abline(v=1, col=transparent_rgb("goldenrod3", 255), lty=3, lwd=.5)
for (j in 1:23) {
	abline(v=CytoBand[j,"end"], col=transparent_rgb("goldenrod3", 255), lty=3, lwd=.5)
}
axis(1, at = .5*(CytoBand[,"start"]+CytoBand[,"end"]), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
dev.off()

pdf(file="../res/rebuttal/MSK-VB-0023_Manhattan-Plot_Low_DP.pdf", width=10, height=2.25)
par(mar=c(5, 5, 4, 2)+.1)
data(CytoBand)
end = NULL
for (i in 1:23) {
	end = c(end, max(CytoBand[CytoBand[,1]==i,"End"]))
}
end = cumsum(end)
start = c(1, end[1:22]+1)
CytoBand = cbind(start, end)
index = NULL
for (i in 1:23) {
	index = c(index, seq(from = CytoBand[i, "start"], to=CytoBand[i, "end"], length=sum(tmp$chrom==i)))
}
indx = tmp$dpnobaq<10000
plot(0, 0, type="n", axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-.07,.1), xlim=c(1,max(index)))
points(index[indx], rep(0, sum(indx)), type="p", pch="|", cex=2, col=cols[tmp$bio_source[indx]])
abline(v=1, col=transparent_rgb("goldenrod3", 255), lty=3, lwd=.5)
for (j in 1:23) {
	abline(v=CytoBand[j,"end"], col=transparent_rgb("goldenrod3", 255), lty=3, lwd=.5)
}
axis(1, at = .5*(CytoBand[,"start"]+CytoBand[,"end"]), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
dev.off()


#==================================================
# Scatter plot of VAF of replicates for MKS-VB-0023
# of mutations at DP > 100000X
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

patient_ids = c("MSK-VB-0023")
vars_rep0 = variants %>%
			filter(patient_id %in% patient_ids)

gdna_params = data_frame(
				subj_type = c("Healthy", "Breast", "Lung", "Prostate"),
				min_p = c(0.8, 0.79, 0.82, 0.79))
clean_target_region = read_tsv(url_target.bed, col_names = c("chrom", "start", "end"))
clean_target_region$in_target = TRUE
techval_repeats = read.csv(url_techval.repeats, header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
				  filter(!patient_id %in% c("MSK-VB-0023_B", "MSK-VL-0028_B", "MSK-VL-0042_B")) %>%
				  mutate(patient_id = gsub(pattern="_A", replacement="", patient_id)) %>%
				  mutate_at(vars(matches("qual")), funs(if_else(is.na(.), 0, .))) %>%
				  mutate(pedge = ifelse(is.na(pedge), 0, pedge), pgtkxgdna = ifelse(is.na(pgtkxgdna), 0, pgtkxgdna)) %>%
				  separate(hgvsp, into = c("p1", "p2", "p3", "p4"), sep = "\\|") %>%
				  separate(is_nonsyn, into = c("t1","t2", "t3", "t4"), sep = "\\|") %>%
				  separate(symbol, into = c("g1", "g2", "g3", "g4"), sep = "\\|") %>%
				  mutate(gene = case_when(
									(.$t1 == TRUE | is.na(.$t2)) ~g1,
									(.$t1 != TRUE & .$t2 == TRUE) ~g2,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 == TRUE) ~g3,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 != TRUE & .$t4 == TRUE) ~g4)) %>%
				  mutate(hgvs_p = case_when(
									(.$t1 == TRUE) ~p1,
									(.$t1 != TRUE & .$t2 == TRUE) ~p2,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 == TRUE) ~p3,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 != TRUE & .$t4 == TRUE) ~p4)) %>%
				  mutate(hgvs_p = sub(".*:", "", hgvs_p),
						 is_nonsyn = ifelse(is.na(hgvs_p), FALSE, TRUE)) %>%
				  dplyr::select(-t1, -t2, -t3, -t4, -g1, -g2, -g3, -g4, -p1, -p2, -p3, -p4) %>%
				  mutate(filter = "",
					  		 subj_type = case_when(
					  		 				grepl("VB", .$patient_id) ~ "Breast",
					  		 				grepl("VL", .$patient_id) ~ "Lung",
					  		 				grepl("VP", .$patient_id) ~ "Prostate")) %>%
				  left_join(gdna_params) %>%
				  mutate_(filter = ~update_filter(
					  	  filter = filter,
					  	  qualnobaq = qualnobaq,
					  	  pgtkxgdna = pgtkxgdna,
					  	  is_edge = isedge,
					  	  min_p = min_p)) %>%
				  dplyr::select(-min_p) %>%
				  mutate(loc = str_c(chrom, ":", pos, "_", ref, ">", alt)) %>%
				  mutate(position_orig = pos, ref_orig = ref, alt_orig = alt, position = pos) %>%
				  mutate(gratio = (adnobaq+2)*(dpgdna+4)/((adgdna+2)*(dpnobaq+4)),
				  		 gzero = adgdna/sqrt(adgdna+2)) %>%
				  mutate(start = position_orig,
                         end = position_orig + 1) %>%
                  genome_left_join(clean_target_region, by = c("chrom", "start", "end")) %>%
                  mutate(chrom = chrom.x) %>%
                  dplyr::select(-c(chrom.x, chrom.y, start.x, end.x, start.y, end.y)) %>%
				  left_join(variants %>% dplyr::select(patient_id, chrom, position_orig, ref_orig, alt_orig, MSK, grail), by=c("patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
				  mutate(MSK = ifelse(is.na(MSK), 0, MSK)) %>%
				  mutate(grail = 1) %>%
				  mutate(study = "TechVal") %>%
				  filter(in_target)
				  
feature_names = intersect(colnames(small_vars_plasma), colnames(techval_repeats))
repeat_variants = bind_rows(small_vars_plasma[,feature_names,drop=FALSE] %>%
							filter(!(patient_id %in% patient_ids)),
							techval_repeats[,feature_names,drop=FALSE])

repeat_variants = label_bio_source(repeat_variants)
repeat_variants = left_join(repeat_variants, msk_anno %>% dplyr::select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
repeat_variants = repeat_variants %>%
		   				mutate(bio_source = case_when(
		   					   MSK == 1 & grail == 1 ~ "biopsy_matched",
		   					   MSK == 1 & grail == 0 ~ "biopsy_only",
		   					   category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
		   					   category %in% c("germline", "germlineish") ~ "germline",
		   					   category %in% c("blood", "bloodier") ~ "WBC_matched",
		   					   category == "somatic" & `IM-T.alt_count` > bam_read_count_cutoff ~ "IMPACT-BAM_matched",
		   					   category == "somatic" ~ "VUSo",
		   					   TRUE ~ "other"),
			   		  	af_nobaq = round(adnobaq / dpnobaq * 100, 2),
			   		  	af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq))

vars_rep1 = repeat_variants %>%
			filter(patient_id %in% patient_ids)
			
techval_repeats = read.csv(url_techval.repeats, header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
				  filter(patient_id %in% c("MSK-VB-0023_B", "MSK-VL-0028_B", "MSK-VL-0042_B")) %>%
				  mutate(patient_id = gsub(pattern="_B", replacement="", patient_id)) %>%
				  mutate_at(vars(matches("qual")), funs(if_else(is.na(.), 0, .))) %>%
				  mutate(pedge = ifelse(is.na(pedge), 0, pedge), pgtkxgdna = ifelse(is.na(pgtkxgdna), 0, pgtkxgdna)) %>%
				  separate(hgvsp, into = c("p1", "p2", "p3", "p4"), sep = "\\|") %>%
				  separate(is_nonsyn, into = c("t1","t2", "t3", "t4"), sep = "\\|") %>%
				  separate(symbol, into = c("g1", "g2", "g3", "g4"), sep = "\\|") %>%
				  mutate(gene = case_when(
									(.$t1 == TRUE | is.na(.$t2)) ~g1,
									(.$t1 != TRUE & .$t2 == TRUE) ~g2,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 == TRUE) ~g3,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 != TRUE & .$t4 == TRUE) ~g4)) %>%
				  mutate(hgvs_p = case_when(
									(.$t1 == TRUE) ~p1,
									(.$t1 != TRUE & .$t2 == TRUE) ~p2,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 == TRUE) ~p3,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 != TRUE & .$t4 == TRUE) ~p4)) %>%
				  mutate(hgvs_p = sub(".*:", "", hgvs_p),
						 is_nonsyn = ifelse(is.na(hgvs_p), FALSE, TRUE)) %>%
				  dplyr::select(-t1, -t2, -t3, -t4, -g1, -g2, -g3, -g4, -p1, -p2, -p3, -p4) %>%
				  mutate(filter = "",
					  		 subj_type = case_when(
					  		 				grepl("VB", .$patient_id) ~ "Breast",
					  		 				grepl("VL", .$patient_id) ~ "Lung",
					  		 				grepl("VP", .$patient_id) ~ "Prostate")) %>%
				  left_join(gdna_params) %>%
				  mutate_(filter = ~update_filter(
					  	  filter = filter,
					  	  qualnobaq = qualnobaq,
					  	  pgtkxgdna = pgtkxgdna,
					  	  is_edge = isedge,
					  	  min_p = min_p)) %>%
				  dplyr::select(-min_p) %>%
				  mutate(loc = str_c(chrom, ":", pos, "_", ref, ">", alt)) %>%
				  mutate(position_orig = pos, ref_orig = ref, alt_orig = alt, position = pos) %>%
				  mutate(gratio = (adnobaq+2)*(dpgdna+4)/((adgdna+2)*(dpnobaq+4)),
				  		 gzero = adgdna/sqrt(adgdna+2)) %>%
				  mutate(start = position_orig,
                         end = position_orig + 1) %>%
                  genome_left_join(clean_target_region, by = c("chrom", "start", "end")) %>%
                  mutate(chrom = chrom.x) %>%
                  dplyr::select(-c(chrom.x, chrom.y, start.x, end.x, start.y, end.y)) %>%
				  left_join(variants %>% dplyr::select(patient_id, chrom, position_orig, ref_orig, alt_orig, MSK, grail), by=c("patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
				  mutate(MSK = ifelse(is.na(MSK), 0, MSK)) %>%
				  mutate(grail = 1) %>%
				  mutate(study = "TechVal") %>%
				  filter(in_target)
				  
feature_names = intersect(colnames(small_vars_plasma), colnames(techval_repeats))
repeat_variants = bind_rows(small_vars_plasma[,feature_names,drop=FALSE] %>%
							filter(!(patient_id %in% patient_ids)),
							techval_repeats[,feature_names,drop=FALSE])

repeat_variants = label_bio_source(repeat_variants)
repeat_variants = left_join(repeat_variants, msk_anno %>% dplyr::select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
repeat_variants = repeat_variants %>%
		   				mutate(bio_source = case_when(
		   					   MSK == 1 & grail == 1 ~ "biopsy_matched",
		   					   MSK == 1 & grail == 0 ~ "biopsy_only",
		   					   category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
		   					   category %in% c("germline", "germlineish") ~ "germline",
		   					   category %in% c("blood", "bloodier") ~ "WBC_matched",
		   					   category == "somatic" & `IM-T.alt_count` > bam_read_count_cutoff ~ "IMPACT-BAM_matched",
		   					   category == "somatic" ~ "VUSo",
		   					   TRUE ~ "other"),
			   		  	af_nobaq = round(adnobaq / dpnobaq * 100, 2),
			   		  	af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq))

vars_rep2 = repeat_variants %>%
			filter(patient_id %in% patient_ids)
			

all_vars = full_join(vars_rep0 %>% mutate(replicate = 1),
					 vars_rep1 %>% mutate(replicate = 2), by=c("study", "subj_type", "patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
		   dplyr::select(patient_id, adnobaq.x, adnobaq.y, dpnobaq.x, dpnobaq.y, bio_source.x, bio_source.y) %>%
		   mutate(af_rep1 = 100*adnobaq.x/dpnobaq.x) %>%
		   mutate(af_rep2 = 100*adnobaq.y/dpnobaq.y) %>%
		   mutate(bio_source.x = ifelse(is.na(bio_source.x), "unmatched", bio_source.x)) %>%
		   mutate(bio_source.y = ifelse(is.na(bio_source.y), "unmatched", bio_source.y)) %>%
		   filter(bio_source.x %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo") | bio_source.y %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
		   mutate(adnobaq.x = ifelse(is.na(adnobaq.x), 0, adnobaq.x)) %>%
		   mutate(adnobaq.y = ifelse(is.na(adnobaq.y), 0, adnobaq.y)) %>%
  		   mutate(dpnobaq.x = ifelse(is.na(dpnobaq.x), 0, dpnobaq.x)) %>%
   		   mutate(dpnobaq.y = ifelse(is.na(dpnobaq.y), 0, dpnobaq.y)) %>%
   		   mutate(afnobaq.x = 100*adnobaq.x/dpnobaq.x) %>%
   		   mutate(afnobaq.x = ifelse(afnobaq.x==0 | is.na(afnobaq.x), 0.01, afnobaq.x)) %>%
   		   mutate(afnobaq.y = 100*adnobaq.y/dpnobaq.y) %>%
   		   mutate(afnobaq.y = ifelse(afnobaq.y==0 | is.na(afnobaq.y), 0.01, afnobaq.y))


cols = c(
		 "Not detected in one replicate"="#D7191C",
		 "Not called in one replicate due\nto low quality"="#FDAE61",
		 "Incorrect assignment between replicates"="#ABDDA4",
		 "Called in both replicates"="#2B83BA"
		 )
		 
tmp_vars = all_vars %>% filter(patient_id == patient_ids)

# fix biopsy_matched in y
index = tmp_vars$bio_source.x == "biopsy_matched" & (tmp_vars$bio_source.y !="biopsy_matched" & tmp_vars$bio_source.y != "unmatched")
if (sum(index)!=0) {
	tmp_vars$bio_source.y[index] = "biopsy_matched"
}

# fix biopsy_matched in y
index = tmp_vars$bio_source.x == "IMPACT-BAM_matched" & (tmp_vars$bio_source.y !="IMPACT-BAM_matched" & tmp_vars$bio_source.y != "unmatched")
if (sum(index)!=0) {
	tmp_vars$bio_source.y[index] = "IMPACT-BAM_matched"
}

tmp_vars = tmp_vars %>%
		   mutate(shape = ifelse(bio_source.x=="biopsy_matched" | bio_source.x=="IMPACT-BAM_matched" | bio_source.y=="biopsy_matched" | bio_source.y=="IMPACT-BAM_matched", "Biopsy matched", "Biopsy unmatched")) %>%
		   mutate(fill = ifelse(bio_source.x!=bio_source.y , "Incorrect assignment between replicates", "Called in both replicates")) %>%
		   mutate(fill = ifelse((bio_source.x=="unmatched" & bio_source.y!="unmatched") | (bio_source.x!="unmatched" & bio_source.y=="unmatched"), "Not detected in one replicate", "Called in both replicates")) %>%
		   mutate(fill = ifelse((bio_source.x=="noise" & (bio_source.y!="noise" & bio_source.y!="other")) | (bio_source.y=="noise" & (bio_source.x!="noise" & bio_source.x!="other")), "Not called in one replicate due\nto low quality", fill)) %>%
		   filter(dpnobaq.x>10000) %>%
		   filter(bio_source.x!="noise")
		   

plot.0 = ggplot(tmp_vars, aes(x = afnobaq.x, y = afnobaq.y, shape = shape, fill = fill)) +
		 geom_abline(linetype = 1, color = "goldenrod3") +
		 geom_point(alpha=1, size=3.5) +
		 scale_fill_manual(values = cols) +
		 scale_shape_manual(values = c(24, 21)) +
		 facet_wrap(~patient_id) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.3, 0.8), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\nReplicate 1 (%)\n", y="Replicate 2 (%)\n") +
		 scale_x_log10(
			breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
			labels = function(x) { c("0", "0.1", "1", "10", "100") }
		 ) + 
		 scale_y_log10(
			breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
			labels = function(x) { c("0", "0.1", "1", "10", "100") }
		 ) +
		 coord_cartesian(xlim = c(0.01, 100), ylim = c(0.01,100)) +
		 annotation_logticks() +
		 guides(shape=guide_legend(title=c("Biopsy concordance"), override.aes=list(fill="black"))) +
		 guides(fill=guide_legend(title=c("Variant category")))
		 
pdf(file=paste0("../res/rebuttal/", patient_ids, "_R1_R2_High_Depth.pdf"), width=5.5, height=6.5)
print(plot.0)
dev.off()

all_vars = full_join(vars_rep0 %>% mutate(replicate = 1),
					 vars_rep2 %>% mutate(replicate = 3), by=c("study", "subj_type", "patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
		   dplyr::select(patient_id, adnobaq.x, adnobaq.y, dpnobaq.x, dpnobaq.y, bio_source.x, bio_source.y) %>%
		   mutate(af_rep1 = 100*adnobaq.x/dpnobaq.x) %>%
		   mutate(af_rep2 = 100*adnobaq.y/dpnobaq.y) %>%
		   mutate(bio_source.x = ifelse(is.na(bio_source.x), "unmatched", bio_source.x)) %>%
		   mutate(bio_source.y = ifelse(is.na(bio_source.y), "unmatched", bio_source.y)) %>%
		   filter(bio_source.x %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo") | bio_source.y %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
		   mutate(adnobaq.x = ifelse(is.na(adnobaq.x), 0, adnobaq.x)) %>%
		   mutate(adnobaq.y = ifelse(is.na(adnobaq.y), 0, adnobaq.y)) %>%
  		   mutate(dpnobaq.x = ifelse(is.na(dpnobaq.x), 0, dpnobaq.x)) %>%
   		   mutate(dpnobaq.y = ifelse(is.na(dpnobaq.y), 0, dpnobaq.y)) %>%
   		   mutate(afnobaq.x = 100*adnobaq.x/dpnobaq.x) %>%
   		   mutate(afnobaq.x = ifelse(afnobaq.x==0 | is.na(afnobaq.x), 0.01, afnobaq.x)) %>%
   		   mutate(afnobaq.y = 100*adnobaq.y/dpnobaq.y) %>%
   		   mutate(afnobaq.y = ifelse(afnobaq.y==0 | is.na(afnobaq.y), 0.01, afnobaq.y))


tmp_vars = all_vars %>% filter(patient_id == patient_ids)
	
# fix biopsy_matched in y
index = tmp_vars$bio_source.x == "biopsy_matched" & (tmp_vars$bio_source.y !="biopsy_matched" & tmp_vars$bio_source.y != "unmatched")
if (sum(index)!=0) {
	tmp_vars$bio_source.y[index] = "biopsy_matched"
}
	
# fix biopsy_matched in y
index = tmp_vars$bio_source.x == "IMPACT-BAM_matched" & (tmp_vars$bio_source.y !="IMPACT-BAM_matched" & tmp_vars$bio_source.y != "unmatched")
if (sum(index)!=0) {
	tmp_vars$bio_source.y[index] = "IMPACT-BAM_matched"
}
	
tmp_vars = tmp_vars %>%
		   mutate(shape = ifelse(bio_source.x=="biopsy_matched" | bio_source.x=="IMPACT-BAM_matched" | bio_source.y=="biopsy_matched" | bio_source.y=="IMPACT-BAM_matched", "Biopsy matched", "Biopsy unmatched")) %>%
		   mutate(fill = ifelse(bio_source.x!=bio_source.y , "Incorrect assignment between replicates", "Called in both replicates")) %>%
		   mutate(fill = ifelse((bio_source.x=="unmatched" & bio_source.y!="unmatched") | (bio_source.x!="unmatched" & bio_source.y=="unmatched"), "Not detected in one replicate", "Called in both replicates")) %>%
		   mutate(fill = ifelse((bio_source.x=="noise" & (bio_source.y!="noise" & bio_source.y!="other")) | (bio_source.y=="noise" & (bio_source.x!="noise" & bio_source.x!="other")), "Not called in one replicate due\nto low quality", fill)) %>%
		   filter(dpnobaq.x>10000) %>%
		   filter(bio_source.x!="noise")
			   	
plot.0 = ggplot(tmp_vars, aes(x = afnobaq.x, y = afnobaq.y, shape = shape, fill = fill)) +
		 geom_abline(linetype = 1, color = "goldenrod3") +
		 geom_point(alpha=1, size=3.5) +
		 scale_fill_manual(values = cols) +
		 scale_shape_manual(values = c(24, 21)) +
		 facet_wrap(~patient_id) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.3, 0.8), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\nReplicate 1 (%)\n", y="Replicate 3 (%)\n") +
		 scale_x_log10(
		 	breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
		 	labels = function(x) { c("0", "0.1", "1", "10", "100") }
		 ) + 
		 scale_y_log10(
		 	breaks = function(x) { c(0.01, 0.1, 1, 10, 100) },
		 	labels = function(x) { c("0", "0.1", "1", "10", "100") }
		 ) +
		 coord_cartesian(xlim = c(0.01, 100), ylim = c(0.01,100)) +
		 annotation_logticks() +
		 guides(shape=guide_legend(title=c("Biopsy concordance"), override.aes=list(fill="black"))) +
		 guides(fill=guide_legend(title=c("Variant category")))
	 
pdf(file=paste0("../res/rebuttal/", patient_ids, "_R1_R3_High_Depth.pdf"), width=5.5, height=6.5)
print(plot.0)
dev.off()
