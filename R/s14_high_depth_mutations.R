#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/figureS14")) {
	dir.create("../res/figureS14")
}

#==================================================
# heatmap of ch-related variants in wbc
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

cols = c("Tumor-matched" = as.character(variant_cols["biopsy_matched"]),
		 "WBC-matched"	 = as.character(variant_cols["WBC_matched"]),
		 "VUSo"			 = as.character(variant_cols["VUSo"]))

plot.0 = ggplot(tmp.1, aes(x=reorder(patient_id, -N), y=n)) +
  		 geom_bar(stat="identity", aes(fill=`Variant category`)) +
  		 theme(axis.text.x=element_text(angle=90)) +
  		 scale_fill_manual(values = cols) +
  		 coord_cartesian(ylim = c(-.1, 120)) +
  		 labs(x="", y="Number of mutations > 10,000X\n") +
  		 guides(fill=guide_legend(title=c("Variant category")))

pdf(file="../res/figureS14/Number_High_Depth_Mutations_by_Patient.pdf", width=10, height=7)
par(mar=c(6.1, 6.5, 4.1, 1.1))
h = as.matrix(table(tmp.1$patient_id, tmp.1$`Variant category`))
for (i in 1:nrow(h)) {
	for (j in 1:ncol(h)) {
		if (h[i,j]==1) {
			sample_name = rownames(h)[i]
			variant_category = colnames(h)[j]
			id = which(tmp.1$patient_id == sample_name & tmp.1$`Variant category` == variant_category)
			h[i,j] = as.numeric(tmp.1[id,"n"])
		}
	}
}
index = order(apply(h, 1, sum), decreasing=TRUE)
h = h[index,,drop=FALSE]
x = barplot(t(h[,3:1,drop=FALSE]), space=.15, beside=FALSE, col=cols[rev(colnames(h))], las=2, ylab="", axes=FALSE)
axis(1, at=x, labels=rep("",length(x)), cex.axis=1.5, padj=0.25, line=.25, tcl=-.35)
axis(2, at=NULL, cex.axis=1.35, las=1, lwd=1.75, lwd.ticks=1.55, line=-.5)
mtext(side=2, text="Number of mutations", line=4, cex=1.5)
legend("topright", pch=22, legend=names(cols), col="black", box.lwd=-1, pt.bg=cols, cex=1, pt.cex=1.65)
dev.off()

tmp = variants %>%
	  filter(bio_source %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
	  mutate(facets = "Collapsed variant level coverage")
	  
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
		 #scale_fill_manual(values = cols) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.25, 0.8), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="Depth in cfDNA\n", x="\nMean collapsed coverage\n") +
		 scale_y_log10(
 		 ) + 
		 annotation_logticks() +
		 guides(fill=FALSE) +
		 facet_wrap(~bio_source)

pdf(file="../res/figureS14/High_Depth_Mutations_by_Mean_Coverage.pdf", width=11, height=11)
par(mar=c(6.1, 6.5, 4.1, 1.1))
z = split.screen(figs=matrix(c(0+.05,.5,.5+.05,1,
							   .5,1-.05,.5+.05,1,
							   0+.05,.5,0+.05,.5,
							   .5,1-.05,0+.05,.5,
							   0, 1, 0, 1), nrow=5, ncol=4, byrow=TRUE))
screen(z[1])
x = tmp$Collapsed_Mean_Coverage[tmp$bio_source=="biopsy_matched"]/10e2
y = tmp$dpnobaq[tmp$bio_source=="biopsy_matched"]
plot(x, y, pch=21, bg=variant_cols["biopsy_matched"], col="black", cex=1.5, axes=FALSE, frame.plot=FALSE, main="", xlab="", ylab="", log="y", xlim=c(0,8), ylim=c(1,100000))
axis(1, at=NULL, cex.axis=1.5, padj=0.25, lwd=1.75, lwd.ticks=1.55)
axis(2, at=c(1, 10, 100, 1000, 10000, 100000), labels=c(1, expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)), cex.axis=1.5, las=1, lwd=1.75, lwd.ticks=1.55)
log10_axis(side=2, at=c(1, 10, 100, 1000, 10000, 100000), lwd=0, lwd.ticks=1)
screen(z[2])
x = tmp$Collapsed_Mean_Coverage[tmp$bio_source=="IMPACT-BAM_matched"]/10e2
y = tmp$dpnobaq[tmp$bio_source=="IMPACT-BAM_matched"]
plot(x, y, pch=21, bg=variant_cols["IMPACT-BAM_matched"], col="black", cex=1.5, axes=FALSE, frame.plot=FALSE, main="", xlab="", ylab="", log="y", xlim=c(0,8), ylim=c(1,100000))
axis(1, at=NULL, cex.axis=1.5, padj=0.25, lwd=1.75, lwd.ticks=1.55)
axis(2, at=c(1, 10, 100, 1000, 10000, 100000), labels=c(1, expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)), cex.axis=1.5, las=1, lwd=1.75, lwd.ticks=1.55)
log10_axis(side=2, at=c(1, 10, 100, 1000, 10000, 100000), lwd=0, lwd.ticks=1)
screen(z[3])
x = tmp$Collapsed_Mean_Coverage[tmp$bio_source=="VUSo"]/10e2
y = tmp$dpnobaq[tmp$bio_source=="VUSo"]
plot(x, y, pch=21, bg=variant_cols["VUSo"], col="black", cex=1.5, axes=FALSE, frame.plot=FALSE, main="", xlab="", ylab="", log="y", xlim=c(0,8), ylim=c(1,100000))
axis(1, at=NULL, cex.axis=1.5, padj=0.25, lwd=1.75, lwd.ticks=1.55)
axis(2, at=c(1, 10, 100, 1000, 10000, 100000), labels=c(1, expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)), cex.axis=1.5, las=1, lwd=1.75, lwd.ticks=1.55)
log10_axis(side=2, at=c(1, 10, 100, 1000, 10000, 100000), lwd=0, lwd.ticks=1)
screen(z[4])
x = tmp$Collapsed_Mean_Coverage[tmp$bio_source=="WBC_matched"]/10e2
y = tmp$dpnobaq[tmp$bio_source=="WBC_matched"]
plot(x, y, pch=21, bg=variant_cols["WBC_matched"], col="black", cex=1.5, axes=FALSE, frame.plot=FALSE, main="", xlab="", ylab="", log="y", xlim=c(0,8), ylim=c(1,100000))
axis(1, at=NULL, cex.axis=1.5, padj=0.25, lwd=1.75, lwd.ticks=1.55)
axis(2, at=c(1, 10, 100, 1000, 10000, 100000), labels=c(1, expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)), cex.axis=1.5, las=1, lwd=1.75, lwd.ticks=1.55)
log10_axis(side=2, at=c(1, 10, 100, 1000, 10000, 100000), lwd=0, lwd.ticks=1)
screen(z[5])
plot(0, 0, type="n", axes=FALSE, frame.plot=FALSE, xlab="", ylab="")
mtext(side=1, text=expression("Mean collapsed depth ("%.%10^3~")"), line=3, cex=1.5)
mtext(side=2, text="Collapsed variant depth in cfDNA", line=3, cex=1.5)
close.screen(all.screens=TRUE)
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

cols = c("tumor_matched" = as.character(variant_cols["biopsy_matched"]),
		 "WBC_matched"	 = as.character(variant_cols["WBC_matched"]),
		 "VUSo"			 = as.character(variant_cols["VUSo"]))
plot.0 = ggplot(tmp, aes(y = dpnobaq, x = cat_2, fill = bio_source)) +
		 geom_violin(trim=FALSE) +
		 scale_fill_manual(values = cols) +
		 stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
		 labs(x="\nInput DNA for library preparation (ng)", y="Depth in cfDNA\n") +
		 scale_y_log10(
		 	breaks = function(x) { c(100, 1000, 10000, 100000) },
		 	labels = function(x) { c(expression(10^2), expression(10^3), expression(10^4), expression(10^5)) }
 		 ) +
 		 annotation_logticks(sides="l") +
 		 coord_cartesian(ylim = c(100,100000)) +
		 guides(fill=FALSE)

pdf(file="../res/figureS14/High_Depth_Mutations_by_input_DNA.pdf", width=10, height=6)
print(plot.0)
dev.off()

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

cols = c("Tumor matched" = as.character(variant_cols["biopsy_matched"]),
		 "WBC matched"	 = as.character(variant_cols["WBC_matched"]),
		 "VUSo"			 = as.character(variant_cols["VUSo"]))
pdf(file="../res/figureS14/High_Depth_Mutations_by_VAF.pdf", width=8.35, height=6)
par(mar=c(6.1, 6.5, 4.1, 1.1))
z = split.screen(figs=matrix(c(0,1/3+.1,0,1,
							   1/3-.05,2/3+.05,0,1,
							   2/3-.1,1,0,1,
							   0, 1, 0, 1), nrow=4, ncol=4, byrow=TRUE))
screen(z[1])
x = tmp$af_nobaq[tmp$bio_source=="Tumor matched"]
y = tmp$dpnobaq[tmp$bio_source=="Tumor matched"]
plot(x, y, pch=21, bg=cols["Tumor matched"], col="black", cex=1.5, axes=FALSE, frame.plot=FALSE, main="", xlab="", ylab="", log="xy", xlim=c(.01,100), ylim=c(1,100000))
axis(1, at=c(.01, .1, 1, 10, 100), labels=c(expression(10^-2), expression(10^-1), expression(1), expression(10^1), expression(10^2)), cex.axis=1.35, padj=0.25, lwd=1.75, lwd.ticks=1.55)
axis(2, at=c(1, 10, 100, 1000, 10000, 100000), labels=c(1, expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)), cex.axis=1.5, las=1, lwd=1.75, lwd.ticks=1.55)
log10_axis(side=2, at=c(1, 10, 100, 1000, 10000, 100000), lwd=0, lwd.ticks=1)
mtext(side=2, text="Collapsed variant depth in cfDNA", line=4, cex=1.5)
title(main="\nTumor matched")
screen(z[2])
x = tmp$af_nobaq[tmp$bio_source=="VUSo"]
y = tmp$dpnobaq[tmp$bio_source=="VUSo"]
plot(x, y, pch=21, bg=cols["VUSo"], col="black", cex=1.5, axes=FALSE, frame.plot=FALSE, main="", xlab="", ylab="", log="xy", xlim=c(.01,100), ylim=c(1,100000))
axis(1, at=c(.01, .1, 1, 10, 100), labels=c(expression(10^-2), expression(10^-1), expression(1), expression(10^1), expression(10^2)), cex.axis=1.35, padj=0.25, lwd=1.75, lwd.ticks=1.55)
axis(2, at=c(1, 10, 100, 1000, 10000, 100000), labels=rep("", 6), cex.axis=1.5, las=1, lwd=1.75, lwd.ticks=1.55)
log10_axis(side=2, at=c(1, 10, 100, 1000, 10000, 100000), lwd=0, lwd.ticks=1)
mtext(side=1, text="VAF (%)", line=4, cex=1.5)
title(main="\nVUSo")
screen(z[3])
x = tmp$af_nobaq[tmp$bio_source=="WBC matched"]
y = tmp$dpnobaq[tmp$bio_source=="WBC matched"]
plot(x, y, pch=21, bg=cols["WBC matched"], col="black", cex=1.5, axes=FALSE, frame.plot=FALSE, main="", xlab="", ylab="", log="xy", xlim=c(.01,100), ylim=c(1,100000))
axis(1, at=c(.01, .1, 1, 10, 100), labels=c(expression(10^-2), expression(10^-1), expression(1), expression(10^1), expression(10^2)), cex.axis=1.35, padj=0.25, lwd=1.75, lwd.ticks=1.55)
axis(2, at=c(1, 10, 100, 1000, 10000, 100000), labels=rep("", 6), cex.axis=1.5, las=1, lwd=1.75, lwd.ticks=1.55)
log10_axis(side=2, at=c(1, 10, 100, 1000, 10000, 100000), lwd=0, lwd.ticks=1)
title(main="\nWBC matched")
close.screen(all.screens=TRUE)
dev.off()
