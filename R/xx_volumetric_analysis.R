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
# Box plot for cfDNA fraction by number of
# metastatic sites
#==================================================
clinical = read_tsv(file=clinical_file_updated, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
cfdna_frac = read.csv(file=url_ctdna_frac, header=TRUE, sep=",", stringsAsFactors=FALSE) %>%
			 filter(!is.na(ctdna_frac)) %>%
			 mutate(index = order_samples(ID)) %>%
 			 arrange(desc(ctdna_frac)) %>%
			 arrange(index)
number_metastatic_sites = unlist(lapply(strsplit(clinical[,"metastatic_sites",drop=TRUE], split=",", fixed=TRUE), function(x) {length(x)}))
number_metastatic_sites[is.na(clinical[,"metastatic_sites",drop=TRUE])] = 0
names(number_metastatic_sites) = clinical[,"patient_id",drop=TRUE]
number_metastatic_sites = number_metastatic_sites[cfdna_frac[,1]]
cfdna_frac = cbind(cfdna_frac, number_metastatic_sites) %>%
			 mutate(Tissue = "") %>%
			 mutate(Tissue = ifelse(grepl("VB", ID), "Breast", Tissue)) %>%
			 mutate(Tissue = ifelse(grepl("VL", ID), "Lung", Tissue)) %>%
			 mutate(Tissue = ifelse(grepl("VP", ID), "Prostate", Tissue)) %>%
			 mutate(cat = case_when(
			 		number_metastatic_sites==1 | number_metastatic_sites==2 ~ 1,
			 		number_metastatic_sites==3 ~ 2,
			 		number_metastatic_sites>=4 ~ 3))

tmp.0 = cfdna_frac %>%
		filter(Tissue=="Breast")
		
p = signif(jonckheere.test(x=tmp.0$ctdna_frac, g=as.numeric(tmp.0$cat), alternative = "increasing")$p.value, 3)

plot.0 = ggplot(tmp.0, aes(x = factor(cat), y = ctdna_frac)) + 
		 geom_boxplot(alpha=1, outlier.size=NA, outlier.shape=NA, color="black", fill="white") +
		 geom_point(alpha=1, size=3.75, shape=21, color="black", fill="salmon", aes(x = jitter(as.numeric(cat)), y = ctdna_frac), data=tmp.0) +
		 facet_wrap(~Tissue) +
		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
 		 labs(x="", y="ctDNA fraction\n") +
 		 coord_cartesian(ylim = c(0,1)) +
		 annotate(geom="text", x=2, y=1, label=paste0("p = ", p), size=4.5)
		 			 
pdf(file="../res/rebuttal/cfDNA_Fraction_by_Ns_Breast.pdf", width=5.5, height=5)
print(plot.0)
dev.off()

tmp.0 = cfdna_frac %>%
		filter(Tissue=="Lung")
		
p = signif(jonckheere.test(x=tmp.0$ctdna_frac, g=as.numeric(tmp.0$cat), alternative = "increasing")$p.value, 3)

plot.0 = ggplot(tmp.0, aes(x = factor(cat), y = ctdna_frac)) + 
		 geom_boxplot(alpha=1, outlier.size=NA, outlier.shape=NA, color="black", fill="white") +
		 geom_point(alpha=1, size=3.75, shape=21, color="black", fill="#FDAE61", aes(x = jitter(as.numeric(cat)), y = ctdna_frac), data=tmp.0) +
		 facet_wrap(~Tissue) +
		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
 		 labs(x="", y="ctDNA fraction\n") +
 		 coord_cartesian(ylim = c(0,1)) +
		 annotate(geom="text", x=2, y=1, label=paste0("p = ", p), size=4.5)
		 			 
pdf(file="../res/rebuttal/cfDNA_Fraction_by_Ns_Lung.pdf", width=5.5, height=5)
print(plot.0)
dev.off()

tmp.0 = cfdna_frac %>%
		filter(Tissue=="Prostate")
		
p = signif(jonckheere.test(x=tmp.0$ctdna_frac, g=as.numeric(tmp.0$cat), alternative = "increasing")$p.value, 3)

plot.0 = ggplot(tmp.0, aes(x = factor(cat), y = ctdna_frac)) + 
		 geom_boxplot(alpha=1, outlier.size=NA, outlier.shape=NA, color="black", fill="white") +
		 geom_point(alpha=1, size=3.75, shape=21, color="black", fill="#ABDDA4", aes(x = jitter(as.numeric(cat)), y = ctdna_frac), data=tmp.0) +
		 facet_wrap(~Tissue) +
		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
 		 labs(x="", y="ctDNA fraction\n") +
 		 coord_cartesian(ylim = c(0,1)) +
		 annotate(geom="text", x=2, y=1, label=paste0("p = ", p), size=4.5)
		 			 
pdf(file="../res/rebuttal/cfDNA_Fraction_by_Ns_Prostate.pdf", width=5.5, height=5)
print(plot.0)
dev.off()

#==================================================
# Box plot for cfDNA fraction by disease
# volume
#==================================================
ctdna_fraction = read_csv(file=url_ctdna_frac, col_types = cols(.default = col_character()))  %>%
 				 type_convert() %>%
 				 rename(GRAIL_ID = ID)

vol_breast = read.xls(xls=url_volumetric_data, sheet=1, stringsAsFactors=FALSE) %>%
			 type_convert() %>%
			 mutate(Date_Tissue = as.character(Date_Tissue)) %>%
			 mutate(Date_Most_Recent_Scan = as.character(Date_Most_Recent_Scan)) %>%
			 mutate(Date_of_Last_Treatment = as.character(Date_of_Last_Treatment))

vol_lung = read.xls(xls=url_volumetric_data, sheet=2, stringsAsFactors=FALSE) %>%
		   type_convert() %>%
		   mutate(Date_Tissue = as.character(Date_Tissue)) %>%
		   mutate(Date_Most_Recent_Scan = as.character(Date_Most_Recent_Scan)) %>%
		   mutate(Date_of_Last_Treatment = as.character(Date_of_Last_Treatment)) %>%
		   mutate(Date_Dx = as.character(Date_Dx)) %>%
		   mutate(Date_Metastatic = as.character(Date_Metastatic)) %>%
		   mutate(Date_Tissue = as.character(Date_Tissue))
		   
vol_prostate = read.xls(xls=url_volumetric_data, sheet=3, stringsAsFactors=FALSE) %>%
		       type_convert() %>%
		       mutate(BSI_Value = ifelse(BSI_Value==">13", 14, BSI_Value)) %>%
		       mutate(BSI_Value = ifelse(BSI_Value=="n/a", NA, BSI_Value)) %>%
		       mutate(BSI_Value = as.numeric(BSI_Value)) %>%
		       rename(BSI_Volume = BSI_Value) %>%
		       mutate(Date_Dx = as.character(Date_Dx)) %>%
		       mutate(Date_Metastatic = as.character(Date_Metastatic)) %>%
		       mutate(Date_Tissue = as.character(Date_Tissue))

volumetric_data = bind_rows(vol_breast, vol_lung) %>%
				  bind_rows(vol_prostate) %>%
				  left_join(ctdna_fraction, by="GRAIL_ID")
				  
volumetric_data = volumetric_data %>%
				  mutate(tot_volu = apply( volumetric_data[,which(grepl('volume', colnames(volumetric_data), ignore.case = T))] , 1 , function(x) sum(as.numeric(x), na.rm=T) )) %>%
				  mutate(tot_mets = apply( volumetric_data[,which(!grepl('volume', colnames(volumetric_data), ignore.case = T) & grepl('Lesions|Lymph_Nodes|Pleural_Disease', colnames(volumetric_data), ignore.case = T))] , 1 , function(x) sum(as.numeric(x), na.rm=T))) %>%
				  filter(!(GRAIL_ID %in% c('MSK-VL-0057','MSK-VL-0003','MSK-VB-0046')))
				  
				  
# triangles = Unmeasurable disease not included in volumetric assessment
# circles = Measured in volumetric assessment

tmp = volumetric_data %>%
 	  mutate(ln_frac = log(ctdna_frac)) %>%
 	  mutate(ln_vol = log(tot_volu)) %>%
 	  mutate(facet = "Breast") %>%
 	  filter(grepl("VB", GRAIL_ID)) %>%
 	  mutate(Tissue = "Breast") %>%
 	  filter(!(is.infinite(ln_frac) | is.na(ln_frac))) %>%
 	  filter(!(is.infinite(ln_vol) | is.na(ln_vol))) %>%
 	  mutate(shape = ifelse(is.na(Bone_Lesions) & is.na(Pleural_Disease), 21, 24))
 	  	  
tmp.0 = tmp %>%
		mutate(cat = case_when(
			tot_volu >= 0 & tot_volu < quantile(tot_volu, 1/3) ~ 1,
			tot_volu >= quantile(tot_volu, 1/3) & tot_volu < quantile(tot_volu, 2/3) ~ 2,
			tot_volu >= quantile(tot_volu, 2/3) ~ 3,
		)) %>%
 		mutate(cat = as.factor(cat))
 
p = signif(jonckheere.test(x=tmp.0$ctdna_frac, g=as.numeric(tmp.0$cat), alternative = "increasing")$p.value, 3)
 
plot.0 = ggplot(tmp.0, aes(x = cat, y = ctdna_frac)) + 
		 geom_boxplot(alpha=1, outlier.size=NA, outlier.shape=NA, color="black", fill="white") +
		 geom_point(alpha = 1, size = 3.75, color = "black", fill = "salmon", aes(x = jitter(as.numeric(cat)), y = ctdna_frac), shape = tmp.0$shape, data=tmp.0, inherit.aes = FALSE) +
		 facet_wrap(~Tissue) +
		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
 		 labs(x="", y="ctDNA fraction\n") +
		 coord_cartesian(ylim = c(0,1)) +
		 annotate(geom="text", x=2, y=1, label=paste0("p = ", p), size=4.5)
		 
pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Breast_Box.pdf", width=5.5, height=5)
print(plot.0)
dev.off()
 
tmp = volumetric_data %>%
	  mutate(ln_frac = log(ctdna_frac)) %>%
	  mutate(ln_vol = log(tot_volu)) %>%
	  mutate(facet = "Lung") %>%
	  filter(grepl("VL", GRAIL_ID)) %>%
	  mutate(Tissue = "Lung") %>%
	  filter(!(is.infinite(ln_frac) | is.na(ln_frac))) %>%
	  filter(!(is.infinite(ln_vol) | is.na(ln_vol))) %>%
	  mutate(shape = ifelse(is.na(Bone_Lesions) & is.na(Pleural_Disease), 21, 24)) %>%
	  # Large left forearm bone and soft-tissue mass 10 by 7 by 17 cm not included in PET-scan and measurements
	  mutate(shape = ifelse(GRAIL_ID=="MSK-VL-0008", 24, shape)) 
	  
tmp.0 = tmp %>%
		mutate(cat = case_when(
			tot_volu >= 0 & tot_volu < quantile(tot_volu, 1/3) ~ 1,
			tot_volu >= quantile(tot_volu, 1/3) & tot_volu < quantile(tot_volu, 2/3) ~ 2,
			tot_volu >= quantile(tot_volu, 2/3) ~ 3,
		)) %>%
 		mutate(cat = as.factor(cat))

p = signif(jonckheere.test(x=tmp.0$ctdna_frac, g=as.numeric(tmp.0$cat), alternative = "two.sided")$p.value, 3)

plot.0 = ggplot(tmp.0, aes(x = cat, y = ctdna_frac)) + 
		 geom_boxplot(alpha=1, outlier.size=NA, outlier.shape=NA, color="black", fill="white") +
		 geom_point(alpha = 1, size = 3.75, color = "black", fill = "#FDAE61", aes(x = jitter(as.numeric(cat)), y = ctdna_frac), shape = tmp.0$shape, data=tmp.0, inherit.aes = FALSE) +
		 facet_wrap(~Tissue) +
		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
 		 labs(x="", y="ctDNA fraction\n") +
 		 coord_cartesian(ylim = c(0,1)) +
		 annotate(geom="text", x=2, y=1, label=paste0("p = ", p), size=4.5)
		 
pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Lung_Box.pdf", width=5.5, height=5)
print(plot.0)
dev.off()

tmp = volumetric_data %>%
	  mutate(tot_volu = ifelse(tot_volu==0, 0.1, tot_volu)) %>%
 	  mutate(ln_frac = log(ctdna_frac)) %>%
 	  mutate(ln_vol = log(tot_volu)) %>%
	  mutate(facet = "Prostate") %>%
	  filter(grepl("VP", GRAIL_ID)) %>%
	  mutate(Tissue = "Prostate") %>%
	  filter(!(is.infinite(ln_frac) | is.na(ln_frac))) %>%
	  filter(!(is.infinite(ln_vol) | is.na(ln_vol)))
tmp.0 = tmp %>%
		mutate(cat = case_when(
			tot_volu >= 0 & tot_volu <= quantile(tot_volu, 1/3) ~ 1,
			tot_volu > quantile(tot_volu, 1/3) & tot_volu <= quantile(tot_volu, 2/3) ~ 2,
			tot_volu > quantile(tot_volu, 2/3) ~ 3,
		)) %>%
 		mutate(cat = as.factor(cat))

p = signif(jonckheere.test(x=tmp.0$ctdna_frac, g=as.numeric(tmp.0$cat), alternative = "increasing")$p.value, 3)

plot.0 = ggplot(tmp.0, aes(x = cat, y = ctdna_frac)) + 
		 geom_boxplot(alpha=1, outlier.size=NA, outlier.shape=NA, color="black", fill="white") +
		 geom_point(alpha=1, size=3.75, shape=24, color="black", fill="#ABDDA4", aes(x = jitter(as.numeric(cat)), y = ctdna_frac), data=tmp.0) +
		 facet_wrap(~Tissue) +
		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
 		 labs(x="", y="ctDNA fraction\n") +
 		 coord_cartesian(ylim = c(0,1)) +
		 annotate(geom="text", x=2, y=1, label=paste0("p = ", p), size=4.5)

pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Prostate_Box.pdf", width=5.5, height=5)
print(plot.0)
dev.off()

#==================================================
# Box plot for cfDNA fraction by disease
# volume combined
#==================================================
pattern = c("VB", "VL", "VP")
tissue = c("Breast", "Lung", "Prostate")

tmp.1 = p.1 = NULL
for (i in 1:length(pattern)) {
	tmp = volumetric_data %>%
		  mutate(tot_volu = ifelse(tot_volu==0, 0.1, tot_volu)) %>%
	 	  mutate(ln_frac = log(ctdna_frac)) %>%
	 	  mutate(ln_vol = log(tot_volu)) %>%
		  filter(grepl(pattern[i], GRAIL_ID)) %>%
		  mutate(Tissue = tissue[i]) %>%
		  filter(!(is.infinite(ln_frac) | is.na(ln_frac))) %>%
		  filter(!(is.infinite(ln_vol) | is.na(ln_vol)))
		
	if (tissue[i]=="Breast") {
		tmp = tmp %>%
			  mutate(shape = ifelse(is.na(Bone_Lesions) & is.na(Pleural_Disease), 21, 24))
	} else if (tissue[i]=="Lung") {
		tmp = tmp %>%
			  mutate(shape = ifelse(is.na(Bone_Lesions) & is.na(Pleural_Disease), 21, 24)) %>%
	  		  # Large left forearm bone and soft-tissue mass 10 by 7 by 17 cm not included in PET-scan and measurements
	  		  mutate(shape = ifelse(GRAIL_ID=="MSK-VL-0008", 24, shape))
	} else if (tissue[i]=="Prostate") {
		tmp = tmp %>%
			  mutate(shape = 24)
	}
	

	tmp.0 = tmp %>%
			mutate(cat = case_when(
				tot_volu >= 0 & tot_volu <= quantile(tot_volu, 1/3) ~ 1,
				tot_volu > quantile(tot_volu, 1/3) & tot_volu <= quantile(tot_volu, 2/3) ~ 2,
				tot_volu > quantile(tot_volu, 2/3) ~ 3,
			)) %>%
 			mutate(cat = as.factor(cat))
 	p = signif(jonckheere.test(x=tmp.0$ctdna_frac, g=as.numeric(tmp.0$cat), alternative = "increasing")$p.value, 3)
 	tmp.1 = rbind(tmp.1, tmp.0)
 	p.1 = c(p.1, p)
}

plot.0 = ggplot(tmp.1, aes(x = Tissue, y = ctdna_frac, color = cat)) + 
		 geom_boxplot(alpha=1, outlier.size=NA, outlier.shape=NA, fill="white", width=.8) +
		 scale_colour_manual(values=c("1"="black", "2"="black", "3"="black")) +
		 geom_jitter(
		 	aes(x = Tissue, y = ctdna_frac, color = cat),
  			position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
  			alpha=1, size=3, shape=tmp.1$shape
  			) +
		 theme_classic(base_size=15) +
 		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10), legend.position="none") +
 		 labs(x="", y="ctDNA fraction\n") +
 		 coord_cartesian(ylim = c(0,1))

pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Combined_Box.pdf", width=8.5, height=5)
print(plot.0)
dev.off()

#==================================================
# Scatter plot of number of metastatic sites & ctDNA
# fraction
#==================================================

# Breast
tmp = volumetric_data %>%
 	  mutate(ln_frac = log(ctdna_frac)) %>%
 	  mutate(ln_vol = log(tot_volu)) %>%
 	  mutate(tissue = "Breast") %>%
 	  filter(grepl("VB", GRAIL_ID)) %>%
 	  filter(!(is.infinite(ln_frac) | is.na(ln_frac))) %>%
 	  filter(!(is.infinite(ln_vol) | is.na(ln_vol)))
 	  
z = cor.test(tmp$ln_frac, tmp$ln_vol, method="kendall")
tau = signif(z$estimate, 3)
p =  signif(z$p.value, 3)
 
plot.0 = ggplot(tmp, aes(x = ln_vol, y = ln_frac)) +
 		 geom_point(alpha = .8, size = 2.5, color = "black", fill = "#2B83BA", shape = 21) +
 		 geom_smooth(method = "lm", color = "goldenrod3", se=TRUE) +
 		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
 		 labs(x="\nDisease volume\n", y="\nctDNA fraction\n") +
 		 facet_wrap(~tissue) +
 		 scale_x_continuous(
 		 	breaks = function(x) { log(c(0, .1, 1, 10, 100, 500)) },
  			labels = function(x) { c(0, .1, 1, 10, 100, 500) }
 		 ) +
  		 scale_y_continuous(
  		 	breaks = function(x) { log(c(1, 0.2, 0.02, 0.002)) },
   			labels = function(x) { c("1", ".2", ".02", ".002") }
  		 ) +
  		 coord_cartesian(xlim=c(0,6.5)) +
  		 annotate(geom="text", x=log(200), y=log(0.0025), label=expression(tau~" = ")) +
  		 annotate(geom="text", x=log(200), y=log(0.0015), label=expression(P~" = "))

pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Breast_loglog.pdf", width=5.5, height=6)
print(plot.0)
dev.off()
 
z = cor.test(tmp$ctdna_frac, tmp$tot_volu, method="kendall")
tau = signif(z$estimate, 3)
p =  signif(z$p.value, 3)

plot.0 = ggplot(tmp, aes(x = tot_volu, y = ctdna_frac)) +
 		 geom_point(alpha = .8, size = 2.5, color = "black", fill = "#2B83BA", shape = 21) +
 		 geom_smooth(method = "lm", color = "goldenrod3", se=TRUE) +
 		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
 		 labs(x="\nDisease volume\n", y="\nctDNA fraction\n") +
 		 facet_wrap(~tissue) +
		 annotate(geom="text", x=450, y=0, label=expression(tau~" = ")) +
  		 annotate(geom="text", x=450, y=-.1, label=expression(P~" = "))
 		 
pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Breast_0log0log.pdf", width=5.5, height=6)
print(plot.0)
dev.off()
 
z = cor.test(tmp$ctdna_frac, tmp$ln_vol, method="kendall")
tau = signif(z$estimate, 3)
p =  signif(z$p.value, 3)
 
plot.0 = ggplot(tmp, aes(x = ln_vol, y = ctdna_frac)) +
 		 geom_point(alpha = .8, size = 2.5, color = "black", fill = "#2B83BA", shape = 21) +
 		 geom_smooth(method = "lm", color = "goldenrod3", se=TRUE) +
 		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
 		 labs(x="\nDisease volume\n", y="\nctDNA fraction\n") +
 		 facet_wrap(~tissue) +
 		 scale_x_continuous(
 		 	breaks = function(x) { log(c(0, .1, 1, 10, 100, 500)) },
  			labels = function(x) { c(0, .1, 1, 10, 100, 500) }
 		 ) +
 		 coord_cartesian(xlim=c(0,6.5), ylim=c(0,1)) +
 		 annotate(geom="text", x=log(200), y=0.05, label=expression(tau~" = ")) +
  		 annotate(geom="text", x=log(200), y=-.01, label=expression(P~" = "))

pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Breast_log0log.pdf", width=5.5, height=6)
print(plot.0)
dev.off()
 
z = cor.test(tmp$ln_frac, tmp$tot_volu, method="kendall")
tau = signif(z$estimate, 3)
p =  signif(z$p.value, 3) 

plot.0 = ggplot(tmp, aes(x = tot_volu, y = ln_frac)) +
 		 geom_point(alpha = .8, size = 2.5, color = "black", fill = "#2B83BA", shape = 21) +
 		 geom_smooth(method = "lm", color = "goldenrod3", se=TRUE) +
 		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
 		 labs(x="\nDisease volume\n", y="\nctDNA fraction\n") +
 		 facet_wrap(~tissue) +
  		 scale_y_continuous(
  		 	breaks = function(x) { log(c(1, 0.2, 0.02, 0.002)) },
   			labels = function(x) { c("1", ".2", ".02", ".002") }
  		 ) +
  		 coord_cartesian(ylim=c(-6.5,1)) +
  		 annotate(geom="text", x=500, y=-6.1, label=expression(tau~" = ")) +
  		 annotate(geom="text", x=500, y=-6.5, label=expression(P~" = "))

pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Breast_0loglog.pdf", width=5.5, height=6)
print(plot.0)
dev.off()
 
# Lung
tmp = volumetric_data %>%
 	  mutate(ln_frac = log(ctdna_frac)) %>%
 	  mutate(ln_vol = log(tot_volu)) %>%
 	  mutate(tissue = "Lung") %>%
 	  filter(grepl("VL", GRAIL_ID)) %>%
 	  filter(!(is.infinite(ln_frac) | is.na(ln_frac))) %>%
 	  filter(!(is.infinite(ln_vol) | is.na(ln_vol)))
 	  
z = cor.test(tmp$ln_frac, tmp$ln_vol, method="kendall")
tau = signif(z$estimate, 3)
p =  signif(z$p.value, 3) 


plot.0 = ggplot(tmp, aes(x = ln_vol, y = ln_frac)) +
 		 geom_point(alpha = .8, size = 2.5, color = "black", fill = "#2B83BA", shape = 21) +
 		 geom_smooth(method = "lm", color = "goldenrod3", se=TRUE) +
 		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
 		 labs(x="\nDisease volume\n", y="\nctDNA fraction\n") +
 		 facet_wrap(~tissue) +
 		 scale_x_continuous(
 		 	breaks = function(x) { log(c(0, .1, 1, 10, 100, 500)) },
  			labels = function(x) { c(0, .1, 1, 10, 100, 500) }
 		 ) +
  		 scale_y_continuous(
  		 	breaks = function(x) { log(c(1, 0.2, 0.02, 0.002, 0.0003)) },
   			labels = function(x) { c("1", ".2", ".02", ".002", ".0003") }
  		 ) +
  		 coord_cartesian(xlim=c(0,6.5), ylim = c(-8.5, 1.1)) +
  		 annotate(geom="text", x=log(200), y=log(.0003), label=expression(tau~" = ")) +
  		 annotate(geom="text", x=log(200), y=log(.0002), label=expression(P~" = "))

pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Lung_loglog.pdf", width=5.5, height=6)
print(plot.0)
dev.off()
 
z = cor.test(tmp$ctdna_frac, tmp$tot_volu, method="kendall")
tau = signif(z$estimate, 3)
p =  signif(z$p.value, 3) 
 
plot.0 = ggplot(tmp, aes(x = tot_volu, y = ctdna_frac)) +
		 geom_point(alpha = .8, size = 2.5, color = "black", fill = "#2B83BA", shape = 21) +
		 geom_smooth(method = "lm", color = "goldenrod3", se=TRUE) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\nDisease volume\n", y="\nctDNA fraction\n") +
		 facet_wrap(~tissue) +
		 annotate(geom="text", x=300, y=0, label=expression(tau~" = ")) +
  		 annotate(geom="text", x=300, y=-.03, label=expression(P~" = "))
  		 
pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Lung_0log0log.pdf", width=5.5, height=6)
print(plot.0)
dev.off()
 
z = cor.test(tmp$ctdna_frac, tmp$ln_vol, method="kendall")
tau = signif(z$estimate, 3)
p =  signif(z$p.value, 3) 

plot.0 = ggplot(tmp, aes(x = ln_vol, y = ctdna_frac)) +
		 geom_point(alpha = .8, size = 2.5, color = "black", fill = "#2B83BA", shape = 21) +
		 geom_smooth(method = "lm", color = "goldenrod3", se=TRUE) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\nDisease volume\n", y="\nctDNA fraction\n") +
		 facet_wrap(~tissue) +
		 scale_x_continuous(
		 	breaks = function(x) { log(c(0, .1, 1, 10, 100, 500)) },
 			labels = function(x) { c(0, .1, 1, 10, 100, 500) }
		 ) +
		 coord_cartesian(xlim=c(0,6.5)) +
		 annotate(geom="text", x=log(1.5), y=0.5, label=expression(tau~" = ")) +
  		 annotate(geom="text", x=log(1.5), y=0.47, label=expression(P~" = "))
 		 
pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Lung_log0log.pdf", width=5.5, height=6)
print(plot.0)
dev.off()
 
z = cor.test(tmp$ln_frac, tmp$tot_volu, method="kendall")
tau = signif(z$estimate, 3)
p =  signif(z$p.value, 3) 

plot.0 = ggplot(tmp, aes(x = tot_volu, y = ln_frac)) +
		 geom_point(alpha = .8, size = 2.5, color = "black", fill = "#2B83BA", shape = 21) +
		 geom_smooth(method = "lm", color = "goldenrod3", se=TRUE) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\nDisease volume\n", y="\nctDNA fraction\n") +
		 facet_wrap(~tissue) +
 		 scale_y_continuous(
 		 	breaks = function(x) { log(c(1, 0.2, 0.02, 0.002, 0.0003)) },
  			labels = function(x) { c("1", ".2", ".02", ".002", ".0003") }
  		 ) +
  		 coord_cartesian(ylim = c(-8.5, 1.1)) +
  		 annotate(geom="text", x=300, y=log(0.0003), label=expression(tau~" = ")) +
  		 annotate(geom="text", x=300, y=log(0.0002), label=expression(P~" = "))
 		 
pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Lung_0loglog.pdf", width=5.5, height=6)
print(plot.0)
dev.off()
 
# Prostate
tmp = volumetric_data %>%
	  mutate(tot_volu = ifelse(tot_volu==0, 0.1, tot_volu)) %>%
	  mutate(ln_frac = log(ctdna_frac)) %>%
	  mutate(ln_vol = log(tot_volu)) %>%
	  mutate(tissue = "Prostate") %>%
	  filter(grepl("VP", GRAIL_ID)) %>%
	  filter(!(is.infinite(ln_frac) | is.na(ln_frac))) %>%
	  filter(!(is.infinite(ln_vol) | is.na(ln_vol)))
 
z = cor.test(tmp$ln_frac, tmp$ln_vol, method="kendall")
tau = signif(z$estimate, 3)
p =  signif(z$p.value, 3) 

plot.0 = ggplot(tmp, aes(x = ln_vol, y = ln_frac)) +
		 geom_point(alpha = .8, size = 2.5, color = "black", fill = "#2B83BA", shape = 21) +
		 geom_smooth(method = "lm", color = "goldenrod3", se=TRUE) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\nBone scan index\n", y="\nctDNA fraction\n") +
		 facet_wrap(~tissue) +
		 scale_x_continuous(
		 	breaks = function(x) { log(c(0, .1, 1, 10, 100, 500)) },
 			labels = function(x) { c(0, .1, 1, 10, 100, 500) }
		 ) +
 		 scale_y_continuous(
 		 	breaks = function(x) { log(c(1, 0.2, 0.02, 0.002, 0.0003)) },
  			labels = function(x) { c("1", ".2", ".02", ".002", ".0003") }
 		 ) +
 		 coord_cartesian(xlim=c(-2.5,3), ylim = c(-8, 1.1)) +
 		 annotate(geom="text", x=-2, y=log(1.8), label=expression(tau~" = ")) +
  		 annotate(geom="text", x=-2, y=log(1), label=expression(P~" = "))
  		 
pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Prostate_loglog.pdf", width=5.5, height=6)
print(plot.0)
dev.off()

z = cor.test(tmp$ctdna_frac, tmp$tot_volu, method="kendall")
tau = signif(z$estimate, 3)
p =  signif(z$p.value, 3) 

plot.0 = ggplot(tmp, aes(x = tot_volu, y = ctdna_frac)) +
		 geom_point(alpha = .8, size = 2.5, color = "black", fill = "#2B83BA", shape = 21) +
		 geom_smooth(method = "lm", color = "goldenrod3", se=TRUE) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
 		 labs(x="\nBone scan index\n", y="\nctDNA fraction\n") +
		 facet_wrap(~tissue) +
		 annotate(geom="text", x=1, y=.7, label=expression(tau~" = ")) +
  		 annotate(geom="text", x=1, y=.65, label=expression(P~" = "))

pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Prostate_0log0log.pdf", width=5.5, height=6)
print(plot.0)
dev.off()
 
z = cor.test(tmp$ctdna_frac, tmp$ln_vol, method="kendall")
tau = signif(z$estimate, 3)
p =  signif(z$p.value, 3) 
 
plot.0 = ggplot(tmp, aes(x = ln_vol, y = ctdna_frac)) +
		 geom_point(alpha = .8, size = 2.5, color = "black", fill = "#2B83BA", shape = 21) +
		 geom_smooth(method = "lm", color = "goldenrod3", se=TRUE) +
		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
 		 labs(x="\nBone scan index\n", y="\nctDNA fraction\n") +
 		 facet_wrap(~tissue) +
 		 scale_x_continuous(
 		 	breaks = function(x) { log(c(0, .1, 1, 10, 100, 500)) },
  			labels = function(x) { c(0, .1, 1, 10, 100, 500) }
 		 ) +
 		 coord_cartesian(xlim=c(-2.5,3)) +
 		 annotate(geom="text", x=-2.4, y=.7, label=expression(tau~" = ")) +
  		 annotate(geom="text", x=-2.4, y=.65, label=expression(P~" = "))
 		 
pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Prostate_log0log.pdf", width=5.5, height=6)
print(plot.0)
dev.off()

z = cor.test(tmp$ln_frac, tmp$tot_volu, method="kendall")
tau = signif(z$estimate, 3)
p =  signif(z$p.value, 3) 

plot.0 = ggplot(tmp, aes(x = tot_volu, y = ln_frac)) +
 		 geom_point(alpha = .8, size = 2.5, color = "black", fill = "#2B83BA", shape = 21) +
 		 geom_smooth(method = "lm", color = "goldenrod3", se=TRUE) +
 		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
 		 labs(x="\nBone scan index\n", y="\nctDNA fraction\n") +
 		 facet_wrap(~tissue) +
  		 scale_y_continuous(
  		 	breaks = function(x) { log(c(1, 0.2, 0.02, 0.002, 0.0003)) },
   			labels = function(x) { c("1", ".2", ".02", ".002", ".0003") }
  		 ) +
  		 coord_cartesian(ylim = c(-8, 1.1)) +
  		 annotate(geom="text", x=1, y=1, label=expression(tau~" = ")) +
  		 annotate(geom="text", x=1, y=.5, label=expression(P~" = "))
 		 
pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Prostate_0loglog.pdf", width=5.5, height=6)
print(plot.0)
dev.off()
 
#==================================================
# Disease volume versus number of metastatic sites
#==================================================
ctdna_fraction = read_csv(file=url_ctdna_frac, col_types = cols(.default = col_character()))  %>%
 				 type_convert() %>%
 				 rename(GRAIL_ID = ID)

vol_breast = read.xls(xls=url_volumetric_data, sheet=1, stringsAsFactors=FALSE) %>%
			 type_convert() %>%
			 mutate(Date_Tissue = as.character(Date_Tissue)) %>%
			 mutate(Date_Most_Recent_Scan = as.character(Date_Most_Recent_Scan)) %>%
			 mutate(Date_of_Last_Treatment = as.character(Date_of_Last_Treatment))

vol_lung = read.xls(xls=url_volumetric_data, sheet=2, stringsAsFactors=FALSE) %>%
		   type_convert() %>%
		   mutate(Date_Tissue = as.character(Date_Tissue)) %>%
		   mutate(Date_Most_Recent_Scan = as.character(Date_Most_Recent_Scan)) %>%
		   mutate(Date_of_Last_Treatment = as.character(Date_of_Last_Treatment)) %>%
		   mutate(Date_Dx = as.character(Date_Dx)) %>%
		   mutate(Date_Metastatic = as.character(Date_Metastatic)) %>%
		   mutate(Date_Tissue = as.character(Date_Tissue))
		   
vol_prostate = read.xls(xls=url_volumetric_data, sheet=3, stringsAsFactors=FALSE) %>%
		       type_convert() %>%
		       mutate(BSI_Value = ifelse(BSI_Value==">13", 14, BSI_Value)) %>%
		       mutate(BSI_Value = ifelse(BSI_Value=="n/a", NA, BSI_Value)) %>%
		       mutate(BSI_Value = as.numeric(BSI_Value)) %>%
		       rename(BSI_Volume = BSI_Value) %>%
		       mutate(Date_Dx = as.character(Date_Dx)) %>%
		       mutate(Date_Metastatic = as.character(Date_Metastatic)) %>%
		       mutate(Date_Tissue = as.character(Date_Tissue))

volumetric_data = bind_rows(vol_breast, vol_lung) %>%
				  bind_rows(vol_prostate) %>%
				  left_join(ctdna_fraction, by="GRAIL_ID")
				  
volumetric_data = volumetric_data %>%
				  mutate(tot_volu = apply( volumetric_data[,which(grepl('volume', colnames(volumetric_data), ignore.case = T))] , 1 , function(x) sum(as.numeric(x), na.rm=T) )) %>%
				  mutate(tot_mets = apply( volumetric_data[,which(!grepl('volume', colnames(volumetric_data), ignore.case = T) & grepl('Lesions|Lymph_Nodes|Pleural_Disease', colnames(volumetric_data), ignore.case = T))] , 1 , function(x) sum(as.numeric(x), na.rm=T))) %>%
				  filter(!(GRAIL_ID %in% c('MSK-VL-0057','MSK-VL-0003','MSK-VB-0046')))
				  
clinical = read_tsv(file=clinical_file_updated, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
cfdna_frac = read.csv(file=url_ctdna_frac, header=TRUE, sep=",", stringsAsFactors=FALSE) %>%
			 filter(!is.na(ctdna_frac)) %>%
			 mutate(index = order_samples(ID)) %>%
 			 arrange(desc(ctdna_frac)) %>%
			 arrange(index)
number_metastatic_sites = unlist(lapply(strsplit(clinical[,"metastatic_sites",drop=TRUE], split=",", fixed=TRUE), function(x) {length(x)}))
number_metastatic_sites[is.na(clinical[,"metastatic_sites",drop=TRUE])] = 0
names(number_metastatic_sites) = clinical[,"patient_id",drop=TRUE]
number_metastatic_sites = number_metastatic_sites[cfdna_frac[,1]]
cfdna_frac = cbind(cfdna_frac, number_metastatic_sites) %>%
			 mutate(Tissue = "") %>%
			 mutate(Tissue = ifelse(grepl("VB", ID), "Breast", Tissue)) %>%
			 mutate(Tissue = ifelse(grepl("VL", ID), "Lung", Tissue)) %>%
			 mutate(Tissue = ifelse(grepl("VP", ID), "Prostate", Tissue)) %>%
			 mutate(cat = case_when(
			 		number_metastatic_sites==1 | number_metastatic_sites==2 ~ 1,
			 		number_metastatic_sites==3 ~ 2,
			 		number_metastatic_sites>=4 ~ 3)) %>%
			 mutate(GRAIL_ID = ID)
			 
tmp = left_join(volumetric_data, cfdna_frac, by="GRAIL_ID") %>%
	  dplyr::select(patient_id = GRAIL_ID, tot_volu, tot_mets, tot_mets_2 = number_metastatic_sites, tissue = Tissue, ctdna_frac = ctdna_frac.x)

# Breast
tmp.0 = tmp %>%
 	    filter(tissue == "Breast") %>%
 	    mutate(ln_frac = log(ctdna_frac)) %>%
 	  	mutate(ln_vol = log(tot_volu)) %>%
 	    filter(!(is.infinite(ln_frac) | is.na(ln_frac))) %>%
 	    filter(!(is.infinite(ln_vol) | is.na(ln_vol)))
 	  
z = cor.test(tmp.0$tot_volu, tmp.0$tot_mets, method="spearman", alternative="greater", exact=FALSE)
rho = signif(z$estimate, 3)
p =  signif(z$p.value, 3)
 
plot.0 = ggplot(tmp.0, aes(x = tot_volu, y = tot_mets)) +
 		 geom_point(alpha = 1, size = 3.5, color = "black", fill = "salmon", shape = 21) +
 		 geom_smooth(method = "lm", color = "goldenrod3", se=TRUE) +
 		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
 		 labs(x="\nDisease volume\n", y="\nNumber of metastases\n") +
 		 facet_wrap(~tissue) +
  		 coord_cartesian(xlim=c(1, 700), ylim=c(0,70)) +
  		 scale_x_log10(
 		 	breaks = c(1, 10, 100, 700),
  			labels = c(1, 10, 100, 700)
 		 ) +
 		 scale_y_continuous(
 		 	breaks = c(0, 10, 20, 30, 40, 50, 60, 70),
  			labels = c(0, 10, 20, 30, 40, 50, 60, "70.0")
 		 ) +
 		 annotate(geom="text", x=log(10), y=65, label=expression(rho~" = ")) +
 		 annotate(geom="text", x=log(15), y=65, label=rho) +
  		 annotate(geom="text", x=log(10), y=60, label=expression(p~" = ")) +
  		 annotate(geom="text", x=log(15), y=60, label=p)

pdf(file="../res/rebuttal/N_mets_vs_Tv_Breast.pdf", width=5.5, height=5.3)
print(plot.0)
dev.off()

# Lung
tmp.0 = tmp %>%
 	    filter(tissue == "Lung") %>%
 	    mutate(ln_frac = log(ctdna_frac)) %>%
 	  	mutate(ln_vol = log(tot_volu)) %>%
 	    filter(!(is.infinite(ln_frac) | is.na(ln_frac))) %>%
 	    filter(!(is.infinite(ln_vol) | is.na(ln_vol)))
 	  
z = cor.test(tmp.0$tot_volu, tmp.0$tot_mets, method="spearman", alternative="greater", exact=FALSE)
rho = signif(z$estimate, 3)
p =  signif(z$p.value, 3)
 
plot.0 = ggplot(tmp.0, aes(x = tot_volu, y = tot_mets)) +
 		 geom_point(alpha = 1, size = 3.5, color = "black", fill = "#FDAE61", shape = 21) +
 		 geom_smooth(method = "lm", color = "goldenrod3", se=TRUE) +
 		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
 		 labs(x="\nDisease volume\n", y="\nNumber of metastases\n") +
 		 facet_wrap(~tissue) +
  		 coord_cartesian(xlim=c(1, 700), ylim=c(0,25)) +
  		 scale_x_log10(
 		 	breaks = c(1, 10, 100, 700),
  			labels = c(1, 10, 100, 700)
 		 ) +
 		 scale_y_continuous(
 		 	breaks = c(0, 5, 10, 15, 20, 25),
  			labels = c(0, 5, 10, 15, 20, "25.0")
 		 ) +
 		 annotate(geom="text", x=log(10), y=25, label=expression(rho~" = ")) +
 		 annotate(geom="text", x=log(15), y=25, label=rho) +
  		 annotate(geom="text", x=log(10), y=20, label=expression(p~" = ")) +
  		 annotate(geom="text", x=log(15), y=20, label=p)

pdf(file="../res/rebuttal/N_mets_vs_Tv_Lung.pdf", width=5.5, height=5.3)
print(plot.0)
dev.off()

# Prostate
tmp.0 = tmp %>%
		mutate(tot_volu = ifelse(tot_volu==0, 0.1, tot_volu)) %>%
 	    filter(tissue == "Prostate") %>%
 	    mutate(ln_frac = log(ctdna_frac)) %>%
 	  	mutate(ln_vol = log(tot_volu)) %>%
	    filter(!(is.infinite(ln_frac) | is.na(ln_frac))) %>%
	    filter(!(is.infinite(ln_vol) | is.na(ln_vol)))
 	  
z = cor.test(tmp.0$tot_volu, tmp.0$tot_mets_2, method="spearman", alternative="two.sided", exact=FALSE)
rho = signif(z$estimate, 3)
p =  signif(z$p.value, 3)
 
plot.0 = ggplot(tmp.0, aes(x = tot_volu, y = tot_mets_2)) +
 		 geom_point(alpha = 1, size = 3.5, color = "black", fill = "#ABDDA4", shape = 21) +
 		 geom_smooth(method = "lm", color = "goldenrod3", se=TRUE) +
 		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
 		 labs(x="\nBone scan index\n", y="\nNumber of metastases\n") +
 		 facet_wrap(~tissue) +
  		 coord_cartesian(xlim=c(0, 15), ylim=c(0,10)) +
  		 scale_y_continuous(
 		 	breaks = c(0, 2, 4, 6, 8, 10),
  			labels = c(0, 2, 4, 6, 8, "10.0")
 		 ) +
 		 annotate(geom="text", x=0, y=10, label=expression(rho~" = ")) +
 		 annotate(geom="text", x=1, y=10, label=rho) +
  		 annotate(geom="text", x=0, y=9, label=expression(p~" = ")) +
  		 annotate(geom="text", x=1, y=9, label=p)

pdf(file="../res/rebuttal/N_mets_vs_Tv_Prostate.pdf", width=5.5, height=5.3)
print(plot.0)
dev.off()

#==================================================
# Breast and Lung combined
#==================================================
tmp = volumetric_data %>%
	  mutate(ln_frac = log(ctdna_frac)) %>%
	  mutate(ln_vol = log(tot_volu)) %>%
	  mutate(facet = "Breast & Lung") %>%
	  filter(grepl("VL", GRAIL_ID) | grepl("VB", GRAIL_ID)) %>%
	  mutate(Tissue = "Breast") %>%
	  mutate(Tissue = ifelse(grepl("VL", GRAIL_ID), "Lung", Tissue)) %>%
	  filter(!(is.infinite(ln_frac) | is.na(ln_frac))) %>%
	  filter(!(is.infinite(ln_vol) | is.na(ln_vol)))
	  
cols = c("Breast"="#D7191C", "Lung"="#2B83BA")
	  
zB = cor.test(tmp$ln_frac[tmp$Tissue=="Breast"], tmp$tot_volu[tmp$Tissue=="Breast"], method="kendall")
tauB = signif(zB$estimate, 3)
pB =  signif(zB$p.value, 3)
zL = cor.test(tmp$ln_frac[tmp$Tissue=="Lung"], tmp$tot_volu[tmp$Tissue=="Lung"], method="kendall")
tauL = signif(zL$estimate, 3)
pL =  signif(zL$p.value, 3) 

plot.0 = ggplot(tmp, aes(x = ln_vol, y = ln_frac, fill = Tissue)) +
 		 geom_point(alpha = .8, size = 2.5, color = "black", shape = 21) +
 		 geom_smooth(method = "lm", se=TRUE, data=tmp, aes(x = ln_vol, y = ln_frac, color = Tissue)) +
 		 scale_fill_manual(values = cols) +
 		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.9), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
 		 labs(x="\nDisease volume\n", y="\nctDNA fraction\n") +
 		 facet_wrap(~facet) +
 		 scale_x_continuous(
 		 	breaks = function(x) { log(c(0, .1, 1, 10, 100, 500)) },
  			labels = function(x) { c(0, .1, 1, 10, 100, 500) }
 		 ) +
  		 scale_y_continuous(
  		 	breaks = function(x) { log(c(1, 0.2, 0.02, 0.002, 0.0003)) },
   			labels = function(x) { c("1", ".2", ".02", ".002", ".0003") }
  		 ) +
  		 coord_cartesian(xlim=c(0,6.5), ylim = c(-8.5, 1.1)) +
  		 annotate(geom="text", x=log(100), y=log(0.002), label=expression(tau~" = ")) +
  		 annotate(geom="text", x=log(100), y=log(0.002)-.5, label=expression(P~" = ")) +
  		 annotate(geom="text", x=log(100), y=log(0.002)-1, label=expression(tau~" = ")) +
  		 annotate(geom="text", x=log(100), y=log(0.002)-1.5, label=expression(P~" = "))
 		 
pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Breast_Lung_loglog.pdf", width=5.5, height=6)
print(plot.0)
dev.off()

zB = cor.test(tmp$ctdna_frac[tmp$Tissue=="Breast"], tmp$tot_volu[tmp$Tissue=="Breast"], method="kendall")
tauB = signif(zB$estimate, 3)
pB =  signif(zB$p.value, 3)
zL = cor.test(tmp$ctdna_frac[tmp$Tissue=="Lung"], tmp$tot_volu[tmp$Tissue=="Lung"], method="kendall")
tauL = signif(zL$estimate, 3)
pL =  signif(zL$p.value, 3) 
 
plot.0 = ggplot(tmp, aes(x = tot_volu, y = ctdna_frac, fill = Tissue)) +
		 geom_point(alpha = .8, size = 2.5, color = "black", shape = 21) +
		 geom_smooth(method = "lm", se=FALSE, data=tmp, aes(x = tot_volu, y = ctdna_frac, color = Tissue)) +
		 scale_fill_manual(values = cols) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.9), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\nDisease volume\n", y="\nctDNA fraction\n") +
		 facet_wrap(~facet) +
		 annotate(geom="text", x=500, y=.25, label=expression(tau~" = ")) +
  		 annotate(geom="text", x=500, y=.20, label=expression(P~" = ")) +
  		 annotate(geom="text", x=500, y=.15, label=expression(tau~" = ")) +
  		 annotate(geom="text", x=500, y=.10, label=expression(P~" = "))
 		 
pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Breast_Lung_0log0log.pdf", width=5.5, height=6)
print(plot.0)
dev.off()
 
zB = cor.test(tmp$ctdna_frac[tmp$Tissue=="Breast"], tmp$ln_vol[tmp$Tissue=="Breast"], method="kendall")
tauB = signif(zB$estimate, 3)
pB =  signif(zB$p.value, 3)
zL = cor.test(tmp$ctdna_frac[tmp$Tissue=="Lung"], tmp$ln_vol[tmp$Tissue=="Lung"], method="kendall")
tauL = signif(zL$estimate, 3)
pL =  signif(zL$p.value, 3) 

plot.0 = ggplot(tmp, aes(x = ln_vol, y = ctdna_frac, fill = Tissue)) +
 		 geom_point(alpha = .8, size = 2.5, color = "black", shape = 21) +
 		 geom_smooth(method = "lm", se=TRUE, data=tmp, aes(x = ln_vol, y = ctdna_frac, color = Tissue)) +
 		 scale_fill_manual(values = cols) +
 		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.9), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
 		 labs(x="\nDisease volume\n", y="\nctDNA fraction\n") +
 		 facet_wrap(~facet) +
 		 scale_x_continuous(
 		 	breaks = function(x) { log(c(0, .1, 1, 10, 100, 500)) },
  			labels = function(x) { c(0, .1, 1, 10, 100, 500) }
 		 ) +
  		 coord_cartesian(xlim=c(0,6.5), ylim=c(-.1, 1)) +
  		 annotate(geom="text", x=log(400), y=.25, label=expression(tau~" = ")) +
  		 annotate(geom="text", x=log(400), y=.20, label=expression(P~" = ")) +
  		 annotate(geom="text", x=log(400), y=.15, label=expression(tau~" = ")) +
  		 annotate(geom="text", x=log(400), y=.10, label=expression(P~" = "))
 		 
pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Breast_Lung_log0log.pdf", width=5.5, height=6)
print(plot.0)
dev.off()
 
zB = cor.test(tmp$ln_frac[tmp$Tissue=="Breast"], tmp$tot_volu[tmp$Tissue=="Breast"], method="kendall")
tauB = signif(zB$estimate, 3)
pB =  signif(zB$p.value, 3)
zL = cor.test(tmp$ln_frac[tmp$Tissue=="Lung"], tmp$tot_volu[tmp$Tissue=="Lung"], method="kendall")
tauL = signif(zL$estimate, 3)
pL =  signif(zL$p.value, 3)

plot.0 = ggplot(tmp, aes(x = tot_volu, y = ln_frac, fill = Tissue)) +
 		 geom_point(alpha = .8, size = 2.5, color = "black", shape = 21) +
 		 geom_smooth(method = "lm", se=TRUE, data=tmp, aes(x = tot_volu, y = ln_frac, color = Tissue)) +
 		 scale_fill_manual(values = cols) +
 		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.9), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
 		 labs(x="\nDisease volume\n", y="\nctDNA fraction\n") +
 		 facet_wrap(~facet) +
  		 scale_y_continuous(
  		 	breaks = function(x) { log(c(1, 0.2, 0.02, 0.002, 0.0003)) },
   			labels = function(x) { c("1", ".2", ".02", ".002", ".0003") }
  		 ) +
  		 coord_cartesian(ylim = c(-8.5, 1.1)) +
  		 annotate(geom="text", x=500, y=log(.002), label=expression(tau~" = ")) +
  		 annotate(geom="text", x=500, y=log(.002)-.5, label=expression(P~" = ")) +
  		 annotate(geom="text", x=500, y=log(.002)-1, label=expression(tau~" = ")) +
  		 annotate(geom="text", x=500, y=log(.002)-1.5, label=expression(P~" = "))
 		 
pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Breast_Lung_0loglog.pdf", width=5.5, height=6)
print(plot.0)
dev.off()
