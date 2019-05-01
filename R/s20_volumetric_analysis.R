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
# Scatter plot of number of metastatic sites & ctDNA
# fraction
#==================================================
ctdna_fraction = read_csv(file=url_ctdna, col_types = cols(.default = col_character()))  %>%
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
				  mutate(tot_volu = apply( volumetric_data[,which(grepl('volume', colnames(volumetric_data), ignore.case = T))] , 1 , function(x) sum(x, na.rm=T) )) %>%
				  mutate(tot_mets = apply( volumetric_data[,which(!grepl('volume', colnames(volumetric_data), ignore.case = T) & grepl('Lesions|Lymph_Nodes|Pleural_Disease', colnames(volumetric_data), ignore.case = T))] , 1 , function(x) sum(x, na.rm=T))) %>%
				  filter(!(GRAIL_ID %in% c('MSK-VL-0057','MSK-VL-0003','MSK-VB-0046')))
 
#==================================================
# Breast
#==================================================
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
 
#==================================================
# Lung
#==================================================
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
 
#==================================================
# Prostate
#==================================================
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
# Box plot for Breast
#==================================================
tmp = volumetric_data %>%
 	  mutate(ln_frac = log(ctdna_frac)) %>%
 	  mutate(ln_vol = log(tot_volu)) %>%
 	  mutate(facet = "Breast") %>%
 	  filter(grepl("VB", GRAIL_ID)) %>%
 	  mutate(Tissue = "Breast") %>%
 	  filter(!(is.infinite(ln_frac) | is.na(ln_frac))) %>%
 	  filter(!(is.infinite(ln_vol) | is.na(ln_vol)))
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
		 geom_point(alpha=1, size=3.75, shape=21, color="black", fill="salmon", aes(x = jitter(as.numeric(cat)), y = ctdna_frac), data=tmp.0) +
		 facet_wrap(~Tissue) +
		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
 		 labs(x="", y="ctDNA fraction\n") +
		 coord_cartesian(ylim = c(0,1)) +
		 annotate(geom="text", x=2, y=1, label=paste0("p = ", p))
		 
pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Breast_Box.pdf", width=5.5, height=6)
print(plot.0)
dev.off()
 
#==================================================
# Box plot for Lung
#==================================================
tmp = volumetric_data %>%
	  mutate(ln_frac = log(ctdna_frac)) %>%
	  mutate(ln_vol = log(tot_volu)) %>%
	  mutate(facet = "Lung") %>%
	  filter(grepl("VL", GRAIL_ID)) %>%
	  mutate(Tissue = "Lung") %>%
	  filter(!(is.infinite(ln_frac) | is.na(ln_frac))) %>%
	  filter(!(is.infinite(ln_vol) | is.na(ln_vol)))
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
		 geom_point(alpha=1, size=3.75, shape=21, color="black", fill="#FDAE61", aes(x = jitter(as.numeric(cat)), y = ctdna_frac), data=tmp.0) +
		 facet_wrap(~Tissue) +
		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
 		 labs(x="", y="ctDNA fraction\n") +
 		 coord_cartesian(ylim = c(0,1)) +
		 annotate(geom="text", x=2, y=1, label=paste0("p = ", p))
		 
pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Lung_Box.pdf", width=5.5, height=6)
print(plot.0)
dev.off()

#==================================================
# Box plot for Prostate
#==================================================
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
		 geom_point(alpha=1, size=3.75, shape=21, color="black", fill="#ABDDA4", aes(x = jitter(as.numeric(cat)), y = ctdna_frac), data=tmp.0) +
		 facet_wrap(~Tissue) +
		 theme_bw(base_size=15) +
 		 theme(axis.text.y = element_text(size=13), axis.text.x = element_text(size=10)) +
 		 labs(x="", y="ctDNA fraction\n") +
 		 coord_cartesian(ylim = c(0,1)) +
		 annotate(geom="text", x=2, y=1, label=paste0("p = ", p))

pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Prostate_Box.pdf", width=5.5, height=6)
print(plot.0)
dev.off()

# #==================================================
# # Breast and lung combined
# #==================================================
# tmp = volumetric_data %>%
# 	  mutate(ln_frac = log(ctdna_frac)) %>%
# 	  mutate(ln_vol = log(tot_volu)) %>%
# 	  mutate(facet = "Breast & Lung") %>%
# 	  filter(grepl("VL", GRAIL_ID) | grepl("VB", GRAIL_ID)) %>%
# 	  mutate(Tissue = "Breast") %>%
# 	  filter(!(is.infinite(ln_frac) | is.na(ln_frac))) %>%
# 	  filter(!(is.infinite(ln_vol) | is.na(ln_vol)))
# tmp[grepl("VL", tmp[,"GRAIL_ID"]),"Tissue"] = "Lung"
# 
# cols = c("Breast"="#D7191C", "Lung"="#2B83BA")
# 
# r2b = signif(summary(lm(ln_frac ~ ln_vol, data=tmp %>% filter(Tissue=="Breast")))$r.squared, 3)
# pb = signif(summary(lm(ln_frac ~ ln_vol, data=tmp %>% filter(Tissue=="Breast")))$coefficients[2,4], 3)
# r2l = signif(summary(lm(ln_frac ~ ln_vol, data=tmp %>% filter(Tissue=="Lung")))$r.squared, 3)
# pl = signif(summary(lm(ln_frac ~ ln_vol, data=tmp %>% filter(Tissue=="Lung")))$coefficients[2,4], 3)
# 
# r2b = signif(cor(tmp$ln_frac[tmp$Tissue=="Breast"], tmp$ln_vol[tmp$Tissue=="Breast"], method="kendall"), 3)
# pb = signif(cor.test(tmp$ln_frac[tmp$Tissue=="Breast"], tmp$ln_vol[tmp$Tissue=="Breast"], method="kendall")$p.value, 3)
# r2l = signif(cor(tmp$ln_frac[tmp$Tissue=="Lung"], tmp$ln_vol[tmp$Tissue=="Lung"], method="kendall"), 3)
# pl = signif(cor.test(tmp$ln_frac[tmp$Tissue=="Lung"], tmp$ln_vol[tmp$Tissue=="Lung"], method="kendall")$p.value, 3)
# 	  
# plot.0 = ggplot(tmp, aes(x = ln_vol, y = ln_frac, fill = Tissue)) +
# 		 geom_point(alpha = .8, size = 2.5, color = "black", shape = 21) +
# 		 geom_smooth(method = "lm", se=FALSE, data=tmp, aes(x = ln_vol, y = ln_frac, color = Tissue)) +
# 		 scale_fill_manual(values = cols) +
# 		 theme_bw(base_size=15) +
# 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.9), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# 		 labs(x="\nDisease volume\n", y="\nctDNA fraction\n") +
# 		 facet_wrap(~facet) +
# 		 scale_x_continuous(
# 		 	breaks = function(x) { log(c(0, .1, 1, 10, 100, 500)) },
#  			labels = function(x) { c(0, .1, 1, 10, 100, 500) }
# 		 ) +
#  		 scale_y_continuous(
#  		 	breaks = function(x) { log(c(1, 0.2, 0.02, 0.002, 0.0003)) },
#   			labels = function(x) { c("1", ".2", ".02", ".002", ".0003") }
#  		 ) +
#  		 coord_cartesian(xlim=c(0,6.5), ylim = c(-8.5, 1.1)) +
#  		 annotate(geom="text", x=log(200), y=log(.001), label=paste0("t = ", r2b, " p = ", pb), color=cols[1]) +
#  		 annotate(geom="text", x=log(200), y=log(.0003), label=paste0("t = ", r2l, " p = ", pl), color=cols[2])
# 		 
# pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Breast_Lung_loglog.pdf", width=5.5, height=6)
# print(plot.0)
# dev.off()
# 
# r2b = signif(summary(lm(ctdna_frac ~ tot_volu, data=tmp %>% filter(Tissue=="Breast")))$r.squared, 3)
# pb = signif(summary(lm(ctdna_frac ~ tot_volu, data=tmp %>% filter(Tissue=="Breast")))$coefficients[2,4], 3)
# r2l = signif(summary(lm(ctdna_frac ~ tot_volu, data=tmp %>% filter(Tissue=="Lung")))$r.squared, 3)
# pl = signif(summary(lm(ctdna_frac ~ tot_volu, data=tmp %>% filter(Tissue=="Lung")))$coefficients[2,4], 3)
# 
# r2b = signif(cor(tmp$ctdna_frac[tmp$Tissue=="Breast"], tmp$tot_volu[tmp$Tissue=="Breast"], method="kendall"), 3)
# pb = signif(cor.test(tmp$ctdna_frac[tmp$Tissue=="Breast"], tmp$tot_volu[tmp$Tissue=="Breast"], method="kendall")$p.value, 3)
# r2l = signif(cor(tmp$ctdna_frac[tmp$Tissue=="Lung"], tmp$tot_volu[tmp$Tissue=="Lung"], method="kendall"), 3)
# pl = signif(cor.test(tmp$ctdna_frac[tmp$Tissue=="Lung"], tmp$tot_volu[tmp$Tissue=="Lung"], method="kendall")$p.value, 3)
# 
# plot.0 = ggplot(tmp, aes(x = tot_volu, y = ctdna_frac, fill = Tissue)) +
# 		 geom_point(alpha = .8, size = 2.5, color = "black", shape = 21) +
# 		 geom_smooth(method = "lm", se=FALSE, data=tmp, aes(x = tot_volu, y = ctdna_frac, color = Tissue)) +
# 		 scale_fill_manual(values = cols) +
# 		 theme_bw(base_size=15) +
# 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.9), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# 		 labs(x="\nDisease volume\n", y="\nctDNA fraction\n") +
# 		 facet_wrap(~facet) +
# 		 annotate(geom="text", x=500, y=.5, label=paste0("t = ", r2b, " p = ", pb), color=cols[1]) +
#  		 annotate(geom="text", x=500, y=.4, label=paste0("t = ", r2l, " p = ", pl), color=cols[2])
# 		 
# pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Breast_Lung_nolognolog.pdf", width=5.5, height=6)
# print(plot.0)
# dev.off()
# 
# r2b = signif(summary(lm(ctdna_frac ~ ln_vol, data=tmp %>% filter(Tissue=="Breast")))$r.squared, 3)
# pb = signif(summary(lm(ctdna_frac ~ ln_vol, data=tmp %>% filter(Tissue=="Breast")))$coefficients[2,4], 3)
# r2l = signif(summary(lm(ctdna_frac ~ ln_vol, data=tmp %>% filter(Tissue=="Lung")))$r.squared, 3)
# pl = signif(summary(lm(ctdna_frac ~ ln_vol, data=tmp %>% filter(Tissue=="Lung")))$coefficients[2,4], 3)
# 
# r2b = signif(cor(tmp$ctdna_frac[tmp$Tissue=="Breast"], tmp$ln_vol[tmp$Tissue=="Breast"], method="kendall"), 3)
# pb = signif(cor.test(tmp$ctdna_frac[tmp$Tissue=="Breast"], tmp$ln_vol[tmp$Tissue=="Breast"], method="kendall")$p.value, 3)
# r2l = signif(cor(tmp$ctdna_frac[tmp$Tissue=="Lung"], tmp$ln_vol[tmp$Tissue=="Lung"], method="kendall"), 3)
# pl = signif(cor.test(tmp$ctdna_frac[tmp$Tissue=="Lung"], tmp$ln_vol[tmp$Tissue=="Lung"], method="kendall")$p.value, 3)
# 
# plot.0 = ggplot(tmp, aes(x = ln_vol, y = ctdna_frac, fill = Tissue)) +
# 		 geom_point(alpha = .8, size = 2.5, color = "black", shape = 21) +
# 		 geom_smooth(method = "lm", se=FALSE, data=tmp, aes(x = ln_vol, y = ctdna_frac, color = Tissue)) +
# 		 scale_fill_manual(values = cols) +
# 		 theme_bw(base_size=15) +
# 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.9), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# 		 labs(x="\nDisease volume\n", y="\nctDNA fraction\n") +
# 		 facet_wrap(~facet) +
# 		 scale_x_continuous(
# 		 	breaks = function(x) { log(c(0, .1, 1, 10, 100, 500)) },
#  			labels = function(x) { c(0, .1, 1, 10, 100, 500) }
# 		 ) +
#  		 coord_cartesian(xlim=c(0,6.5)) +
#  		 annotate(geom="text", x=log(5), y=.70, label=paste0("t = ", r2b, " p = ", pb), color=cols[1]) +
#  		 annotate(geom="text", x=log(5), y=.65, label=paste0("t = ", r2l, " p = ", pl), color=cols[2])
# 		 
# pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Breast_Lung_lognolog.pdf", width=5.5, height=6)
# print(plot.0)
# dev.off()
# 
# r2b = signif(summary(lm(ln_frac ~ tot_volu, data=tmp %>% filter(Tissue=="Breast")))$r.squared, 3)
# pb = signif(summary(lm(ln_frac ~ tot_volu, data=tmp %>% filter(Tissue=="Breast")))$coefficients[2,4], 3)
# r2l = signif(summary(lm(ln_frac ~ tot_volu, data=tmp %>% filter(Tissue=="Lung")))$r.squared, 3)
# pl = signif(summary(lm(ln_frac ~ tot_volu, data=tmp %>% filter(Tissue=="Lung")))$coefficients[2,4], 3)
# 
# r2b = signif(cor(tmp$ln_frac[tmp$Tissue=="Breast"], tmp$tot_volu[tmp$Tissue=="Breast"], method="kendall"), 3)
# pb = signif(cor.test(tmp$ln_frac[tmp$Tissue=="Breast"], tmp$tot_volu[tmp$Tissue=="Breast"], method="kendall")$p.value, 3)
# r2l = signif(cor(tmp$ln_frac[tmp$Tissue=="Lung"], tmp$tot_volu[tmp$Tissue=="Lung"], method="kendall"), 3)
# pl = signif(cor.test(tmp$ln_frac[tmp$Tissue=="Lung"], tmp$tot_volu[tmp$Tissue=="Lung"], method="kendall")$p.value, 3)
# 
# plot.0 = ggplot(tmp, aes(x = tot_volu, y = ln_frac, fill = Tissue)) +
# 		 geom_point(alpha = .8, size = 2.5, color = "black", shape = 21) +
# 		 geom_smooth(method = "lm", se=FALSE, data=tmp, aes(x = tot_volu, y = ln_frac, color = Tissue)) +
# 		 scale_fill_manual(values = cols) +
# 		 theme_bw(base_size=15) +
# 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.9), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# 		 labs(x="\nDisease volume\n", y="\nctDNA fraction\n") +
# 		 facet_wrap(~facet) +
#  		 scale_y_continuous(
#  		 	breaks = function(x) { log(c(1, 0.2, 0.02, 0.002, 0.0003)) },
#   			labels = function(x) { c("1", ".2", ".02", ".002", ".0003") }
#  		 ) +
#  		 coord_cartesian(ylim = c(-8.5, 1.1)) +
#  		 annotate(geom="text", x=500, y=log(0.001), label=paste0("t = ", r2b, " p = ", pb), color=cols[1]) +
#  		 annotate(geom="text", x=500, y=log(0.0003), label=paste0("t = ", r2l, " p = ", pl), color=cols[2])
# 		 
# pdf(file="../res/rebuttal/cfDNA_Fraction_vs_Tv_Breast_Lung_nologlog.pdf", width=5.5, height=6)
# print(plot.0)
# dev.off()
