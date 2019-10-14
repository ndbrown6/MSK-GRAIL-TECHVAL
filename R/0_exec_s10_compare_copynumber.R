#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/figureS10")) {
	dir.create("../res/figureS10")
}

if (!dir.exists("../res/etc/Source_Data_Extended_Data_Fig_9")) {
	dir.create("../res/etc/Source_Data_Extended_Data_Fig_9")
}

'absolute_' <- function(rho, psi, gamma=1, x) {
	rho = ifelse(is.na(rho), 1, rho)
	psi = ifelse(is.na(psi), 2, psi)
	return(invisible(((((2^(x/gamma))*(rho*psi+(1-rho)*2)) - ((1-rho)*2))/rho)))
}
 
'prune_' <- function(x, n=10)
{
	cnm = matrix(NA, nrow=nrow(x), ncol=nrow(x))
	for (j in 1:nrow(x)) {
		cnm[,j] = abs(2^x[j,"log2"] - 2^x[,"log2"])
	}
	cnt = hclust(as.dist(cnm), "average")
	cnc = cutree(tree=cnt, k=n)
	for (j in unique(cnc)) {
		indx = which(cnc==j)
		if (length(indx)>2) {
 			mcl = mean(x[indx,"log2"])
			scl = sd(x[indx,"log2"])
			ind = which(x[indx,"log2"]<(mcl+1.96*scl) & x[indx,"log2"]>(mcl-1.96*scl))
			x[indx[ind],"log2"] = mean(x[indx[ind],"log2"])
		} else {
			x[indx,"log2"] = mean(x[indx,"log2"])
		}
	}
	return(x)
}

'plot_log3_' <- function(x, y, axis = TRUE, ylim=c(-2.5, 2.5))
{
   	par(mar=c(6.1, 9.5, 4.1, 1.1))
   	data(CytoBand)
   	end = NULL
	for (j in 1:23) {
		end = c(end, max(CytoBand$End[CytoBand$Chromosome==j]))
	}
	end = cumsum(end)
	start = rep(0, 23)
	start[2:23] = end[1:22]+1
	for (j in 1:23) {
		y[y[,"Chromosome"]==j,"Start"] = y[y[,"Chromosome"]==j,"Start"] + start[j]
		y[y[,"Chromosome"]==j,"End"] = y[y[,"Chromosome"]==j,"End"] + start[j]
		x[x[,"Chromosome"]==j,"Position"] = x[x[,"Chromosome"]==j,"Position"] + start[j]
	}
	plot(x[,"Position"], x[,"Log2Ratio"], type="p", pch=".", cex=1, col="grey75", axes=FALSE, frame=FALSE, xlab="", ylab="", main="", ylim=ylim)
	for (j in 1:nrow(y)) {
 		lines(x=c(y[j,"Start"], y[j,"End"]), y=rep(y[j,"Log2Ratio"],2), lty=1, lwd=2.5, col="red")
 	}
  	axis(2, at = c(-2, -1, 0, 1, 2), labels = c(-2, -1, 0, 1, 2), cex.axis = 1, las = 1)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3, cex = 1.15, las=1)
	points(c(-10000000, max(x[,"Position"])), c(0,0), type="l", col="black", lty=1, lwd=1)
	if (axis) {
		axis(side=1, at=c(start, end[length(end)]), labels=rep("", length(start)+1), tcl=.5)
		axis(1, at = .5*(start+end), labels=c(1:22, "X"), tcl=-.5, lwd=0, lwd.ticks=1, tcl=-.25)
	}

}

'undo_' <- function(x, n=10)
{
	cnm = matrix(NA, nrow=nrow(x), ncol=nrow(x))
	for (j in 1:nrow(x)) {
		cnm[,j] = abs(2^x[j,"Log2Ratio"] - 2^x[,"Log2Ratio"])
	}
	cnt = hclust(as.dist(cnm), "average")
	cnc = cutree(tree=cnt, k=n)
	for (j in unique(cnc)) {
		indx = which(cnc==j)
		if (length(indx)>2) {
 			mcl = mean(x[indx,"Log2Ratio"])
			scl = sd(x[indx,"Log2Ratio"])
			ind = which(x[indx,"Log2Ratio"]<(mcl+1.96*scl) & x[indx,"Log2Ratio"]>(mcl-1.96*scl))
			x[indx[ind],"Log2Ratio"] = mean(x[indx[ind],"Log2Ratio"])
		} else {
			x[indx,"Log2Ratio"] = mean(x[indx,"Log2Ratio"])
		}
	}
	return(x)
}

'ploidy_' <- function(x, y, w)
{
	psi = seq(from=1.5, to=3.5, length=100)
	sse = vector(mode="numeric", length=length(psi))
	for (i in 1:length(psi)) {
		z = absolute_(rho=y, psi=psi[i], gamma=.85, x=x)
		sse[i] = sum(((w/sum(w)) * (z-round(z)))^2)
	}
	return(sse)
}

#==================================================
# % agreement based on log2 ratios
#==================================================
key_file = read_tsv(file=url_master_key, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
 		   dplyr::select(PATIENT_ID, GRAIL_ID, DMP_ID, TUMOR_ID, NORMAL_ID, GRAIL_alpha, GRAIL_psi, IMPACT_alpha, IMPACT_psi)
		   
ctDNA_fraction = read_csv(file=url_ctdna_frac, col_types = cols(.default = col_character())) %>%
		   		 type_convert() %>%
		   		 rename(GRAIL_ID = ID)
		   		 
key_file = left_join(key_file, ctDNA_fraction, by="GRAIL_ID")
 
res = foreach (i=1:nrow(key_file)) %dopar% {
	cat(key_file$GRAIL_ID[i], "\n")
	impact_path = paste0("../res/rebuttal/msk_impact/cnvkit/totalcopy/", key_file$TUMOR_ID[i], ".RData")
	impact_data = new.env()
	load(impact_path, envir=impact_data)

	impact_seg = impact_data$tmp %>%
				 dplyr::select(chrom=Chromosome, start = Start, end = End, log2 = Log2Ratio, n=N) %>%
				 filter(chrom<23)
	impact_seg = prune_(x=impact_seg) %>%
				 bind_cols(cn = absolute_(rho=key_file$IMPACT_alpha[i],
										  psi=key_file$IMPACT_psi[i],
										  x=impact_seg$log2)) %>%
				 filter(n>=50) %>%
				 mutate(n = cumsum(n))
		
	Chromosome = impact_seg[,"chrom"]
	Start = impact_seg[,"start"]
	End = impact_seg[,"end"]
	Calls = impact_seg[,"log2"]
	im = data.frame(Chromosome, Start, End, Calls)

	grail_path = paste0("../res/rebuttal/uncollapsed_bam/cnvkit/totalcopy/", key_file$GRAIL_ID[i], "-T.RData")
	grail_data = new.env()
	load(grail_path, envir=grail_data)

	grail_seg = grail_data$tmp %>%
				dplyr::select(chrom=Chromosome, start = Start, end = End, log2 = Log2Ratio, n=N) %>%
				filter(chrom<23)


	grail_seg = prune_(x=grail_seg) %>%
				bind_cols(cn = absolute_(rho=key_file$GRAIL_alpha[i],
										 psi=key_file$GRAIL_psi[i],
										 x=grail_seg$log2)) %>%
				filter(n>=50) %>%
				mutate(n = cumsum(n))

	Chromosome = grail_seg[,"chrom"]
	Start = grail_seg[,"start"]
	End = grail_seg[,"end"]
	Calls = grail_seg[,"log2"]
	gr = data.frame(Chromosome, Start, End, Calls)

	biopsy_gr <- im %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Calls = Calls)
	cfdna_gr <- gr %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Calls = Calls)
	fo <- findOverlaps(biopsy_gr, cfdna_gr)
	df <- data.frame(cfDNA=mcols(cfdna_gr)[subjectHits(fo),], Biopsy=mcols(biopsy_gr)[queryHits(fo),])
	return(invisible(df))
}

r = unlist(foreach (i=1:nrow(key_file)) %dopar% {
	cat(key_file$GRAIL_ID[i], "\n")
	x = res[[i]][,"cfDNA"]
	y = res[[i]][,"Biopsy"]
	k = cor(x, y, method="pearson")
	return(k)
})

tmp.0 = data.frame(correlation_coefficient = r,
				   ctdna_fraction = key_file$ctdna_frac,
				   sample_id = key_file$GRAIL_ID,
				   tumor_purity = key_file$IMPACT_alpha) %>%
				   mutate(tissue = case_when(
						grepl("VB", sample_id) ~ "Breast",
						grepl("VL", sample_id) ~ "Lung",
						grepl("VP", sample_id) ~ "Prostate",
				   )) %>%
				   mutate(ctdna_fraction_cat = case_when(
						is.na(ctdna_fraction) ~ "NE",
						ctdna_fraction >= 0  & ctdna_fraction < .1 ~ "1-9",
						ctdna_fraction >= .1 & ctdna_fraction < .3 ~ "10-29",
						ctdna_fraction >= .3 & ctdna_fraction < .6 ~ "30-59",
						ctdna_fraction >= .6 & ctdna_fraction <= 1 ~ "60-100"
				   )) %>%
				   mutate(tissue = as.factor(tissue)) %>%
				   mutate(ctdna_fraction_cat = factor(ctdna_fraction_cat, levels=c("NE", "1-9", "10-29", "30-59", "60-100"), ordered=TRUE)) %>%
				   mutate(ctdna_fraction_w_ne = ifelse(is.na(ctdna_fraction), -.1, ctdna_fraction)) %>%
				   mutate(ctdna_fraction_cat_w_ne = ifelse(is.na(ctdna_fraction), "No", "Yes")) %>%
				   mutate(ctdna_fraction_cat_w_ne = factor(ctdna_fraction_cat_w_ne)) %>%
				   mutate(tumor_purity_cat = case_when(
						tumor_purity >= 0  & tumor_purity < .3 ~ "0-29",
						tumor_purity >= .3 & tumor_purity < .4 ~ "30-39",
						tumor_purity >= .4 & tumor_purity < .6 ~ "40-59",
						tumor_purity >= .6 & tumor_purity < .8 ~ "60-79",
						tumor_purity >= .8 ~ "80-100"
				   )) %>%
				   mutate(tumor_purity_cat = factor(tumor_purity_cat, ordered=TRUE, levels=c("0-29", "30-39", "40-59", "60-79", "80-100")))
		   
p1 = jonckheere.test(x=tmp.0$correlation_coefficient, g=tmp.0$ctdna_fraction_cat, alternative = "increasing", nperm=10000)

plot.0 = ggplot(tmp.0, aes(x = ctdna_fraction_cat, y = correlation_coefficient)) + 
		 geom_boxplot(outlier.shape=NA, width=.5, color="#de2d26", fill="white") +
		 theme_classic(base_size=15) +
		 coord_cartesian(ylim = c(-.25,1.1)) +
		 scale_y_continuous(
			breaks = function(x) { c(-.2, 0, .2, .4, .6, .8, 1) },
			labels = function(x) { c(-.2, 0, .2, .4, .6, .8, 1) }
		 ) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15)) +
		 labs(x="\nctDNA fraction (%)", y="\n\nCorrelation coefficient\n") +
		 annotate(geom="text", x=3, y=1, label=paste0("p = ",p1$p.value), size=5)
 
pdf(file="../res/figureS10/ctDNA_Fraction_versus_Correlation_Coefficient_Log2.pdf", width=6.5, height=6)
print(plot.0)
dev.off()
			   
p2 = jonckheere.test(x=tmp.0$correlation_coefficient, g=tmp.0$tumor_purity_cat, alternative = "increasing", nperm=10000)

plot.0 = ggplot(tmp.0, aes(x = tumor_purity_cat, y = correlation_coefficient)) + 
		 geom_boxplot(outlier.shape=NA, width=.5, color="#756bb1", fill="white") +
		 theme_classic(base_size=15) +
		 coord_cartesian(ylim = c(-.25,1.1)) +
		 scale_y_continuous(
			breaks = function(x) { c(-.2, 0, .2, .4, .6, .8, 1) },
			labels = function(x) { c(-.2, 0, .2, .4, .6, .8, 1) }
		 ) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15)) +
		 labs(x="\nTumor purity (%)", y="\n\nCorrelation coefficient\n") +
		 annotate(geom="text", x=3, y=1, label=paste0("p = ",p2$p.value), size=5)
 
pdf(file="../res/figureS10/Tumor_Purity_versus_Correlation_Coefficient_Log2.pdf", width=6.5, height=6)
print(plot.0)
dev.off()

export_x = tmp.0 %>%
		   dplyr::select(`patient_id` = `sample_id`,
		   				 `tissue` = `tissue`,
		   				 `corr_coeff` = `correlation_coefficient`,
		   				 `tumor_purity` = `tumor_purity`,
		   				 `ctdna_fraction` = `ctdna_fraction`,
		   				 `ctdna_fraction_cat` = `ctdna_fraction_cat`,
		   				 `ctdna_fraction_ne` = `ctdna_fraction_w_ne`,
		   				 `ctdna_fraction_cat_ne` = `ctdna_fraction_cat_w_ne`)
write_tsv(export_x, path="../res/etc/Source_Data_Extended_Data_Fig_9/Extended_Data_Fig_9f.tsv", append=FALSE, col_names=TRUE)

#==================================================
# log2 ratio plots grail cfdna tumor samples
#==================================================
load("../res/rebuttal/uncollapsed_bam/cnvkit/totalcopy/MSK-VB-0008-T.RData")
tmp2 = winsorize(CN, method="mad", tau=4.5, verbose=FALSE)
colnames(tmp2) = c("Chromosome","Position","Log2Ratio")
tmp2[tmp2[,"Chromosome"]==11,"Log2Ratio"] = CN[CN[,"Chromosome"]==11,"Log2Ratio"]
tmp = list()
for (i in 1:23) {
	if (i==11) {
		tmp[[i]] = pcf(data=tmp2[tmp2$Chromosome==i,,drop=FALSE], kmin=50, gamma=50, verbose=FALSE)[,2:7,drop=FALSE]
		colnames(tmp[[i]]) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
	} else {
		tmp[[i]] = pcf(data=tmp2[tmp2$Chromosome==i,,drop=FALSE], kmin=100, gamma=150, verbose=FALSE)[,2:7,drop=FALSE]
		colnames(tmp[[i]]) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
	}
}
tmp = do.call(rbind, tmp)
tmp = undo_(tmp, n=5)
pdf(file="../res/figureS10/MSK-VB-0008_cfDNA.pdf", width=9.25, height=4)
plot_log3_(x=tmp2, y=tmp, axis=TRUE)
dev.off()

export_x = tmp2 %>%
		   dplyr::select(chromosome = `Chromosome`,
		   				 position = `Position`,
		   				 log2ratio = `Log2Ratio`)
export_y = tmp %>%
		   dplyr::select(chromosome = `Chromosome`,
		   				 arm = `Arm`,
		   				 start = `Start`,
		   				 end = `End`,
		   				 n = `N`,
		   				 log2ratio = `Log2Ratio`)
write_tsv(export_x, path="../res/etc/Source_Data_Extended_Data_Fig_9/Extended_Data_Fig_9c_3.tsv", append=FALSE, col_names=TRUE)
write_tsv(export_y, path="../res/etc/Source_Data_Extended_Data_Fig_9/Extended_Data_Fig_9c_4.tsv", append=FALSE, col_names=TRUE)

load("../res/rebuttal/uncollapsed_bam/cnvkit/totalcopy/MSK-VL-0056-T.RData")
tmp2 = winsorize(CN, method="mad", tau=3.5, verbose=FALSE)
colnames(tmp2) = c("Chromosome","Position","Log2Ratio")
tmp = undo_(tmp, n=2)
pdf(file="../res/figureS10/MSK-VL-0056_cfDNA.pdf", width=9.25, height=4)
plot_log3_(x=tmp2, y=tmp, axis=TRUE)
dev.off()

export_x = tmp2 %>%
		   dplyr::select(`chromosome` = `Chromosome`,
		   				 `position` = `Position`,
		   				 `log2ratio` = `Log2Ratio`)
export_y = tmp %>%
		   dplyr::select(`chromosome` = `Chromosome`,
		   				 `arm` = `Arm`,
		   				 `start` = `Start`,
		   				 `end` = `End`,
		   				 `n` = `N`,
		   				 `log2ratio` = `Log2Ratio`)
write_tsv(export_x, path="../res/etc/Source_Data_Extended_Data_Fig_9/Extended_Data_Fig_9d_3.tsv", append=FALSE, col_names=TRUE)
write_tsv(export_y, path="../res/etc/Source_Data_Extended_Data_Fig_9/Extended_Data_Fig_9d_4.tsv", append=FALSE, col_names=TRUE)

load("../res/rebuttal/uncollapsed_bam/cnvkit/totalcopy/MSK-VP-0004-T.RData")
tmp2 = winsorize(CN, method="mad", tau=4.5, verbose=FALSE)
colnames(tmp2) = c("Chromosome","Position","Log2Ratio")
tmp = list()
for (i in 1:23) {
	if (i==13) {
		tmp[[i]] = pcf(data=tmp2[tmp2$Chromosome==i,,drop=FALSE], kmin=10, gamma=40, verbose=FALSE)[,2:7,drop=FALSE]
		colnames(tmp[[i]]) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
	} else {
		tmp[[i]] = pcf(data=tmp2[tmp2$Chromosome==i,,drop=FALSE], kmin=50, gamma=50, verbose=FALSE)[,2:7,drop=FALSE]
		colnames(tmp[[i]]) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
	}
}
tmp = do.call(rbind, tmp)
tmp = undo_(tmp, n=9)
pdf(file="../res/figureS10/MSK-VP-0004_cfDNA.pdf", width=9.25, height=4)
plot_log3_(x=tmp2, y=tmp, axis=TRUE)
dev.off()

export_x = tmp2 %>%
		   dplyr::select(`chromosome` = `Chromosome`,
		   				 `position` = `Position`,
		   				 `log2ratio` = `Log2Ratio`)
export_y = tmp %>%
		   dplyr::select(`chromosome` = `Chromosome`,
		   				 `arm` = `Arm`,
		   				 `start` = `Start`,
		   				 `end` = `End`,
		   				 `n` = `N`,
		   				 `log2ratio` = `Log2Ratio`)
write_tsv(export_x, path="../res/etc/Source_Data_Extended_Data_Fig_9/Extended_Data_Fig_9e_3.tsv", append=FALSE, col_names=TRUE)
write_tsv(export_y, path="../res/etc/Source_Data_Extended_Data_Fig_9/Extended_Data_Fig_9e_4.tsv", append=FALSE, col_names=TRUE)

#==================================================
# log2 ratio plots msk-impact tumor samples
#==================================================
key_file = read_tsv(file=url_master_key, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   dplyr::select(PATIENT_ID, GRAIL_ID, DMP_ID, TUMOR_ID, NORMAL_ID, GRAIL_alpha, GRAIL_psi, IMPACT_alpha, IMPACT_psi) %>%
		   filter(GRAIL_ID %in% c("MSK-VB-0008", "MSK-VL-0056", "MSK-VP-0004"))

cat(key_file$GRAIL_ID[1], "\n")
load(paste0("../res/rebuttal/msk_impact/cnvkit/totalcopy/", key_file$TUMOR_ID[1], ".RData"))
tmp2 = winsorize(CN, method="mad", tau=4.5, verbose=FALSE)
colnames(tmp2) = c("Chromosome","Position","Log2Ratio")
tmp = list()
for (ii in 1:23) {
	if (ii==13) {
		tmp[[ii]] = pcf(data=tmp2[tmp2$Chromosome==ii,,drop=FALSE], kmin=10, gamma=40, verbose=FALSE)[,2:7,drop=FALSE]
		colnames(tmp[[ii]]) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
	} else {
		tmp[[ii]] = pcf(data=tmp2[tmp2$Chromosome==ii,,drop=FALSE], kmin=10, gamma=40, verbose=FALSE)[,2:7,drop=FALSE]
		colnames(tmp[[ii]]) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
	}
}
tmp = do.call(rbind, tmp)
tmp = undo_(tmp, n=9)
pdf(file=paste0("../res/figureS10/", key_file$GRAIL_ID[1], "_Tumor.pdf"), width=9.25, height=4)
plot_log3_(x=tmp2, y=tmp, axis=FALSE, ylim=c(-4,4))
dev.off()

export_x = tmp2 %>%
		   dplyr::select(`chromosome` = `Chromosome`,
		   				 `position` = `Position`,
		   				 `log2ratio` = `Log2Ratio`)
export_y = tmp %>%
		   dplyr::select(`chromosome` = `Chromosome`,
		   				 `arm` = `Arm`,
		   				 `start` = `Start`,
		   				 `end` = `End`,
		   				 `n` = `N`,
		   				 `log2ratio` = `Log2Ratio`)
write_tsv(export_x, path="../res/etc/Source_Data_Extended_Data_Fig_9/Extended_Data_Fig_9c_1.tsv", append=FALSE, col_names=TRUE)
write_tsv(export_y, path="../res/etc/Source_Data_Extended_Data_Fig_9/Extended_Data_Fig_9c_2.tsv", append=FALSE, col_names=TRUE)

cat(key_file$GRAIL_ID[2], "\n")
load(paste0("../res/rebuttal/msk_impact/cnvkit/totalcopy/", key_file$TUMOR_ID[2], ".RData"))
tmp2 = winsorize(CN, method="mad", tau=4.5, verbose=FALSE)
colnames(tmp2) = c("Chromosome","Position","Log2Ratio")
tmp = list()
for (ii in 1:23) {
	if (ii==13) {
		tmp[[ii]] = pcf(data=tmp2[tmp2$Chromosome==ii,,drop=FALSE], kmin=10, gamma=40, verbose=FALSE)[,2:7,drop=FALSE]
		colnames(tmp[[ii]]) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
	} else {
		tmp[[ii]] = pcf(data=tmp2[tmp2$Chromosome==ii,,drop=FALSE], kmin=10, gamma=40, verbose=FALSE)[,2:7,drop=FALSE]
		colnames(tmp[[ii]]) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
	}
}
tmp = do.call(rbind, tmp)
tmp = undo_(tmp, n=9)
pdf(file=paste0("../res/figureS10/", key_file$GRAIL_ID[2], "_Tumor.pdf"), width=9.25, height=4)
plot_log3_(x=tmp2, y=tmp, axis=FALSE, ylim=c(-4,4))
dev.off()

export_x = tmp2 %>%
		   dplyr::select(`chromosome` = `Chromosome`,
		   				 `position` = `Position`,
		   				 `log2ratio` = `Log2Ratio`)
export_y = tmp %>%
		   dplyr::select(`chromosome` = `Chromosome`,
		   				 `arm` = `Arm`,
		   				 `start` = `Start`,
		   				 `end` = `End`,
		   				 `n` = `N`,
		   				 `log2ratio` = `Log2Ratio`)
write_tsv(export_x, path="../res/etc/Source_Data_Extended_Data_Fig_9/Extended_Data_Fig_9d_1.tsv", append=FALSE, col_names=TRUE)
write_tsv(export_y, path="../res/etc/Source_Data_Extended_Data_Fig_9/Extended_Data_Fig_9d_2.tsv", append=FALSE, col_names=TRUE)

cat(key_file$GRAIL_ID[3], "\n")
load(paste0("../res/rebuttal/msk_impact/cnvkit/totalcopy/", key_file$TUMOR_ID[3], ".RData"))
tmp2 = winsorize(CN, method="mad", tau=4.5, verbose=FALSE)
colnames(tmp2) = c("Chromosome","Position","Log2Ratio")
tmp = list()
for (ii in 1:23) {
	if (ii==13) {
		tmp[[ii]] = pcf(data=tmp2[tmp2$Chromosome==ii,,drop=FALSE], kmin=10, gamma=40, verbose=FALSE)[,2:7,drop=FALSE]
		colnames(tmp[[ii]]) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
	} else {
		tmp[[ii]] = pcf(data=tmp2[tmp2$Chromosome==ii,,drop=FALSE], kmin=10, gamma=40, verbose=FALSE)[,2:7,drop=FALSE]
		colnames(tmp[[ii]]) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
	}
}
tmp = do.call(rbind, tmp)
tmp = undo_(tmp, n=9)
pdf(file=paste0("../res/figureS10/", key_file$GRAIL_ID[3], "_Tumor.pdf"), width=9.25, height=4)
plot_log3_(x=tmp2, y=tmp, axis=FALSE, ylim=c(-4,4))
dev.off()

export_x = tmp2 %>%
		   dplyr::select(`chromosome` = `Chromosome`,
		   				 `position` = `Position`,
		   				 `log2ratio` = `Log2Ratio`)
export_y = tmp %>%
		   dplyr::select(`chromosome` = `Chromosome`,
		   				 `arm` = `Arm`,
		   				 `start` = `Start`,
		   				 `end` = `End`,
		   				 `n` = `N`,
		   				 `log2ratio` = `Log2Ratio`)
write_tsv(export_x, path="../res/etc/Source_Data_Extended_Data_Fig_9/Extended_Data_Fig_9e_1.tsv", append=FALSE, col_names=TRUE)
write_tsv(export_y, path="../res/etc/Source_Data_Extended_Data_Fig_9/Extended_Data_Fig_9e_2.tsv", append=FALSE, col_names=TRUE)
