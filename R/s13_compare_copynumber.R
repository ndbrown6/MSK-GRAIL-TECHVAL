#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/rebuttal")) {
	dir.create("../res/rebuttal")
}

'smooth_' <- function(data)
{
	return(invisible(winsorize(data, method="pcf", k=50, gamma=50, verbose=FALSE)))
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

'fix_6' <- function(x, y, ix = NULL)
{
	if (!is.null(ix)) {
		start = y$start[ix]
		end = y$end[ix]
		y$log2[ix] = 0
		indx = which(x$chrom==6 & x$pos>=start & x$pos<=end)
		z = x$log2[indx]
		z = z + abs(median(z, na.rm=TRUE))
		x$log2[indx] = z
	} else {
		index = which(y$chrom==6)
		if (length(index)==2) {
			start = y$start[index[1]]
			end = y$end[index[1]]
			y$log2[index[1]] = 0
			indx = which(x$chrom==6 & x$pos>=start & x$pos<=end)
			z = x$log2[indx]
			z = z + abs(median(z, na.rm=TRUE))
			x$log2[indx] = z
		} else if (length(index)>=3) {
			start = y$start[index[2]]
			end = y$end[index[2]]
			y$log2[index[2]] = 0
			indx = which(x$chrom==6 & x$pos>=start & x$pos<=end)
			z = x$log2[indx]
			z = z + abs(median(z, na.rm=TRUE))
			x$log2[indx] = z
		}
	}
	return(list(x, y))
}

'fix_7' <- function(x)
{
	index = x$chrom==7
	x$cnlr[index] = rnorm(n=sum(index), mean=median(x$cnlr[index]), sd=sd(x$cnlr)/1.5)
	return(x)
}

'plot_log2_' <- function(x, y, n=10, purity, ploidy, fix_6=FALSE, ix=NULL, title = "")
{

	cn = x$jointseg %>%
		 select(chrom, pos = maploc, log2 = cnlr)
	seg = y$cncf %>%
		  select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
	if (fix_6) {
		fixed_cn = fix_6(cn, seg, ix)
		cn = fixed_cn[[1]]
		seg = fixed_cn[[2]]
	}
	seg = prune_(x=seg, n) %>%
		  mutate(n = cumsum(n))
	
   	par(mar=c(5, 5, 4, 2)+.1)
   	load("../res/etc/CytoBand.RData")
   	end = NULL
   	for (i in 1:23) {
   		end = c(end, max(CytoBand[CytoBand[,1]==i,"End"]))
   	}
   	end = cumsum(end)
   	start = c(1, end[1:22]+1)
   	CytoBand = cbind(start, end)
   	index = NULL
   	for (i in 1:23) {
   		index = c(index, seq(from = CytoBand[i, "start"], to=CytoBand[i, "end"], length=sum(cn$chrom==i)))
   	}
	plot(index, cn$log2, type="p", pch=".", cex=1.95, col="grey80", axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,5))
 	for (j in 1:nrow(seg)) {
 		if (j == 1) {
 			lines(x=c(1, index[seg[j,"n"]]), y=rep(seg[j,"log2"],2), lty=1, lwd=2.75, col="red")
 		} else {
 			lines(x=c(index[seg[j-1,"n"]], index[seg[j,"n"]]), y=rep(seg[j,"log2"],2), lty=1, lwd=2.75, col="red")
 		}
  	}
  	axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1, las = 1)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	abline(v=1, col=transparentRgb("goldenrod3", 255), lty=3, lwd=.5)
	for (j in 1:23) {
		abline(v=CytoBand[j,"end"], col=transparentRgb("goldenrod3", 255), lty=3, lwd=.5)
	}
	axis(1, at = .5*(CytoBand[,"start"]+CytoBand[,"end"]), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
    for (k in c(1, 2, 4, 8, 14)) {
		abline(h=(log2(((purity)*k + (1-purity)*2)/((purity)*ploidy + (1-purity)*2))), col=transparentRgb("brown", 155), lty=3)
	}
	rect(xleft=1-1e10, xright=CytoBand[23,"end"]+1e10, ybottom=4, ytop=6, col="lightgrey", border="black", lwd=1.5)
	title(main = title, line=-1, cex.main=.75, font.main=1)
    box(lwd=1.5)
}

#==================================================
# Update IMPACT tumor (alpha, psi)
#==================================================
key_file = read_tsv(file="../res/etc/master_sample_key.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   select(PATIENT_ID, GRAIL_ID, DMP_ID, TUMOR_ID, NORMAL_ID, GRAIL_alpha, GRAIL_psi, IMPACT_alpha, IMPACT_psi)

if (FALSE) { foreach (i=1:nrow(key_file)) %dopar% {
	print(key_file$GRAIL_ID[i])
	impact_path = paste0("../res/rebuttal/MSK-IMPACT/facets/cncf/", key_file$TUMOR_ID[i], "_", key_file$NORMAL_ID[i], ".Rdata")
	impact_data = new.env()
	load(impact_path, envir=impact_data)
	
	impact_cn = impact_data$out2$jointseg %>%
			    select(chrom, pos = maploc, log2 = cnlr)
	impact_seg = impact_data$fit$cncf %>%
				 select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
	impact_seg = prune_(x=impact_seg) %>%
				 bind_cols(cn = absolute_(rho=key_file$IMPACT_alpha[i],
										  psi=key_file$IMPACT_psi[i],
										  x=impact_seg$log2)) %>%
				mutate(n = cumsum(n))
				
	file_path = paste0("../res/rebuttal/MSK-IMPACT/facets/plots/ext/", key_file$GRAIL_ID[i], ".pdf")
										
	pdf(file_path, width=10, height=4)
	par(mar=c(5, 5, 4, 2)+.1)	
	plot(impact_cn$log2, type="p", pch=".", cex=1.5, col="grey70", axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-2.5,2.5))
	for (j in 1:nrow(impact_seg)) {
		if (j == 1) {
			lines(x=c(1, impact_seg[j,"n"]), y=rep(impact_seg[j,"log2"],2), lty=1, lwd=2.75, col="red")
		} else {
			lines(x=c(impact_seg[j-1,"n"], impact_seg[j,"n"]), y=rep(impact_seg[j,"log2"],2), lty=1, lwd=2.75, col="red")
		}
 	}
 	axis(2, at = NULL, cex.axis = 1.15, las = 1)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	abline(v=1, col=transparentRgb("goldenrod3", 255), lty=3, lwd=.5)
	abline(h=0, col="black")
	for (j in 1:23) {
		v = max(which(impact_seg[,"chrom"]==j))
		abline(v=impact_seg[v,"n"], col=transparentRgb("goldenrod3", 255), lty=3, lwd=.5)
	}
	start = end = 1:23
	for (j in 2:23) {
		start[j] = impact_seg[max(which(impact_seg[,"chrom"]==j-1)),"n"]
		end[j] = impact_seg[max(which(impact_seg[,"chrom"]==j)),"n"]
	}
	start[1] = 1
	end[1] = start[2]
	axis(1, at = .5*(start+end), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
	
    purity = key_file$IMPACT_alpha[i]
    ploidy = key_file$IMPACT_psi[i]
    for (k in c(1, 2, 3, 4, 8, 14)) {
		abline(h=(log2(((purity)*k + (1-purity)*2)/((purity)*ploidy + (1-purity)*2))), col=transparentRgb("brown", 155), lty=3)
	}
    box(lwd=1.5)
	dev.off()
} }

#==================================================
# Update GRAIL cfDNA (alpha, psi)
#==================================================
key_file = read_tsv(file="../res/etc/master_sample_key.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   select(PATIENT_ID, GRAIL_ID, DMP_ID, TUMOR_ID, NORMAL_ID, GRAIL_alpha, GRAIL_psi, IMPACT_alpha, IMPACT_psi)

if (FALSE) { foreach (i=1:nrow(key_file)) %dopar% {
	print(key_file$GRAIL_ID[i])
	grail_path = paste0("../res/rebuttal/GRAIL/facets/cncf/", key_file$GRAIL_ID[i], "_", key_file$GRAIL_ID[i], "-N.Rdata")
	grail_data = new.env()
	load(grail_path, envir=grail_data)
	
	grail_cn = grail_data$out2$jointseg %>%
			   select(chrom, pos = maploc, log2 = cnlr)
	grail_seg = grail_data$fit$cncf %>%
				select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
	
	fixed_cn = fix_6(grail_cn, grail_seg)
	grail_cn = fixed_cn[[1]]
	grail_seg = fixed_cn[[2]]
	
	grail_seg = prune_(x=grail_seg) %>%
				bind_cols(cn = absolute_(rho=key_file$GRAIL_alpha[i],
										 psi=key_file$GRAIL_psi[i],
										 x=grail_seg$log2)) %>%
				mutate(n = cumsum(n))
	
				
	file_path = paste0("../res/rebuttal/GRAIL/facets/plots/ext/", key_file$GRAIL_ID[i], ".pdf")
										
	pdf(file_path, width=10, height=4)
	par(mar=c(5, 5, 4, 2)+.1)	
	plot(grail_cn$log2, type="p", pch=".", cex=1.5, col="grey70", axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-2.5,2.5))
	for (j in 1:nrow(grail_seg)) {
		if (j == 1) {
			lines(x=c(1, grail_seg[j,"n"]), y=rep(grail_seg[j,"log2"],2), lty=1, lwd=2.75, col="red")
		} else {
			lines(x=c(grail_seg[j-1,"n"], grail_seg[j,"n"]), y=rep(grail_seg[j,"log2"],2), lty=1, lwd=2.75, col="red")
		}
 	}
 	axis(2, at = NULL, cex.axis = 1.15, las = 1)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	abline(v=1, col=transparentRgb("goldenrod3", 255), lty=3, lwd=.5)
	abline(h=0, col="black")
	for (j in 1:23) {
		v = max(which(grail_seg[,"chrom"]==j))
		abline(v=grail_seg[v,"n"], col=transparentRgb("goldenrod3", 255), lty=3, lwd=.5)
	}
	start = end = 1:23
	for (j in 2:23) {
		start[j] = grail_seg[max(which(grail_seg[,"chrom"]==j-1)),"n"]
		end[j] = grail_seg[max(which(grail_seg[,"chrom"]==j)),"n"]
	}
	start[1] = 1
	end[1] = start[2]
	axis(1, at = .5*(start+end), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
	
    purity = key_file$GRAIL_alpha[i]
    ploidy = key_file$GRAIL_psi[i]
    for (k in c(1, 2, 3, 4, 8, 14)) {
		abline(h=(log2(((purity)*k + (1-purity)*2)/((purity)*ploidy + (1-purity)*2))), col=transparentRgb("brown", 155), lty=3)
	}
    box(lwd=1.5)
	dev.off()
} }

#==================================================
# Updated MSK-IMPACT tumor Log2 Ratios
#==================================================
key_file = read_tsv(file="../res/etc/master_sample_key.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   select(PATIENT_ID, GRAIL_ID, DMP_ID, TUMOR_ID, NORMAL_ID, GRAIL_alpha, GRAIL_psi, IMPACT_alpha, IMPACT_psi)
		   
foreach (i=1:nrow(key_file)) %dopar% {
	print(key_file$GRAIL_ID[i])
	load(paste0("../res/rebuttal/MSK-IMPACT/facets/cncf/", key_file$TUMOR_ID[i], "_", key_file$NORMAL_ID[i], ".Rdata"))
	purity = key_file$IMPACT_alpha[i]
	ploidy = key_file$IMPACT_psi[i]
	pdf(file=paste0("../res/rebuttal/MSK-IMPACT/facets/plots/ext/", key_file$GRAIL_ID[i], ".pdf"), width=10, height=4.25)
	plot_log2_(x=out2, y=fit, n=10, purity, ploidy, FALSE, ix=NULL, title = paste0(key_file$GRAIL_ID[i], " | Tumor | purity = ", signif(purity,2), " | ploidy = ", signif(ploidy,3)))
	dev.off()
	return(invisible(1))
}

#==================================================
# Updated GRAIL cfDNA Log2 Ratios
#==================================================
key_file = read_tsv(file="../res/etc/master_sample_key.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   select(PATIENT_ID, GRAIL_ID, DMP_ID, TUMOR_ID, NORMAL_ID, GRAIL_alpha, GRAIL_psi, IMPACT_alpha, IMPACT_psi)
		   
foreach (i=1:nrow(key_file)) %dopar% {
	print(key_file$GRAIL_ID[i])
	load(paste0("../res/rebuttal/GRAIL/facets/cncf/", key_file$GRAIL_ID[i], "_", key_file$GRAIL_ID[i], "-N.Rdata"))
	purity = key_file$GRAIL_alpha[i]
	ploidy = key_file$GRAIL_psi[i]
	pdf(file=paste0("../res/rebuttal/GRAIL/facets/plots/ext/", key_file$GRAIL_ID[i], ".pdf"), width=10, height=4.25)
	plot_log2_(x=out2, y=fit, n=4, purity, ploidy, TRUE, ix=NULL, title = paste0(key_file$GRAIL_ID[i], " | cfDNA | purity = ", signif(purity,2), " | ploidy = ", signif(ploidy,3)))
	dev.off()
	return(invisible(1))
}

#==================================================
# % agreement based on copy number calls
#==================================================
key_file = read_tsv(file="../res/etc/master_sample_key.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   select(PATIENT_ID, GRAIL_ID, DMP_ID, TUMOR_ID, NORMAL_ID, GRAIL_alpha, GRAIL_psi, IMPACT_alpha, IMPACT_psi)
		   
ctDNA_fraction = read_csv(file=url_ctdna, col_types = cols(.default = col_character())) %>%
		   		 type_convert() %>%
		   		 rename(GRAIL_ID = ID)
		   		 
key_file = left_join(key_file, ctDNA_fraction, by="GRAIL_ID")
		   
res = foreach (i=1:nrow(key_file)) %dopar% {
 	cat(key_file$GRAIL_ID[i], "\n")
 	impact_path = paste0("../res/rebuttal/MSK-IMPACT/facets/cncf/", key_file$TUMOR_ID[i], "_", key_file$NORMAL_ID[i], ".Rdata")
	impact_data = new.env()
	load(impact_path, envir=impact_data)
	
	impact_seg = impact_data$fit$cncf %>%
				 select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
	impact_seg = prune_(x=impact_seg) %>%
				 bind_cols(cn = absolute_(rho=key_file$IMPACT_alpha[i],
										  psi=key_file$IMPACT_psi[i],
										  x=impact_seg$log2)) %>%
				 filter(n>=50) %>%
				 mutate(n = cumsum(n))
				
 	Chromosome = impact_seg[,"chrom"]
 	Start = impact_seg[,"start"]
 	End = impact_seg[,"end"]
 	Calls = impact_seg[,"cn"]
 	im = data.frame(Chromosome, Start, End, Calls)
 	
 	grail_path = paste0("../res/rebuttal/GRAIL/facets/cncf/", key_file$GRAIL_ID[i], "_", key_file$GRAIL_ID[i], "-N.Rdata")
	grail_data = new.env()
	load(grail_path, envir=grail_data)
	
	grail_cn = grail_data$out2$jointseg %>%
			   select(chrom, pos = maploc, log2 = cnlr)
	grail_seg = grail_data$fit$cncf %>%
				select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
	
	fixed_cn = fix_6(grail_cn, grail_seg)
	grail_cn = fixed_cn[[1]]
	grail_seg = fixed_cn[[2]]
	
	grail_seg = prune_(x=grail_seg) %>%
				bind_cols(cn = absolute_(rho=key_file$GRAIL_alpha[i],
										 psi=key_file$GRAIL_psi[i],
										 x=grail_seg$log2)) %>%
			 	filter(n>=50) %>%
				mutate(n = cumsum(n))
 	
 	Chromosome = grail_seg[,"chrom"]
 	Start = grail_seg[,"start"]
 	End = grail_seg[,"end"]
 	Calls = grail_seg[,"cn"]
 	gr = data.frame(Chromosome, Start, End, Calls)

 	biopsy_gr <- im %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Calls = Calls)
 	cfdna_gr <- gr %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Calls = Calls)
 	fo <- findOverlaps(biopsy_gr, cfdna_gr)
 	df <- data.frame(cfDNA=mcols(cfdna_gr)[subjectHits(fo),], Biopsy=mcols(biopsy_gr)[queryHits(fo),])
 	return(invisible(round(df)))
}

pa = unlist(foreach (i=1:nrow(key_file)) %dopar% {
 	cat(key_file$GRAIL_ID[i], "\n")
 	x = res[[i]][,"cfDNA"] 
 	y = res[[i]][,"Biopsy"]
 	
 	x2 = rep(0, length(x))
 	y2 = rep(0, length(y))
 	
 	psi = round(key_file$IMPACT_psi[i])
 	psi = ifelse(psi==1, 2, psi)
 	if (psi==2) {
	 	x2[x<1] = -1
		x2[x>6] = 1
	} else if (psi==3) {
		x2[x<1] = -1
		x2[x>7] = 1
	} else {
		x2[x<2] = -1
		x2[x>8] = 1
	}
	
	psi = round(key_file$GRAIL_psi[i])
 	psi = ifelse(psi==1, 2, psi)
	if (psi==2) {
	 	y2[y<1] = -1
		y2[y>6] = 1
	} else if (psi==3) {
		y2[y<1] = -1
		y2[y>7] = 1
	} else {
		y2[y<2] = -1
		y2[y>8] = 1
	}
	
	class = c(-1, 0, 1)
	z = matrix(0, ncol=length(class), nrow=length(class))
	
	for (ii in 1:length(class)) {
		for (jj in 1:length(class)) {
			z[ii,jj] = sum(x2==class[ii] & y2==class[jj])
		}
	}
	k = sum(diag(z))/sum(z)
	return(k)
})

tmp.0 = data.frame(pa = 100*pa,
				   ctdna_f = key_file$ctdna_frac,
				   sample_id = key_file$GRAIL_ID) %>%
				   mutate(tissue = case_when(
				   		grepl("VB", sample_id) ~ "Breast",
				   		grepl("VL", sample_id) ~ "Lung",
				   		grepl("VP", sample_id) ~ "Prostate",
				   )) %>%
				   mutate(ctdna_cat = case_when(
				   		ctdna_f > 0 & ctdna_f < .2 ~ 1,
				   		ctdna_f >= .2 & ctdna_f < .4 ~ 2,
				   		ctdna_f >= .4 & ctdna_f < .6 ~ 3,
				   		ctdna_f >= .6 & ctdna_f < .8 ~ 4,
				   		ctdna_f >= .8 & ctdna_f < 1 ~ 5
				   )) %>%
				   mutate(tissue = as.factor(tissue)) %>%
				   mutate(ctdna_cat = as.factor(ctdna_cat))

plot.0 = ggplot(tmp.0, aes(x=pa)) + 
  		 geom_density(fill="lightgrey") +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="\nDensity\n", x="\nAgreement (%)\n") +
		 coord_cartesian(xlim=c(0,100))
pdf(file="../res/rebuttal/Distribution_Percent_Agreement_All_Calls.pdf", width=5, height=6)
print(plot.0)
dev.off()

plot.0 = ggplot(tmp.0, aes(x=ctdna_f, y=pa, fill = tissue)) +
		 geom_point(alpha=1, size = 2.5, shape=21) +
		 scale_fill_manual(values = c("salmon", "#FDAE61", "#ABDDA4")) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\nctDNA fraction\n\n", y="Agreement (%)\n") +
		 coord_cartesian(ylim=c(0, 100)) +
		 theme(legend.justification = c(1, 0),
		 	   legend.position = c(1, 0),
		 	   legend.title = element_blank(),
		 	   legend.background = element_blank(),
		 	   legend.text=element_text(size=8)) +
		 guides(fill=guide_legend(title=c("Cancer type"), override.aes=list(shape=21)))
pdf(file="../res/rebuttal/Percent_Agreement_ctDNA_Fraction_Calls.pdf", width=6, height=6)
print(plot.0)
dev.off()

plot.0 = ggplot(tmp.0 %>% filter(tissue == "Breast"), aes(x=pa)) + 
  		 geom_density(fill="salmon") +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="\nDensity\n", x="\nAgreement (%)\n") +
		 coord_cartesian(xlim=c(0,100))
pdf(file="../res/rebuttal/Distribution_Percent_Agreement_Breast_Calls.pdf", width=5, height=6)
print(plot.0)
dev.off()

plot.0 = ggplot(tmp.0 %>% filter(tissue == "Lung"), aes(x=pa)) + 
  		 geom_density(fill="#FDAE61") +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="\nDensity\n", x="\nAgreement (%)\n") +
		 coord_cartesian(xlim=c(0,100))
pdf(file="../res/rebuttal/Distribution_Percent_Agreement_Lung_Calls.pdf", width=5, height=6)
print(plot.0)
dev.off()

plot.0 = ggplot(tmp.0 %>% filter(tissue == "Prostate"), aes(x=pa)) + 
  		 geom_density(fill="#ABDDA4") +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="\nDensity\n", x="\nAgreement (%)\n") +
		 coord_cartesian(xlim=c(0,100))
pdf(file="../res/rebuttal/Distribution_Percent_Agreement_Prostate_Calls.pdf", width=5, height=6)
print(plot.0)
dev.off()

#==================================================
# % agreement based on absolute copy numbers
#==================================================
key_file = read_tsv(file="../res/etc/master_sample_key.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   select(PATIENT_ID, GRAIL_ID, DMP_ID, TUMOR_ID, NORMAL_ID, GRAIL_alpha, GRAIL_psi, IMPACT_alpha, IMPACT_psi)
		   
ctDNA_fraction = read_csv(file=url_ctdna, col_types = cols(.default = col_character())) %>%
		   		 type_convert() %>%
		   		 rename(GRAIL_ID = ID)
		   		 
key_file = left_join(key_file, ctDNA_fraction, by="GRAIL_ID")
		   
res = foreach (i=1:nrow(key_file)) %dopar% {
 	cat(key_file$GRAIL_ID[i], "\n")
 	impact_path = paste0("../res/rebuttal/MSK-IMPACT/facets/cncf/", key_file$TUMOR_ID[i], "_", key_file$NORMAL_ID[i], ".Rdata")
	impact_data = new.env()
	load(impact_path, envir=impact_data)
	
	impact_seg = impact_data$fit$cncf %>%
				 select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
	impact_seg = prune_(x=impact_seg) %>%
				 bind_cols(cn = absolute_(rho=key_file$IMPACT_alpha[i],
										  psi=key_file$IMPACT_psi[i],
										  x=impact_seg$log2)) %>%
				 filter(n>=50) %>%
				 mutate(n = cumsum(n))
				
 	Chromosome = impact_seg[,"chrom"]
 	Start = impact_seg[,"start"]
 	End = impact_seg[,"end"]
 	Calls = impact_seg[,"cn"]
 	im = data.frame(Chromosome, Start, End, Calls)
 	
 	grail_path = paste0("../res/rebuttal/GRAIL/facets/cncf/", key_file$GRAIL_ID[i], "_", key_file$GRAIL_ID[i], "-N.Rdata")
	grail_data = new.env()
	load(grail_path, envir=grail_data)
	
	grail_cn = grail_data$out2$jointseg %>%
			   select(chrom, pos = maploc, log2 = cnlr)
	grail_seg = grail_data$fit$cncf %>%
				select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
	
	fixed_cn = fix_6(grail_cn, grail_seg)
	grail_cn = fixed_cn[[1]]
	grail_seg = fixed_cn[[2]]
	
	grail_seg = prune_(x=grail_seg) %>%
				bind_cols(cn = absolute_(rho=key_file$GRAIL_alpha[i],
										 psi=key_file$GRAIL_psi[i],
										 x=grail_seg$log2)) %>%
				filter(n>=50) %>%
				mutate(n = cumsum(n))
 	
 	Chromosome = grail_seg[,"chrom"]
 	Start = grail_seg[,"start"]
 	End = grail_seg[,"end"]
 	Calls = grail_seg[,"cn"]
 	gr = data.frame(Chromosome, Start, End, Calls)

 	biopsy_gr <- im %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Calls = Calls)
 	cfdna_gr <- gr %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Calls = Calls)
 	fo <- findOverlaps(biopsy_gr, cfdna_gr)
 	df <- data.frame(cfDNA=mcols(cfdna_gr)[subjectHits(fo),], Biopsy=mcols(biopsy_gr)[queryHits(fo),])
 	return(invisible(round(df)))
}

pa = unlist(foreach (i=1:nrow(key_file)) %dopar% {
 	cat(key_file$GRAIL_ID[i], "\n")
 	x = res[[i]][,"cfDNA"]
 	x[x<0] = 0
 	y = res[[i]][,"Biopsy"]
 	y[y<0] = 0
 	
 	class = 0:max(c(x, y))
 	z = matrix(0, ncol=length(class), nrow=length(class))

	for (ii in 1:length(class)) {
 		for (jj in 1:length(class)) {
 			z[ii,jj] = sum(x==class[ii] & y==class[jj])
 		}
 	}
 	k = sum(diag(z))/sum(z)
 	return(k)
})

tmp.0 = data.frame(pa = 100*pa,
				   ctdna_f = key_file$ctdna_frac,
				   sample_id = key_file$GRAIL_ID) %>%
				   mutate(tissue = case_when(
				   		grepl("VB", sample_id) ~ "Breast",
				   		grepl("VL", sample_id) ~ "Lung",
				   		grepl("VP", sample_id) ~ "Prostate",
				   )) %>%
				   mutate(ctdna_cat = case_when(
				   		ctdna_f > 0 & ctdna_f < .2 ~ 1,
				   		ctdna_f >= .2 & ctdna_f < .4 ~ 2,
				   		ctdna_f >= .4 & ctdna_f < .6 ~ 3,
				   		ctdna_f >= .6 & ctdna_f < .8 ~ 4,
				   		ctdna_f >= .8 & ctdna_f < 1 ~ 5
				   )) %>%
				   mutate(tissue = as.factor(tissue)) %>%
				   mutate(ctdna_cat = as.factor(ctdna_cat))
				   
plot.0 = ggplot(tmp.0, aes(x=pa)) + 
  		 geom_density(fill="lightgrey") +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="\nDensity\n", x="\nAgreement (%)\n") +
		 coord_cartesian(xlim=c(0,100))
pdf(file="../res/rebuttal/Distribution_Percent_Agreement_All_CN.pdf", width=5, height=6)
print(plot.0)
dev.off()

plot.0 = ggplot(tmp.0, aes(x=ctdna_f, y=pa, fill = tissue)) +
		 geom_point(alpha=1, size = 2.5, shape=21) +
		 scale_fill_manual(values = c("salmon", "#FDAE61", "#ABDDA4")) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\nctDNA fraction\n\n", y="Agreement (%)\n") +
		 coord_cartesian(ylim=c(0, 100)) +
		 theme(legend.justification = c(1, 0),
		 	   legend.position = c(1, 0),
		 	   legend.title = element_blank(),
		 	   legend.background = element_blank(),
		 	   legend.text=element_text(size=8)) +
		 guides(fill=guide_legend(title=c("Cancer type"), override.aes=list(shape=21)))
pdf(file="../res/rebuttal/Percent_Agreement_ctDNA_Fraction_CN.pdf", width=6, height=6)
print(plot.0)
dev.off()

plot.0 = ggplot(tmp.0 %>% filter(tissue == "Breast"), aes(x=pa)) + 
  		 geom_density(fill = "salmon") +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="\nDensity\n", x="\nAgreement (%)\n") +
		 coord_cartesian(xlim=c(0,100))
pdf(file="../res/rebuttal/Distribution_Percent_Agreement_Breast_CN.pdf", width=5, height=6)
print(plot.0)
dev.off()

plot.0 = ggplot(tmp.0 %>% filter(tissue == "Lung"), aes(x=pa)) + 
  		 geom_density(fill = "#FDAE61") +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="\nDensity\n", x="\nAgreement (%)\n") +
		 coord_cartesian(xlim=c(0,100))
pdf(file="../res/rebuttal/Distribution_Percent_Agreement_Lung_CN.pdf", width=5, height=6)
print(plot.0)
dev.off()

plot.0 = ggplot(tmp.0 %>% filter(tissue == "Prostate"), aes(x=pa)) + 
  		 geom_density(fill = "#ABDDA4") +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="\nDensity\n", x="\nAgreement (%)\n") +
		 coord_cartesian(xlim=c(0,100))
pdf(file="../res/rebuttal/Distribution_Percent_Agreement_Prostate_CN.pdf", width=5, height=6)
print(plot.0)
dev.off()

#==================================================
# % agreement based on Log2 ratios
#==================================================
key_file = read_tsv(file="../res/etc/master_sample_key.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   select(PATIENT_ID, GRAIL_ID, DMP_ID, TUMOR_ID, NORMAL_ID, GRAIL_alpha, GRAIL_psi, IMPACT_alpha, IMPACT_psi)
		   
ctDNA_fraction = read_csv(file=url_ctdna, col_types = cols(.default = col_character())) %>%
		   		 type_convert() %>%
		   		 rename(GRAIL_ID = ID)
		   		 
key_file = left_join(key_file, ctDNA_fraction, by="GRAIL_ID")

res = foreach (i=1:nrow(key_file)) %dopar% {
 	cat(key_file$GRAIL_ID[i], "\n")
 	impact_path = paste0("../res/rebuttal/MSK-IMPACT/facets/cncf/", key_file$TUMOR_ID[i], "_", key_file$NORMAL_ID[i], ".Rdata")
	impact_data = new.env()
	load(impact_path, envir=impact_data)
	
	impact_seg = impact_data$fit$cncf %>%
				 select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
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
 	
 	grail_path = paste0("../res/rebuttal/GRAIL/facets/cncf/", key_file$GRAIL_ID[i], "_", key_file$GRAIL_ID[i], "-N.Rdata")
	grail_data = new.env()
	load(grail_path, envir=grail_data)
	
	grail_cn = grail_data$out2$jointseg %>%
			   select(chrom, pos = maploc, log2 = cnlr)
	grail_seg = grail_data$fit$cncf %>%
				select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
	
	fixed_cn = fix_6(grail_cn, grail_seg)
	grail_cn = fixed_cn[[1]]
	grail_seg = fixed_cn[[2]]
	
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

pa = unlist(foreach (i=1:nrow(key_file)) %dopar% {
 	cat(key_file$GRAIL_ID[i], "\n")
 	x = res[[i]][,"cfDNA"]
 	y = res[[i]][,"Biopsy"]
 	k = cor(x, y, method="pearson")
 	return(k)
})


tmp.0 = data.frame(pa = pa,
				   ctdna_f = key_file$ctdna_frac,
				   sample_id = key_file$GRAIL_ID) %>%
				   mutate(tissue = case_when(
				   		grepl("VB", sample_id) ~ "Breast",
				   		grepl("VL", sample_id) ~ "Lung",
				   		grepl("VP", sample_id) ~ "Prostate",
				   )) %>%
				   mutate(ctdna_cat = case_when(
				   		is.na(ctdna_f) ~ "0",
				   		ctdna_f > 0 & ctdna_f < .2 ~ "0-19",
				   		ctdna_f >= .2 & ctdna_f < .4 ~ "20-39",
				   		ctdna_f >= .4 & ctdna_f < .6 ~ "40-59",
				   		ctdna_f >= .6 & ctdna_f < .8 ~ "60-79",
				   		ctdna_f >= .8 & ctdna_f < 1 ~ "80-100"
				   )) %>%
				   mutate(tissue = as.factor(tissue)) %>%
				   mutate(ctdna_cat = as.factor(ctdna_cat))

plot.0 = ggplot(tmp.0, aes(x=pa)) + 
  		 geom_density(fill = "lightgrey") +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="\nDensity\n", x="\nPearson's r\n") +
		 coord_cartesian(xlim=c(-1,1))
pdf(file="../res/rebuttal/Distribution_Percent_Agreement_All_Log2.pdf", width=5, height=6)
print(plot.0)
dev.off()

plot.0 = ggplot(tmp.0, aes(x=ctdna_f, y=pa, fill = tissue)) +
		 geom_point(alpha=1, size = 2.5, shape=21) +
		 scale_fill_manual(values = c("salmon", "#FDAE61", "#ABDDA4")) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\nctDNA fraction\n\n", y="Pearson's r\n") +
		 coord_cartesian(ylim=c(-1,1)) +
		 theme(legend.justification = c(1, 0),
		 	   legend.position = c(1, 0),
		 	   legend.title = element_blank(),
		 	   legend.background = element_blank(),
		 	   legend.text=element_text(size=8)) +
		 guides(fill=guide_legend(title=c("Cancer type"), override.aes=list(shape=21)))
pdf(file="../res/rebuttal/Percent_Agreement_ctDNA_Fraction_Log2.pdf", width=6, height=6)
print(plot.0)
dev.off()

plot.0 = ggplot(tmp.0, aes(x=ctdna_cat, y=pa)) +
		 geom_boxplot(alpha=1, outlier.size=2.5, outlier.shape=21, fill="grey80") +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(x="\nctDNA fraction (%)\n\n", y=expression("Pearson"~r)) +
		 coord_cartesian(ylim=c(-1,1)) +
		 theme(legend.justification = c(1, 0),
		 	   legend.position = c(1, 0),
		 	   legend.title = element_blank(),
		 	   legend.background = element_blank(),
		 	   legend.text=element_text(size=8))
pdf(file="../res/rebuttal/Percent_Agreement_ctDNA_Fraction_Log2_b.pdf", width=6, height=6)
print(plot.0)
dev.off()

plot.0 = ggplot(tmp.0 %>% filter(tissue == "Breast"), aes(x=pa)) + 
  		 geom_density(fill = "salmon") +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="\nDensity\n", x="\nPearson's r\n") +
		 coord_cartesian(xlim=c(-1,1))
pdf(file="../res/rebuttal/Distribution_Percent_Agreement_Breast_Log2.pdf", width=5, height=6)
print(plot.0)
dev.off()

plot.0 = ggplot(tmp.0 %>% filter(tissue == "Lung"), aes(x=pa)) + 
  		 geom_density(fill = "#FDAE61") +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="\nDensity\n", x="\nPearson's r\n") +
		 coord_cartesian(xlim=c(-1,1))
pdf(file="../res/rebuttal/Distribution_Percent_Agreement_Lung_Log2.pdf", width=5, height=6)
print(plot.0)
dev.off()

plot.0 = ggplot(tmp.0 %>% filter(tissue == "Prostate"), aes(x=pa)) + 
  		 geom_density(fill = "#ABDDA4") +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="\nDensity\n", x="\nPearson's r\n") +
		 coord_cartesian(xlim=c(-1,1))
pdf(file="../res/rebuttal/Distribution_Percent_Agreement_Prostate_Log2.pdf", width=5, height=6)
print(plot.0)
dev.off()

# #==================================================
# # Comparison of copy number aberrations by gene
# #==================================================
# tracker = read_tsv(file="../res/etc/master_sample_key.tsv", col_types = cols(.default = col_character())) %>%
# 		  type_convert() %>%
# 		  select(PATIENT_ID, GRAIL_ID, DMP_ID, TUMOR_ID, NORMAL_ID, GRAIL_alpha, GRAIL_psi, IMPACT_alpha, IMPACT_psi)
# 
# i_bygene = foreach (i=1:nrow(tracker)) %dopar% {
#  	cat(tracker$GRAIL_ID[i], "\n")
#  	impact_path = paste0("../res/rebuttal/MSK-IMPACT/facets/cncf/", key_file$TUMOR_ID[i], "_", key_file$NORMAL_ID[i], ".Rdata")
# 	impact_data = new.env()
# 	load(impact_path, envir=impact_data)
# 	
# 	impact_seg = impact_data$fit$cncf %>%
# 				 select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
# 	impact_seg = prune_(x=impact_seg) %>%
# 				 bind_cols(cn = absolute_(rho=key_file$IMPACT_alpha[i],
# 										  psi=key_file$IMPACT_psi[i],
# 										  x=impact_seg$log2)) %>%
# 				 mutate(n = cumsum(n))
# 				
#  	Chromosome = impact_seg[,"chrom"]
#  	Start = impact_seg[,"start"]
#  	End = impact_seg[,"end"]
#  	Calls = impact_seg[,"cn"]
#  	res = data.frame(Chromosome, Start, End, Calls)
#  	annot = read.csv(file="~/share/reference/IMPACT410_genes_for_copynumber.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
#  			select(hgnc_symbol, chr, start_position, end_position) %>%
#  			rename(Hugo_Symbol = hgnc_symbol,
#  				   Chromosome = chr,
#  				   Start = start_position,
#  				   End = end_position) %>%
#  			mutate(Chromosome = ifelse(Chromosome %in% "X", 23, Chromosome)) %>%
#  			arrange(as.numeric(Chromosome), as.numeric(Start))
#  				   
#  	annot_by_gene <- annot %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Hugo_Symbol = Hugo_Symbol)
#  	res_by_segment <- res %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Calls = Calls)
#  	fo <- findOverlaps(res_by_segment, annot_by_gene)
#  
#  	df <- data.frame(Hugo_Symbol=mcols(annot_by_gene)[subjectHits(fo),], Calls=mcols(res_by_segment)[queryHits(fo),])
#  	Hugo_Symbol = which(duplicated(df$Hugo_Symbol))
#  	for (j in 1:length(Hugo_Symbol)) {
#  		index = which(as.character(df$Hugo_Symbol)==as.character(df$Hugo_Symbol[Hugo_Symbol[j]]))
#  		df[index,2] = mean(df[index,2], na.rm=TRUE)
#  	}
#  	df = df %>% filter(!duplicated(Hugo_Symbol))
#  	df[,2] = round(df[,2])
#  	res = rep(0, nrow(annot))
#  	names(res) = annot[,"Hugo_Symbol"]
#  	res[as.character(df$Hugo_Symbol)] = df$Calls
#  	return(invisible(res))
# }
# i_bygene = do.call(cbind, i_bygene)
# colnames(i_bygene) = tracker$GRAIL_ID
# 
#  
# g_bygene = foreach (i=1:nrow(tracker)) %dopar% {
#  	cat(tracker$GRAIL_ID[i], "\n")
# 	grail_path = paste0("../res/rebuttal/GRAIL/facets/cncf/", key_file$GRAIL_ID[i], "_", key_file$GRAIL_ID[i], "-N.Rdata")
# 	grail_data = new.env()
# 	load(grail_path, envir=grail_data)
# 	
# 	grail_cn = grail_data$out2$jointseg %>%
# 			   select(chrom, pos = maploc, log2 = cnlr)
# 	grail_seg = grail_data$fit$cncf %>%
# 				select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
# 	
# 	fixed_cn = fix_6(grail_cn, grail_seg)
# 	grail_cn = fixed_cn[[1]]
# 	grail_seg = fixed_cn[[2]]
# 	
# 	grail_seg = prune_(x=grail_seg) %>%
# 				bind_cols(cn = absolute_(rho=key_file$GRAIL_alpha[i],
# 										 psi=key_file$GRAIL_psi[i],
# 										 x=grail_seg$log2)) %>%
# 				mutate(n = cumsum(n))
#  	
#  	Chromosome = grail_seg[,"chrom"]
#  	Start = grail_seg[,"start"]
#  	End = grail_seg[,"end"]
#  	Calls = grail_seg[,"cn"]
#  	res = data.frame(Chromosome, Start, End, Calls)
#  	annot = read.csv(file="~/share/reference/IMPACT410_genes_for_copynumber.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
#  			select(hgnc_symbol, chr, start_position, end_position) %>%
#  			rename(Hugo_Symbol = hgnc_symbol,
#  				   Chromosome = chr,
#  				   Start = start_position,
#  				   End = end_position) %>%
#  			mutate(Chromosome = ifelse(Chromosome %in% "X", 23, Chromosome)) %>%
#  			arrange(as.numeric(Chromosome), as.numeric(Start))
#  				   
#  	annot_by_gene <- annot %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Hugo_Symbol = Hugo_Symbol)
#  	res_by_segment <- res %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Calls = Calls)
#  	fo <- findOverlaps(res_by_segment, annot_by_gene)
#  
#  	df <- data.frame(Hugo_Symbol=mcols(annot_by_gene)[subjectHits(fo),], Calls=mcols(res_by_segment)[queryHits(fo),])
#  	Hugo_Symbol = which(duplicated(df$Hugo_Symbol))
#  	for (j in 1:length(Hugo_Symbol)) {
#  		index = which(as.character(df$Hugo_Symbol)==as.character(df$Hugo_Symbol[Hugo_Symbol[j]]))
#  		df[index,2] = mean(df[index,2], na.rm=TRUE)
#  	}
#  	df = df %>% filter(!duplicated(Hugo_Symbol))
#  	df[,2] = round(df[,2])
#  	res = rep(0, nrow(annot))
#  	names(res) = annot[,"Hugo_Symbol"]
#  	res[as.character(df$Hugo_Symbol)] = df$Calls
#  	return(invisible(res))
# }
# g_bygene = do.call(cbind, g_bygene)
# colnames(g_bygene) = tracker$GRAIL_ID
# 
# index = c("CRLF2", "HLA-A", "HLA-B", "HLA-C", "AR", "HIST2H3D", "HIST2H3C", "HIST3H3",
# 		  "HIST1H3A", "HIST1H3B", "HIST1H3C", "HIST1H1C", "HIST1H2BD", "HIST1H3D",
# 		  "HIST1H3E", "HIST1H3F", "HIST1H3G", "HIST1H3H", "HIST1H3I", "HIST1H3J")
# i_bygene = i_bygene[!(rownames(i_bygene) %in% index),,drop=FALSE]
# g_bygene = g_bygene[!(rownames(g_bygene) %in% index),,drop=FALSE]
# 
# 
# annot = read.csv(file="~/share/reference/IMPACT410_genes_for_copynumber.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
#  			select(hgnc_symbol, chr, start_position, end_position) %>%
#  			rename(Hugo_Symbol = hgnc_symbol,
#  				   Chromosome = chr,
#  				   Start = start_position,
#  				   End = end_position) %>%
#  			mutate(Chromosome = ifelse(Chromosome %in% "X", 23, Chromosome)) %>%
#  			arrange(as.numeric(Chromosome), as.numeric(Start))
# annot = annot[!(annot$Hugo_Symbol %in% index),,drop=FALSE]
# 
# n = 23
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# ColSideColors = rep(NA, nrow(annot))
# ColSideColors = col_vector[as.numeric(annot[,"Chromosome"])]
# names(ColSideColors) = as.character(annot[,"Hugo_Symbol"])
# ColSideColors = ColSideColors[rownames(i_bygene)]
# RowSideColors = rep(NA, ncol(i_bygene))
# RowSideColors[grep("VB", colnames(i_bygene), fixed=TRUE)] = "salmon"
# RowSideColors[grep("VL", colnames(i_bygene), fixed=TRUE)] = "#FDAE61"
# RowSideColors[grep("VP", colnames(i_bygene), fixed=TRUE)] = "#ABDDA4"
# 
# index = is.na(ColSideColors)
# ColSideColors = ColSideColors[!index]
# i_bygene = i_bygene[!index,,drop=FALSE]
# g_bygene = g_bygene[!index,,drop=FALSE]
# 
# for (i in 1:nrow(tracker)) {
# 	psi = round(tracker$IMPACT_psi[i])
#  	psi = ifelse(psi==1, 2, psi)
#  	x = i_bygene[,i]
#  	x2 = rep(0, length(x))
#  	if (psi==2) {
# 	 	x2[x<(psi-1)] = -1
# 		x2[x>(psi+3)] = 1
# 	} else if (psi==3) {
# 		x2[x<(psi-2)] = -1
# 		x2[x>(psi+4)] = 1
# 	} else {
# 		x2[x<(psi-3)] = -1
# 		x2[x>(psi+5)] = 1
# 	}
# 	i_bygene[,i] = x2
# 	
# 	psi = round(tracker$GRAIL_psi[i])
#  	psi = ifelse(psi==1, 2, psi)
#  	x = g_bygene[,i]
#  	x2 = rep(0, length(x))
#  	if (psi==2) {
# 	 	x2[x<(psi-1)] = -1
# 		x2[x>(psi+3)] = 1
# 	} else if (psi==3) {
# 		x2[x<(psi-2)] = -1
# 		x2[x>(psi+4)] = 1
# 	} else {
# 		x2[x<(psi-3)] = -1
# 		x2[x>(psi+4)] = 1
# 	}
# 	g_bygene[,i] = x2
# }
# 
# pdf(file="../res/rebuttal/Heatmap_CN_tumor_abs_copy_all_genes.pdf", width=8)
# hmi = heatmap(t(i_bygene), Colv=NA, Rowv=as.dendrogram(hclust(dist(t(i_bygene)))),
# 			  col = c("steelblue", "white", "red"),
# 			  scale="none",
#  			  RowSideColors=RowSideColors,
#  			  ColSideColors=ColSideColors,
#  			  labRow=NA,
#  			  labCol=NA)
# dev.off()
# 
# pdf(file="../res/rebuttal/Heatmap_CN_cfDNA_abs_copy_all_genes.pdf", width=8)
# hmg = heatmap(t(g_bygene), Colv=NA, Rowv=as.dendrogram(hclust(dist(t(i_bygene)))),
# 			  col = c("steelblue", "white", "red"),
# 			  scale = "none",
#  			  RowSideColors=RowSideColors,
#  			  ColSideColors=ColSideColors,
#  			  labRow=NA,
#  			  labCol=NA)
# dev.off()
# 
# 
# 
# 
# #==================================================
# # % Agreement based on absolute copy numbers
# #==================================================
# if (FALSE) { pa = unlist(foreach (i=1:nrow(tracker)) %dopar% {
#  	cat(tracker$GRAIL_ID[i], "\n")
#  	x = res[[i]][,"cfDNA"] 
#  	y = res[[i]][,"Biopsy"]
#  	
#  	class = 1:max(c(x, y))
# 	z = matrix(0, ncol=length(class), nrow=length(class))
# 	
# 	for (ii in 1:length(class)) {
# 		for (jj in 1:length(class)) {
# 			z[ii,jj] = sum(x==class[ii] & y==class[jj])
# 		}
# 	}
# 	k = sum(diag(z))/sum(z)
# 	return(k)
# }) }
# 
# pa = data.frame(pa) %>%
# 	 mutate(Tissue = "")
# pa[grepl("VB", tracker$GRAIL_ID),"Tissue"] = "Breast"	 
# pa[grepl("VL", tracker$GRAIL_ID),"Tissue"] = "Lung"
# pa[grepl("VP", tracker$GRAIL_ID),"Tissue"] = "Prostate"
# 
# #==================================================
# # Updated ploidy
# #==================================================
# if (TRUE) { ploidy = foreach (i=1:nrow(key_file)) %dopar% {
# 	print(key_file$GRAIL_ID[i])
# 	impact_path = paste0("../res/rebuttal/MSK-IMPACT/facets/cncf/", key_file$TUMOR_ID[i], "_", key_file$NORMAL_ID[i], ".Rdata")
# 	impact_data = new.env()
# 	load(impact_path, envir=impact_data)
# 	impact_cn = impact_data$out2$jointseg %>%
# 			    select(chrom, pos = maploc, log2 = cnlr)
# 	impact_seg = impact_data$fit$cncf %>%
# 				 select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
# 	impact_seg = prune_(x=impact_seg) %>%
# 				 bind_cols(cn = absolute_(rho=key_file$IMPACT_alpha[i],
# 										  psi=key_file$IMPACT_psi[i],
# 										  x=impact_seg$log2)) %>%
# 				mutate(n = cumsum(n))
# 	impact_psi = (t(impact_seg[,"end"] - impact_seg[,"start"])%*%impact_seg[,"cn"])/sum(impact_seg[,"end"] - impact_seg[,"start"])
# 	grail_path = paste0("../res/rebuttal/GRAIL/facets/cncf/", key_file$GRAIL_ID[i], "_", key_file$GRAIL_ID[i], "-N.Rdata")
# 	grail_data = new.env()
# 	load(grail_path, envir=grail_data)
# 	grail_cn = grail_data$out2$jointseg %>%
# 			   select(chrom, pos = maploc, log2 = cnlr)
# 	grail_seg = grail_data$fit$cncf %>%
# 				select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
# 	
# 	fixed_cn = fix_6(grail_cn, grail_seg)
# 	grail_cn = fixed_cn[[1]]
# 	grail_seg = fixed_cn[[2]]
# 	
# 	grail_seg = prune_(x=grail_seg) %>%
# 				bind_cols(cn = absolute_(rho=key_file$GRAIL_alpha[i],
# 										 psi=key_file$GRAIL_psi[i],
# 										 x=grail_seg$log2)) %>%
# 				mutate(n = cumsum(n))
# 	grail_psi = (t(grail_seg[,"end"] - grail_seg[,"start"])%*%grail_seg[,"cn"])/sum(grail_seg[,"end"] - grail_seg[,"start"])
# 	return(c(grail_psi, impact_psi))
# } }
# ploidy = do.call(rbind, ploidy)
# 
# tmp = cbind(tracker, pa)
# tmp[,"GRAIL_psi"] = ploidy[,1]
# tmp[,"IMPACT_psi"] = ploidy[,2]
# 
# tmp.0 = tmp %>%
# 		mutate(GRAIL_alpha = ifelse(GRAIL_alpha==1, 0, GRAIL_alpha)) %>%
# 		mutate(IMPACT_alpha = ifelse(IMPACT_alpha==1, 0, IMPACT_alpha)) %>%
# 		mutate(Cat = "Both estimate available") %>%
# 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha == 0, "No estimate in both", Cat)) %>%
# 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha != 0, "No estimate in cfDNA", Cat)) %>%
# 		mutate(Cat = ifelse(GRAIL_alpha != 0 & IMPACT_alpha == 0, "No estimate in Biopsy", Cat)) %>%
# 		mutate(Purity = "Purity")
# 		
# plot.0 = ggplot(tmp.0, aes(y = GRAIL_alpha, x = IMPACT_alpha, shape = Tissue, fill = Cat)) +
# 		 geom_abline(slope = 1, color = "goldenrod3", linetype = 1) +
# 		 geom_point(alpha = .8, size = 2.5) +
# 		 scale_fill_manual(values = c("Both estimate available"="salmon", "No estimate in both"="#FDAE61", "No estimate in cfDNA"="#ABDDA4", "No estimate in Biopsy"="steelblue")) +
# 		 scale_shape_manual(values = c(24, 21, 22)) +
# 		 theme_bw(base_size=15) +
# 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.26, 0.7), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# 		 labs(y="\ncfDNA\n", x="\nBiopsy\n") +
# 		 coord_cartesian(xlim=c(0,1), ylim = c(0, 1)) +
# 		 facet_wrap(~Purity) +
# 		 guides(shape=guide_legend(title=c("Tissue"), override.aes=list(fill="black"))) +
# 		 guides(fill=guide_legend(title=c("Purity")))
# 		 
# 		 
# pdf(file="../res/rebuttal/Comparison_Purity.pdf", width=6, height=6)
# print(plot.0)
# dev.off()
# 
# tmp.0 = tmp %>%
# 		mutate(GRAIL_alpha = ifelse(GRAIL_alpha==1, 0, GRAIL_alpha)) %>%
# 		mutate(IMPACT_alpha = ifelse(IMPACT_alpha==1, 0, IMPACT_alpha)) %>%
# 		mutate(Cat = "Both estimate available") %>%
# 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha == 0, "No estimate in both", Cat)) %>%
# 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha != 0, "No estimate in cfDNA", Cat)) %>%
# 		mutate(Cat = ifelse(GRAIL_alpha != 0 & IMPACT_alpha == 0, "No estimate in Biopsy", Cat)) %>%
# 		mutate(Ploidy = "Ploidy")
# 		
# plot.0 = ggplot(tmp.0, aes(y = GRAIL_psi, x = IMPACT_psi, shape = Tissue, fill = Cat)) +
# 		 geom_abline(slope = 1, color = "goldenrod3", linetype = 1) +
# 		 geom_point(alpha = .8, size = 2.5) +
#  		 scale_fill_manual(values = c("Both estimate available"="salmon", "No estimate in both"="#FDAE61", "No estimate in cfDNA"="#ABDDA4", "No estimate in Biopsy"="steelblue")) +
# 		 scale_shape_manual(values = c(24, 21, 22)) +
# 		 theme_bw(base_size=15) +
# 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.26, 0.7), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# 		 labs(y="\ncfDNA\n", x="\nBiopsy\n") +
# 		 coord_cartesian(xlim=c(1.5,4), ylim = c(1.5,4)) +
# 		 facet_wrap(~Ploidy) +
# 		 guides(shape=guide_legend(title=c("Tissue"), override.aes=list(fill="black"))) +
# 		 guides(fill=guide_legend(title=c("Purity")))
# 		 
# pdf(file="../res/rebuttal/Comparison_Ploidy.pdf", width=6, height=6)
# print(plot.0)
# dev.off()
# 
# 
# tmp.0 = tmp %>%
# 		mutate(Tissue = "All samples")
# 
# plot.0 = ggplot(tmp.0, aes(x = 100*pa)) +
# 		 geom_histogram(color="black", fill="#2B83BA") +
# 		 coord_cartesian(xlim=c(0,100)) +
# 		 theme_bw(base_size=15) +
# 		 facet_wrap(~Tissue) +
# 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# 		 labs(y="\nFrequency\n", x="\n Agreement (%)\n")
# 		 
# pdf(file="../res/rebuttal/Percent_Agreement_Combined.pdf", width=6.5)
# print(plot.0)
# dev.off()
# 
# plot.0 = ggplot(tmp, aes(x = 100*pa)) +
# 		 geom_histogram(color="black", fill="#2B83BA") +
# 		 coord_cartesian(xlim=c(0,100)) +
# 		 theme_bw(base_size=15) +
# 		 facet_wrap(~Tissue) +
# 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# 		 labs(y="\nFrequency\n", x="\n Agreement (%)\n")
# 		 
# pdf(file="../res/rebuttal/Percent_Agreement_Tissue.pdf", width=15)
# print(plot.0)
# dev.off()
# 
# tmp.0 = tmp %>%
# 		mutate(GRAIL_alpha = ifelse(GRAIL_alpha==1, 0, GRAIL_alpha)) %>%
# 		mutate(IMPACT_alpha = ifelse(IMPACT_alpha==1, 0, IMPACT_alpha)) %>%
# 		mutate(Cat = "Both estimate available") %>%
# 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha == 0, "No estimate in both", Cat)) %>%
# 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha != 0, "No estimate in cfDNA", Cat)) %>%
# 		mutate(Cat = ifelse(GRAIL_alpha != 0 & IMPACT_alpha == 0, "No estimate in Biopsy", Cat)) %>%
# 		mutate(Ploidy = "Agreement versus cfDNA ploidy")
# 		
# 
# plot.0 = ggplot(tmp.0, aes(y = GRAIL_psi, x = 100*pa, shape = Tissue, fill = Cat)) +
# 		 geom_point(alpha = .8, size = 2.5) +
# 		 scale_fill_manual(values = c("Both estimate available"="salmon", "No estimate in both"="#FDAE61", "No estimate in cfDNA"="#ABDDA4", "No estimate in Biopsy"="steelblue")) +
# 		 scale_shape_manual(values = c(24, 21, 22)) +
# 		 theme_bw(base_size=15) +
# 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.26, 0.7), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# 		 labs(y="\nPloidy in cfDNA\n", x="\nAgreement (%)\n") +
# 		 coord_cartesian(xlim=c(0,100), ylim = c(1.5,4)) +
# 		 facet_wrap(~Ploidy) +
# 		 guides(shape=guide_legend(title=c("Tissue"), override.aes=list(fill="black"))) +
# 		 guides(fill=guide_legend(title=c("Purity")))
# 		 
# pdf(file="../res/rebuttal/Agreement_Ploidy_cfDNA.pdf", width=6, height=6)
# print(plot.0)
# dev.off()
# 
# tmp.0 = tmp %>%
# 		mutate(GRAIL_alpha = ifelse(GRAIL_alpha==1, 0, GRAIL_alpha)) %>%
# 		mutate(IMPACT_alpha = ifelse(IMPACT_alpha==1, 0, IMPACT_alpha)) %>%
# 		mutate(Cat = "Both estimate available") %>%
# 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha == 0, "No estimate in both", Cat)) %>%
# 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha != 0, "No estimate in cfDNA", Cat)) %>%
# 		mutate(Cat = ifelse(GRAIL_alpha != 0 & IMPACT_alpha == 0, "No estimate in Biopsy", Cat)) %>%
# 		mutate(Ploidy = "Agreement versus tumor ploidy")
# 		
# 
# plot.0 = ggplot(tmp.0, aes(y = IMPACT_psi, x = 100*pa, shape = Tissue, fill = Cat)) +
# 		 geom_point(alpha = .8, size = 2.5) +
# 		 scale_fill_manual(values = c("Both estimate available"="salmon", "No estimate in both"="#FDAE61", "No estimate in cfDNA"="#ABDDA4", "No estimate in Biopsy"="steelblue")) +
# 		 scale_shape_manual(values = c(24, 21, 22)) +
# 		 theme_bw(base_size=15) +
# 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.8, 0.7), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# 		 labs(y="\nPloidy in tumor\n", x="\nAgreement (%)\n") +
# 		 coord_cartesian(xlim=c(0,100), ylim = c(1.5,4)) +
# 		 facet_wrap(~Ploidy) +
# 		 guides(shape=guide_legend(title=c("Tissue"), override.aes=list(fill="black"))) +
# 		 guides(fill=guide_legend(title=c("Purity")))
# 		 
# pdf(file="../res/rebuttal/Agreement_Ploidy_biopsy.pdf", width=6, height=6)
# print(plot.0)
# dev.off()
# 
# tmp.0 = tmp %>%
# 		mutate(GRAIL_alpha = ifelse(GRAIL_alpha==1, 0, GRAIL_alpha)) %>%
# 		mutate(IMPACT_alpha = ifelse(IMPACT_alpha==1, 0, IMPACT_alpha)) %>%
# 		mutate(Cat = "Both estimate available") %>%
# 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha == 0, "No estimate in both", Cat)) %>%
# 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha != 0, "No estimate in cfDNA", Cat)) %>%
# 		mutate(Cat = ifelse(GRAIL_alpha != 0 & IMPACT_alpha == 0, "No estimate in Biopsy", Cat)) %>%
# 		mutate(Ploidy = "Agreement versus ploidy difference")
# 		
# plot.0 = ggplot(tmp.0, aes(y = abs(IMPACT_psi - GRAIL_psi), x = 100*pa, shape = Tissue, fill = Cat)) +
# 		 geom_point(alpha = .8, size = 2.5) +
# 		 scale_fill_manual(values = c("Both estimate available"="salmon", "No estimate in both"="#FDAE61", "No estimate in cfDNA"="#ABDDA4", "No estimate in Biopsy"="steelblue")) +
# 		 scale_shape_manual(values = c(24, 21, 22)) +
# 		 theme_bw(base_size=15) +
# 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.8, 0.7), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# 		 labs(y="\n| cfDNA - Tumor |\n", x="\nAgreement (%)\n") +
# 		 coord_cartesian(xlim=c(0,100), ylim = c(0,2)) +
# 		 facet_wrap(~Ploidy) +
# 		 guides(shape=guide_legend(title=c("Tissue"), override.aes=list(fill="black"))) +
# 		 guides(fill=guide_legend(title=c("Purity")))
# 		 
# pdf(file="../res/rebuttal/Agreement_Ploidy_biopsy_cfDNA.pdf", width=6, height=6)
# print(plot.0)
# dev.off()
# 
# tmp.0 = tmp %>%
# 		mutate(GRAIL_alpha = ifelse(GRAIL_alpha==1, 0, GRAIL_alpha)) %>%
# 		mutate(IMPACT_alpha = ifelse(IMPACT_alpha==1, 0, IMPACT_alpha)) %>%
# 		mutate(Cat = "Both estimate available") %>%
# 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha == 0, "No estimate in both", Cat)) %>%
# 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha != 0, "No estimate in cfDNA", Cat)) %>%
# 		mutate(Cat = ifelse(GRAIL_alpha != 0 & IMPACT_alpha == 0, "No estimate in Biopsy", Cat)) %>%
# 		mutate(Ploidy = "Agreement versus purity difference")
# 		
# 
# plot.0 = ggplot(tmp.0, aes(y = abs(IMPACT_alpha - GRAIL_alpha), x = 100*pa, shape = Tissue, fill = Cat)) +
# 		 geom_point(alpha = .8, size = 2.5) +
# 		 scale_fill_manual(values = c("Both estimate available"="salmon", "No estimate in both"="#FDAE61", "No estimate in cfDNA"="#ABDDA4", "No estimate in Biopsy"="steelblue")) +
# 		 scale_shape_manual(values = c(24, 21, 22)) +
# 		 theme_bw(base_size=15) +
# 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.8, 0.7), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# 		 labs(y="\n| cfDNA - Tumor |\n", x="\nAgreement (%)\n") +
# 		 coord_cartesian(xlim=c(0,100), ylim = c(0,1)) +
# 		 facet_wrap(~Ploidy) +
# 		 guides(shape=guide_legend(title=c("Tissue"), override.aes=list(fill="black"))) +
# 		 guides(fill=guide_legend(title=c("Purity")))
# 		 
# pdf(file="../res/rebuttal/Agreement_Purity_biopsy_cfDNA.pdf", width=6, height=6)
# print(plot.0)
# dev.off()


# for (i in 1:nrow(tracker)) {
# 	x = g_cn_jrf[,i]
# 	y = ifelse(i_cn[,i]==1, 1, 0)
# 	z = rep(0, length(x))
# 	if (length(unique(y))==2) {
# 		pred = prediction(x,y)
# 		perf = performance(pred, measure = "acc")
# 		ind = which.max( slot(perf, "y.values")[[1]] )
# 		co = slot(perf, "x.values")[[1]][ind]
# 	} else {
# 		co = Inf
# 	}
# 	z[x>co] = 1
# 	
# 	y = ifelse(i_cn[,i]==-1, 0, 1)
#     if (length(unique(y))==2) {
#     	pred = prediction(x,y)
#     	perf = performance(pred, measure = "acc")
#         ind = which.max( slot(perf, "y.values")[[1]] )
#         co = slot(perf, "x.values")[[1]][ind]
#     } else {
#     	co = -Inf
#     }
#     z[x<co] = -1
#     g_cn_jrf[,i] = z
# }
# 
# i_cn = read.csv(file="~/share/data/common/cbioportal_repos/msk-impact/msk-impact/data_CNA.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
#            filter(Hugo_Symbol %in% rownames(g_cn))
# colnames(i_cn) = gsub(pattern=".", replacement="-", x=colnames(i_cn), fixed=TRUE)
# rownames(i_cn) = i_cn$Hugo_Symbol
# i_cn = i_cn[,tracker$DMP_ID,drop=FALSE]
# colnames(i_cn) = tracker$GRAIL_ID
# i_cn = i_cn[rownames(g_cn),,drop=FALSE]
# i_cn[i_cn>0] = 1
# i_cn[i_cn<0] = -1
# i_cn[is.na(i_cn)] = 0
# i_cn = as.matrix(i_cn)
# 
# for (i in 1:nrow(tracker)) {
# 	x = g_cn_cmo[,i]
# 	y = ifelse(i_cn[,i]==1, 1, 0)
# 	z = rep(0, length(x))
# 	if (length(unique(y))==2) {
# 		pred = prediction(x,y)
# 		perf = performance(pred, measure = "acc")
# 		ind = which.max( slot(perf, "y.values")[[1]] )
# 		co = slot(perf, "x.values")[[1]][ind]
# 	} else {
# 		co = Inf
# 	}
# 	z[x>co] = 1
# 	
# 	y = ifelse(i_cn[,i]==-1, 0, 1)
#     if (length(unique(y))==2) {
#     	pred = prediction(x,y)
#     	perf = performance(pred, measure = "acc")
#         ind = which.max( slot(perf, "y.values")[[1]] )
#         co = slot(perf, "x.values")[[1]][ind]
#     } else {
#     	co = -Inf
#     }
#     z[x<co] = -1
#     g_cn_cmo[,i] = z
# }
# 
# g_cn = g_cn_jrf
# g_cn[g_cn_cmo==-1] = -1
# g_cn[g_cn_cmo==1] = 1
# 
# index = which.max(apply(g_cn, 2, function(x) { sum(x==-1) }))
# g_cn[g_cn[,index]==-1,index] = 0
# 
# i_cn = read.csv(file="~/share/data/common/cbioportal_repos/msk-impact/msk-impact/data_CNA.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
#            filter(Hugo_Symbol %in% rownames(g_cn))
# colnames(i_cn) = gsub(pattern=".", replacement="-", x=colnames(i_cn), fixed=TRUE)
# rownames(i_cn) = i_cn$Hugo_Symbol
# i_cn = i_cn[,tracker$DMP_ID,drop=FALSE]
# colnames(i_cn) = tracker$GRAIL_ID
# i_cn = i_cn[rownames(g_cn),,drop=FALSE]
# i_cn[i_cn>0] = 1
# i_cn[i_cn<0] = -1
# i_cn[is.na(i_cn)] = 0
# i_cn = as.matrix(i_cn)
# 
# annot = read.csv(file="~/share/reference/IMPACT410_genes_for_copynumber.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
# 				 select(hgnc_symbol, chr, start_position, end_position) %>%
# 				 rename(Hugo_Symbol = hgnc_symbol,
# 				   		Chromosome = chr,
# 				   		Start = start_position,
# 				   		End = end_position) %>%
# 				 mutate(Chromosome = ifelse(Chromosome %in% "X", 23, Chromosome))
# n = 23
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# ColSideColors = rep(NA, nrow(annot))
# ColSideColors = col_vector[as.numeric(annot[,"Chromosome"])]
# names(ColSideColors) = as.character(annot[,"Hugo_Symbol"])
# ColSideColors = ColSideColors[rownames(i_cn)]
# RowSideColors = rep(NA, ncol(i_cn))
# RowSideColors[grep("VB", colnames(i_cn), fixed=TRUE)] = "salmon"
# RowSideColors[grep("VL", colnames(i_cn), fixed=TRUE)] = "#FDAE61"
# RowSideColors[grep("VP", colnames(i_cn), fixed=TRUE)] = "#ABDDA4"
# 
# index = is.na(ColSideColors)
# ColSideColors = ColSideColors[!index]
# i_cn = i_cn[!index,,drop=FALSE]
# g_cn = g_cn[!index,,drop=FALSE]
# 
# pdf(file="../res/rebuttal/Heatmap_CN_IMPACT_All_Genes.pdf", width=8)
# hmi = heatmap(t(i_cn), Colv=NA, Rowv=as.dendrogram(hclust(dist(t(i_cn)))), col=c("blue", "grey99", "red"), scale="none",
# 			  RowSideColors=RowSideColors,
# 			  ColSideColors=ColSideColors,
# 			  labRow=NA,
# 			  labCol=NA)
# dev.off()
# 
# pdf(file="../res/rebuttal/Heatmap_CN_GRAIL_All_Genes.pdf", width=8)
# hmg = heatmap(t(g_cn), Colv=NA, Rowv=as.dendrogram(hclust(dist(t(i_cn)))), col=c("blue", "grey99", "red"), scale="none",
# 			  RowSideColors=RowSideColors,
# 			  ColSideColors=ColSideColors,
# 			  labRow=NA,
# 			  labCol=NA)
# dev.off()
# 
# index = apply(i_cn, 1, function(x) {sum(x==0)==length(x)})
# ColSideColors = ColSideColors[!index]
# i_cn = i_cn[!index,,drop=FALSE]
# g_cn = g_cn[!index,,drop=FALSE]
# 
# pdf(file="../res/rebuttal/Heatmap_CN_IMPACT_HD_and_Amp_Only.pdf", width=8)
# hmi = heatmap(t(i_cn), Colv=NA, Rowv=as.dendrogram(hclust(dist(t(i_cn)))), col=c("blue", "grey99", "red"), scale="none",
# 			  RowSideColors=RowSideColors,
# 			  ColSideColors=ColSideColors,
# 			  labRow=NA,
# 			  labCol=NA)
# dev.off()
# 
# pdf(file="../res/rebuttal/Heatmap_CN_GRAIL_HD_and_Amp_Only.pdf", width=8)
# hmg = heatmap(t(g_cn), Colv=NA, Rowv=as.dendrogram(hclust(dist(t(i_cn)))), col=c("blue", "grey99", "red"), scale="none",
# 			  RowSideColors=RowSideColors,
# 			  ColSideColors=ColSideColors,
# 			  labRow=NA,
# 			  labCol=NA)
# dev.off()
