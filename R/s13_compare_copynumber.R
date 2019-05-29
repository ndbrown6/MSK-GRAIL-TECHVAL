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
		 dplyr::select(chrom, pos = maploc, log2 = cnlr)
	seg = y$cncf %>%
		  dplyr::select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
	if (fix_6) {
		fixed_cn = fix_6(cn, seg, ix)
		cn = fixed_cn[[1]]
		seg = fixed_cn[[2]]
	}
	seg = prune_(x=seg, n) %>%
		  mutate(n = cumsum(n))
	
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
	abline(v=1, col=transparent_rgb("goldenrod3", 255), lty=3, lwd=.5)
	for (j in 1:23) {
		abline(v=CytoBand[j,"end"], col=transparent_rgb("goldenrod3", 255), lty=3, lwd=.5)
	}
	axis(1, at = .5*(CytoBand[,"start"]+CytoBand[,"end"]), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
    for (k in c(1, 2, 4, 8, 14)) {
		abline(h=(log2(((purity)*k + (1-purity)*2)/((purity)*ploidy + (1-purity)*2))), col=transparent_rgb("brown", 155), lty=3)
	}
	rect(xleft=1-1e10, xright=CytoBand[23,"end"]+1e10, ybottom=4, ytop=6, col="lightgrey", border="black", lwd=1.5)
	title(main = title, line=-1, cex.main=.75, font.main=1)
    box(lwd=1.5)
}

'plot_log3_' <- function(x, y, title = "")
{
   	par(mar=c(5, 5, 4, 2)+.1)
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
	plot(x[,"Position"], x[,"Log2Ratio"], type="p", pch=".", cex=1, col="grey75", axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,5))
	for (j in 1:nrow(y)) {
 		lines(x=c(y[j,"Start"], y[j,"End"]), y=rep(y[j,"Log2Ratio"],2), lty=1, lwd=2.5, col="red")
 	}
  	axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1, las = 1)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	abline(v=1, col="goldenrod3", lty=3, lwd=.5)
	abline(h=0, col="red", lty=1, lwd=1)
	for (j in 2:23) {
		v = start[j]
		abline(v=v, col="goldenrod3", lty=3, lwd=.5)
	}
	abline(v=max(x[,"Position"]), col="goldenrod3", lty=3, lwd=.5)
	axis(1, at = .5*(start+end), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)	
	rect(xleft=1-1e10, xright=x[nrow(x),"Position"]+1e10, ybottom=4, ytop=6, col="lightgrey", border="black", lwd=1.5)
	title(main = title, line=-1, cex.main=.75, font.main=1)
    box(lwd=1.5)
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

#==================================================
# update msk_impact tumor (alpha, psi)
#==================================================
key_file = read_tsv(file=url_master_key, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   dplyr::select(PATIENT_ID, GRAIL_ID, DMP_ID, TUMOR_ID, NORMAL_ID, GRAIL_alpha, GRAIL_psi, IMPACT_alpha, IMPACT_psi)

if (TRUE) { foreach (i=1:nrow(key_file)) %dopar% {
	print(key_file$GRAIL_ID[i])
	impact_path = paste0("../res/rebuttal/msk_impact/facets/cncf/", key_file$TUMOR_ID[i], "_", key_file$NORMAL_ID[i], ".Rdata")
	impact_data = new.env()
	load(impact_path, envir=impact_data)
	
	impact_cn = impact_data$out2$jointseg %>%
			    dplyr::select(chrom, pos = maploc, log2 = cnlr)
	impact_seg = impact_data$fit$cncf %>%
				 dplyr::select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
	impact_seg = prune_(x=impact_seg) %>%
				 bind_cols(cn = absolute_(rho=key_file$IMPACT_alpha[i],
										  psi=key_file$IMPACT_psi[i],
										  x=impact_seg$log2)) %>%
				mutate(n = cumsum(n))
				
	file_path = paste0("../res/rebuttal/msk_impact/facets/plots/ext/", key_file$GRAIL_ID[i], ".pdf")
										
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
	abline(v=1, col=transparent_rgb("goldenrod3", 255), lty=3, lwd=.5)
	abline(h=0, col="black")
	for (j in 1:23) {
		v = max(which(impact_seg[,"chrom"]==j))
		abline(v=impact_seg[v,"n"], col=transparent_rgb("goldenrod3", 255), lty=3, lwd=.5)
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
		abline(h=(log2(((purity)*k + (1-purity)*2)/((purity)*ploidy + (1-purity)*2))), col=transparent_rgb("brown", 155), lty=3)
	}
    box(lwd=1.5)
	dev.off()
} }

#==================================================
# update grail cfDNA uncollapsed (alpha, psi)
#==================================================
key_file = read_tsv(file=url_master_key, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   dplyr::select(PATIENT_ID, GRAIL_ID, DMP_ID, TUMOR_ID, NORMAL_ID, GRAIL_alpha, GRAIL_psi, IMPACT_alpha, IMPACT_psi)

if (TRUE) { foreach (i=1:nrow(key_file)) %dopar% {
 	print(key_file$GRAIL_ID[i])
 	grail_path = paste0("../res/rebuttal/uncollapsed_bam/facets/cncf/", key_file$GRAIL_ID[i], "-T_", key_file$GRAIL_ID[i], "-N.Rdata")
 	grail_data = new.env()
 	load(grail_path, envir=grail_data)
 	
 	grail_cn = grail_data$out2$jointseg %>%
 			   dplyr::select(chrom, pos = maploc, log2 = cnlr)
 	grail_seg = grail_data$fit$cncf %>%
 				dplyr::select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
 	
 	fixed_cn = fix_6(grail_cn, grail_seg)
 	grail_cn = fixed_cn[[1]]
 	grail_seg = fixed_cn[[2]]
 	
 	grail_seg = prune_(x=grail_seg) %>%
 				bind_cols(cn = absolute_(rho=key_file$GRAIL_alpha[i],
 										 psi=key_file$GRAIL_psi[i],
 										 x=grail_seg$log2)) %>%
 				mutate(n = cumsum(n))
 	
 				
 	file_path = paste0("../res/rebuttal/uncollapsed_bam/facets/plots/ext/", key_file$GRAIL_ID[i], ".pdf")
 										
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
 	abline(v=1, col=transparent_rgb("goldenrod3", 255), lty=3, lwd=.5)
 	abline(h=0, col="black")
 	for (j in 1:23) {
 		v = max(which(grail_seg[,"chrom"]==j))
 		abline(v=grail_seg[v,"n"], col=transparent_rgb("goldenrod3", 255), lty=3, lwd=.5)
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
 		abline(h=(log2(((purity)*k + (1-purity)*2)/((purity)*ploidy + (1-purity)*2))), col=transparent_rgb("brown", 155), lty=3)
 	}
    box(lwd=1.5)
 	dev.off()
} }

#==================================================
# update grail cfDNA collapsed (alpha, psi)
#==================================================
if (TRUE) { foreach (i=1:nrow(key_file)) %dopar% {
 	print(key_file$GRAIL_ID[i])
 	grail_path = paste0("../res/rebuttal/collapsed_bam/facets/cncf/", key_file$GRAIL_ID[i], "_", key_file$GRAIL_ID[i], "-N.Rdata")
 	grail_data = new.env()
 	load(grail_path, envir=grail_data)
 	
 	grail_cn = grail_data$out2$jointseg %>%
 			   dplyr::select(chrom, pos = maploc, log2 = cnlr)
 	grail_seg = grail_data$fit$cncf %>%
 				dplyr::select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
 	
 	fixed_cn = fix_6(grail_cn, grail_seg)
 	grail_cn = fixed_cn[[1]]
 	grail_seg = fixed_cn[[2]]
 	
 	grail_seg = prune_(x=grail_seg) %>%
 				bind_cols(cn = absolute_(rho=key_file$GRAIL_alpha[i],
 										 psi=key_file$GRAIL_psi[i],
 										 x=grail_seg$log2)) %>%
 				mutate(n = cumsum(n))
 	
 				
 	file_path = paste0("../res/rebuttal/collapsed_bam/facets/plots/ext/", key_file$GRAIL_ID[i], ".pdf")
 										
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
 	abline(v=1, col=transparent_rgb("goldenrod3", 255), lty=3, lwd=.5)
 	abline(h=0, col="black")
 	for (j in 1:23) {
 		v = max(which(grail_seg[,"chrom"]==j))
 		abline(v=grail_seg[v,"n"], col=transparent_rgb("goldenrod3", 255), lty=3, lwd=.5)
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
 		abline(h=(log2(((purity)*k + (1-purity)*2)/((purity)*ploidy + (1-purity)*2))), col=transparent_rgb("brown", 155), lty=3)
 	}
    box(lwd=1.5)
 	dev.off()
} }
 
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
				   mutate(ctdna_fraction_cat_w_ne = factor(ctdna_fraction_cat_w_ne))
				   
d = density(tmp.0$correlation_coefficient, from=-.5, to=1)
d$x = c(d$x, rev(d$x))
d$y = c(d$y, rep(0, length(d$y)))
d = data_frame(x=d$x, y=d$y, facets=paste0("All samples (n = ", nrow(tmp.0), ")"))

pdf(file="../res/rebuttal/distribution_correlation_coefficient_all_log2_ratio.pdf", width=5, height=6)
plot(d$x, d$y, type="n", xlab="", ylab="", frame.plot=FALSE, axes=FALSE, ylim=c(0,1.6))
polygon(x=d$x, y=d$y, col="grey80", border="grey50", lwd=1.5)
abline(v=0, lty=3, lwd=1, col="grey10")
axis(1, at = NULL, cex.axis = 1.15, las = 1)
axis(2, at = NULL, cex.axis = 1.15, las = 1)
mtext(side = 1, text = "Correlation coefficient", line = 3.15, cex = 1.25)
mtext(side = 2, text = "Density", line = 3.15, cex = 1.25)
rect(-10, 1.5, 10, 1.75, col="lightgrey", border="black", lwd=1.5)
title(main=d[1,3], line=-1.35, cex.main=1, font.main=1)
box(lwd=1.5)
dev.off()

d = density(tmp.0$correlation_coefficient[tmp.0$tissue=="Breast"], from=-.5, to=1)
d$x = c(d$x, rev(d$x))
d$y = c(d$y, rep(0, length(d$y)))
d = data_frame(x=d$x, y=d$y, facets=paste0("Breast (n = ", sum(tmp.0$tissue=="Breast"), ")"))

pdf(file="../res/rebuttal/distribution_correlation_coefficient_breast_log2_ratio.pdf", width=5, height=6)
plot(d$x, d$y, type="n", xlab="", ylab="", frame.plot=FALSE, axes=FALSE, ylim=c(0,1.6))
polygon(x=d$x, y=d$y, col=transparent_rgb("salmon", 205), border="salmon", lwd=1.5)
abline(v=0, lty=3, lwd=1, col="grey10")
axis(1, at = NULL, cex.axis = 1.15, las = 1)
axis(2, at = NULL, cex.axis = 1.15, las = 1)
mtext(side = 1, text = "Correlation coefficient", line = 3.15, cex = 1.25)
mtext(side = 2, text = "Density", line = 3.15, cex = 1.25)
rect(-10, 1.5, 10, 1.75, col="lightgrey", border="black", lwd=1.5)
title(main=d[1,3], line=-1.35, cex.main=1, font.main=1)
box(lwd=1.5)
dev.off()

d = density(tmp.0$correlation_coefficient[tmp.0$tissue=="Lung"], from=-.6, to=1)
d$x = c(d$x, rev(d$x))
d$y = c(d$y, rep(0, length(d$y)))
d = data_frame(x=d$x, y=d$y, facets=paste0("Lung (n = ", sum(tmp.0$tissue=="Lung"), ")"))

pdf(file="../res/rebuttal/distribution_correlation_coefficient_lung_log2_ratio.pdf", width=5, height=6)
plot(d$x, d$y, type="n", xlab="", ylab="", frame.plot=FALSE, axes=FALSE, ylim=c(0,1.6))
polygon(x=d$x, y=d$y, col=transparent_rgb("#FDAE61", 205), border="#FDAE61", lwd=1.5)
abline(v=0, lty=3, lwd=1, col="grey10")
axis(1, at = NULL, cex.axis = 1.15, las = 1)
axis(2, at = NULL, cex.axis = 1.15, las = 1)
mtext(side = 1, text = "Correlation coefficient", line = 3.15, cex = 1.25)
mtext(side = 2, text = "Density", line = 3.15, cex = 1.25)
rect(-10, 1.5, 10, 1.75, col="lightgrey", border="black", lwd=1.5)
title(main=d[1,3], line=-1.35, cex.main=1, font.main=1)
box(lwd=1.5)
dev.off()

d = density(tmp.0$correlation_coefficient[tmp.0$tissue=="Prostate"], from=-.6, to=1)
d$x = c(d$x, rev(d$x))
d$y = c(d$y, rep(0, length(d$y)))
d = data_frame(x=d$x, y=d$y, facets=paste0("Prostate (n = ", sum(tmp.0$tissue=="Prostate"), ")"))

pdf(file="../res/rebuttal/distribution_correlation_coefficient_prostate_log2_ratio.pdf", width=5, height=6)
plot(d$x, d$y, type="n", xlab="", ylab="", frame.plot=FALSE, axes=FALSE, ylim=c(0,1.6))
polygon(x=d$x, y=d$y, col=transparent_rgb("#ABDDA4", 205), border="#ABDDA4", lwd=1.5)
abline(v=0, lty=3, lwd=1, col="grey10")
axis(1, at = NULL, cex.axis = 1.15, las = 1)
axis(2, at = NULL, cex.axis = 1.15, las = 1)
mtext(side = 1, text = "Correlation coefficient", line = 3.15, cex = 1.25)
mtext(side = 2, text = "Density", line = 3.15, cex = 1.25)
rect(-10, 1.5, 10, 1.75, col="lightgrey", border="black", lwd=1.5)
title(main=d[1,3], line=-1.35, cex.main=1, font.main=1)
box(lwd=1.5)
dev.off()

colours = c("Breast"="salmon", "Lung"="#FDAE61", "Prostate"="#ABDDA4")
shapes = c("Yes"=17, "No"=19)
pdf(file="../res/rebuttal/ctDNA_Fraction_versus_Correlation_Coefficient_Log2.pdf", width=5, height=6)
plot(tmp.0$ctdna_fraction_w_ne*100, tmp.0$correlation_coefficient, type="p", xlab="", ylab="", col=unlist(lapply(colours[tmp.0$tissue], transparent_rgb, alpha=205)), pch=shapes[tmp.0$ctdna_fraction_cat_w_ne], cex=1.15, frame.plot=FALSE, axes=FALSE, ylim=c(-.5,1.09375))
axis(1, at = NULL, cex.axis = 1.15, las = 1)
axis(2, at = NULL, cex.axis = 1.15, las = 1)
mtext(side = 1, text = "ctDNA fraction (%)", line = 3.15, cex = 1.25)
mtext(side = 2, text = "Correlation coefficient", line = 3.15, cex = 1.25)
rect(-100, 1, 100, 1.75, col="lightgrey", border="black", lwd=1.5)
title(main="All samples", line=-1.35, cex.main=1, font.main=1)
fit = lm(correlation_coefficient ~ log(ctdna_fraction_w_ne), data = tmp.0 %>% filter(ctdna_fraction_cat!="NE") %>% filter(tissue=="Breast"))
x = seq(from=0, to=.9, length=1000)
y = predict(fit, newdata = list(ctdna_fraction_w_ne=x))
points(x*100, y, type="l", col=colours[1], lwd=3, lty=1)
fit = lm(correlation_coefficient ~ log(ctdna_fraction_w_ne), data = tmp.0 %>% filter(ctdna_fraction_cat!="NE") %>% filter(tissue=="Lung"))
x = seq(from=0, to=.9, length=1000)
y = predict(fit, newdata = list(ctdna_fraction_w_ne=x))
points(x*100, y, type="l", col=colours[2], lwd=3, lty=1)
fit = lm(correlation_coefficient ~ log(ctdna_fraction_w_ne), data = tmp.0 %>% filter(ctdna_fraction_cat!="NE") %>% filter(tissue=="Prostate"))
x = seq(from=0, to=.9, length=1000)
y = predict(fit, newdata = list(ctdna_fraction_w_ne=x))
points(x*100, y, type="l", col=colours[3], lwd=3, lty=1)
box(lwd=1.5)
legend(x=60, y=.1, legend=c("No", "Yes"), pch=shapes, title="ctDNA fraction", box.lwd=-1, cex=.85)
legend(x=60, y=-.2, legend=names(colours), pch=19, col=colours, title="Cancer type", box.lwd=-1, cex=.85)
dev.off()

pdf(file="../res/rebuttal/Tumor_Purity_versus_Correlation_Coefficient_Log2.pdf", width=5, height=6)
plot(tmp.0$tumor_purity*100, tmp.0$correlation_coefficient, type="p", xlab="", ylab="", col=unlist(lapply(colours[tmp.0$tissue], transparent_rgb, alpha=205)), pch=19, cex=1.15, frame.plot=FALSE, axes=FALSE, ylim=c(-.5,1.09375), xlim=c(10,100))
axis(1, at = NULL, cex.axis = 1.15, las = 1)
axis(2, at = NULL, cex.axis = 1.15, las = 1)
mtext(side = 1, text = "Tumor purity (%)", line = 3.15, cex = 1.25)
mtext(side = 2, text = "Correlation coefficient", line = 3.15, cex = 1.25)
rect(-100, 1, 1000, 1.75, col="lightgrey", border="black", lwd=1.5)
title(main="All samples", line=-1.35, cex.main=1, font.main=1)
fit = lm(correlation_coefficient ~ tumor_purity, data = tmp.0 %>% filter(tissue=="Breast"))
x = seq(from=0, to=1, length=1000)
y = predict(fit, newdata = list(tumor_purity=x))
points(x*100, y, type="l", col=colours[1], lwd=3, lty=1)
fit = lm(correlation_coefficient ~ tumor_purity, data = tmp.0 %>% filter(tissue=="Lung"))
x = seq(from=0, to=1, length=1000)
y = predict(fit, newdata = list(tumor_purity=x))
points(x*100, y, type="l", col=colours[2], lwd=3, lty=1)
fit = lm(correlation_coefficient ~ tumor_purity, data = tmp.0 %>% filter(tissue=="Prostate"))
x = seq(from=0, to=1, length=1000)
y = predict(fit, newdata = list(tumor_purity=x))
points(x*100, y, type="l", col=colours[3], lwd=3, lty=1)
box(lwd=1.5)
legend(x=70, y=-.2, legend=names(colours), pch=19, col=colours, title="Cancer type", box.lwd=-1, cex=.85)
dev.off()

p = jonckheere.test(x=tmp.0$correlation_coefficient, g=tmp.0$ctdna_fraction_cat, alternative = "increasing", nperm=10000)

pdf(file="../res/rebuttal/ctDNA_Fraction_versus_Correlation_Coefficient_Log2_BoxPlot.pdf", width=8, height=5)
boxplot(correlation_coefficient ~ ctdna_fraction_cat, data=tmp.0, axes=FALSE, frame.plot=FALSE, ylim=c(-.5,1.14375), col="salmon", lwd=1.5)
axis(1, at = 1:5, labels=c("NE", "1-9", "10-29", "30-59", "60-100"), cex.axis = 1.15, las = 1)
axis(2, at = NULL, cex.axis = 1.15, las = 1)
mtext(side = 1, text = "ctDNA fraction (%)", line = 3.15, cex = 1.25)
mtext(side = 2, text = "Correlation coefficient", line = 3.15, cex = 1.25)
rect(-100, 1, 1000, 1.75, col="lightgrey", border="black", lwd=1.5)
text(x=3, y=-.3, labels=paste0("p = ", signif(p$p.value, 3)), cex=1.25)
box(lwd=1.5)
dev.off()

#==================================================
# Log2 ratio plots grail cfdna tumor samples
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
pdf(file="../res/rebuttal/MSK-VB-0008_cfDNA.pdf", width=9.25, height=4)
plot_log3_(x=tmp2, y=tmp, title = "MSK-VB-0008 | cfDNA")
dev.off()

load("../res/rebuttal/uncollapsed_bam/cnvkit/totalcopy/MSK-VL-0056-T.RData")
tmp2 = winsorize(CN, method="mad", tau=3.5, verbose=FALSE)
colnames(tmp2) = c("Chromosome","Position","Log2Ratio")
tmp = undo_(tmp, n=2)
pdf(file="../res/rebuttal/MSK-VL-0056_cfDNA.pdf", width=9.25, height=4)
plot_log3_(x=tmp2, y=tmp, title = "MSK-VL-0056 | cfDNA")
dev.off()

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
pdf(file="../res/rebuttal/MSK-VP-0004_cfDNA.pdf", width=9.25, height=4)
plot_log3_(x=tmp2, y=tmp, title = "MSK-VP-0004 | cfDNA")
dev.off()
 
load("../res/rebuttal/uncollapsed_bam/cnvkit/totalcopy/MSK-VB-0023-T.RData")
tmp2 = winsorize(CN, method="mad", tau=4.5, verbose=FALSE)
colnames(tmp2) = c("Chromosome","Position","Log2Ratio")
tmp = undo_(tmp, n=2)
pdf(file="../res/rebuttal/MSK-VB-0023_cfDNA.pdf", width=9.25, height=4)
plot_log3_(x=tmp2, y=tmp, title = "MSK-VB-0023 | cfDNA")
dev.off()

#==================================================
# Histogram of log2 ratios grail cfdna healthy
# controls
#==================================================
fileNames = dir(path="../res/rebuttal/uncollapsed_bam/cnvkit/totalcopy", pattern="W", full.names=TRUE)
res = foreach (i=1:length(fileNames)) %dopar% {
	load(fileNames[i])
	sex = ifelse(mean(CN[CN[,"Chromosome"]==23,"Log2Ratio"])<(-.5), "M", "F")
	if (sex=="M") {
		fit = density(undo_(tmp, n=2)[,"Log2Ratio"], adjust=.5, from=-1, to=1)
	} else {
		fit = density(undo_(tmp, n=1)[,"Log2Ratio"], adjust=5, from=-1, to=1)
	}
	res = list(x = fit$x,
			   y = fit$y,
			   z = sex)
	return(invisible(res))
}

pdf(file="../res/rebuttal/Histogram_Log2_Ratio_Healthy_M.pdf", width=4, height=6)
plot(0, 0, type="n", xlab="", ylab="", frame.plot=FALSE, axes=FALSE, ylim=c(0,10.625), xlim=c(-2,2))
for (i in 1:length(fileNames)) {
	if (res[[i]]$z=="M") {
		points(res[[i]]$x, res[[i]]$y, type="l", col=transparent_rgb("salmon"), lwd=1.5)
	}
}
abline(v=0, lty=3, col="grey10", lwd=1)
abline(v=-.85, lty=3, col="grey10", lwd=1)
axis(1, at = NULL, cex.axis = 1.15, las = 1)
axis(2, at = NULL, cex.axis = 1.15, las = 1)
mtext(side = 1, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
mtext(side = 2, text = "Density", line = 3.15, cex = 1.25)
rect(-100, 10, 1000, 12, col="lightgrey", border="black", lwd=1.5)
title(main=paste0("Healthy male controls (n=", sum(unlist(lapply(res, function(x) {x$z} ))=="M"), ")"), line=-1.35, cex.main=1, font.main=1)
box(lwd=1.5)
dev.off()

pdf(file="../res/rebuttal/Histogram_Log2_Ratio_Healthy_F.pdf", width=4, height=6)
plot(0, 0, type="n", xlab="", ylab="", frame.plot=FALSE, axes=FALSE, ylim=c(0,10.625), xlim=c(-2,2))
for (i in 1:length(fileNames)) {
	if (res[[i]]$z=="F") {
		points(res[[i]]$x, res[[i]]$y, type="l", col=transparent_rgb("steelblue"), lwd=1.5)
	}
}
abline(v=0, lty=3, col="grey10", lwd=1)
axis(1, at = NULL, cex.axis = 1.15, las = 1)
axis(2, at = NULL, cex.axis = 1.15, las = 1)
mtext(side = 1, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
mtext(side = 2, text = "Density", line = 3.15, cex = 1.25)
rect(-100, 10, 1000, 12, col="lightgrey", border="black", lwd=1.5)
title(main=paste0("Healthy female controls (n=", sum(unlist(lapply(res, function(x) {x$z} ))=="F"), ")"), line=-1.35, cex.main=1, font.main=1)
box(lwd=1.5)
dev.off()


#==================================================
# Breast HER2 amplified cases
#==================================================
'prunesegments.cn' <- function(x, n=10)
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
	
'plotIdeogram' <- function (chrom, cyto.text = FALSE, cex = 0.6, cyto.data, cyto.unit = "bp", unit)
{
	if (chrom == 23) {
			chrom.cytoband <- cyto.data[cyto.data[, 1] == "chrX", ]
	} else {
		if (chrom == 24) {
			chrom.cytoband <- cyto.data[cyto.data[, 1] == "chrY", ]
		} else {
			chrom.cytoband <- cyto.data[cyto.data[, 1] == paste("chr", chrom, sep = ""), ]
		}
	}
	cyto.start <- chrom.cytoband[, 2]
	cyto.end <- chrom.cytoband[, 3]
	scale <- copynumber:::convert.unit(unit1 = unit, unit2 = cyto.unit)
	xleft <- cyto.start * scale
	xright <- cyto.end * scale
	n <- length(xleft)
	chrom.length <- xright[n] - xleft[1]
	stain <- chrom.cytoband[, 5]
	sep.stain <- c("gpos", "gneg", "acen", "gvar", "stalk")
	g <- sapply(sep.stain, grep, x = stain, fixed = TRUE)
	centromere <- g$acen
	stalk <- g$stalk
	col <- rep("", n)
	col[stain == "gneg"] <- "white"
	col[stain == "gpos100"] <- "black"
	col[stain == "gpos75"] <- "gray25"
	col[stain == "gpos50"] <- "gray50"
	col[stain == "gpos25"] <- "gray75"
	col[stain == "stalk"] <- "gray90"
	col[stain == "gvar"] <- "grey"
	col[stain == "acen"] <- "yellow"
	density <- rep(NA, n)
	angle <- rep(45, n)
	density[stain == "gvar"] <- 15
	ylow <- 0
	yhigh <- 1
	plot(x = c(0, max(xright)), y = c(ylow, yhigh), type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, max(xright)), ylim = c(0, 1), xaxs = "i")
	skip.rect <- c(1, n, stalk)
	rect(xleft[-skip.rect], rep(ylow, n - length(skip.rect)), xright[-skip.rect], rep(yhigh, n - length(skip.rect)), 
	col = col[-skip.rect], border = "black", density = density[-skip.rect], 
	angle = angle[-skip.rect])
	draw.roundEdge(start = xleft[1], stop = xright[1], y0 = ylow, y1 = yhigh, col = col[1], bow = "left", density = density[1], angle = angle[1], chrom.length = chrom.length)
	draw.roundEdge(start = xleft[n], stop = xright[n], y0 = ylow, y1 = yhigh, col = col[n], bow = "right", density = density[n], angle = angle[n], chrom.length = chrom.length)
	if (length(stalk) > 0) {
		for (i in 1:length(stalk)) {
			copynumber:::drawStalk(xleft[stalk[i]], xright[stalk[i]], ylow, yhigh, col = col[stalk[i]])
		}
	}
	if (cyto.text) {
		mtext(text = paste(chrom.cytoband[, 4], "-", sep = " "), side = 1, at = (xleft + (xright - xleft)/2), cex = cex, las = 2, adj = 1, xpd = NA)
	}
}

'draw.roundEdge' <- function (start, stop, y0, y1, col, bow, density = NA, angle = 45, lwd = 1, chrom.length)
{
	f <- rep(0, 0)
	f[1] <- 0.001
	i = 1
	half <- y0 + (y1 - y0)/2
	while (f[i] < half) {
		f[i + 1] <- f[i] * 1.3
		i <- i + 1
	}
	f <- f[-length(f)]
	Y <- c(y1, y1, y1 - f, half, y0 + rev(f), y0, y0)
	cyto.length <- stop - start
	share <- cyto.length/chrom.length
	if (share > 0.2) {
		share <- 0.2
	}
	if (bow == "left") {
		round.start <- start + cyto.length * (1 - share)^20
		x <- seq(round.start, start, length.out = (length(f) + 2))
		revx <- rev(x[-length(x)])
		x <- c(x, revx)
		X <- c(stop, x, stop)
	} else {
		if (bow == "right") {
			round.start <- stop - cyto.length * (1 - share)^20
			x <- seq(round.start, stop, length.out = (length(f) + 2))
			revx <- rev(x[-length(x)])
			x <- c(x, revx)
			X <- c(start, x, start)
		}
	}
	polygon(x = X, y = Y, col = col, border = "black", density = density, angle = angle, lwd = lwd)
}

tracker = read_tsv(file=url_master_key, col_types = cols(.default = col_character())) %>%
 		  type_convert() %>%
 		  dplyr::select(PATIENT_ID, GRAIL_ID, DMP_ID, TUMOR_ID, NORMAL_ID, GRAIL_alpha, GRAIL_psi, IMPACT_alpha, IMPACT_psi)
i_cn = read.csv(file="~/share/data/common/cbioportal_repos/msk-impact/msk-impact/data_CNA.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(i_cn) = gsub(pattern=".", replacement="-", x=colnames(i_cn), fixed=TRUE)
rownames(i_cn) = i_cn$Hugo_Symbol
i_cn = i_cn[,tracker$DMP_ID,drop=FALSE]
colnames(i_cn) = tracker$GRAIL_ID

## MSK-IMPACT
index = i_cn["ERBB2",]==2 & grepl("VB", colnames(i_cn))
sample_names = tracker$TUMOR_ID[index]
for (j in 1:length(sample_names)) {
	load(paste0("../res/rebuttal/msk_impact/cnvkit/totalcopy/", sample_names[j], ".RData"))
	pdf(file=paste0("../res/rebuttal/Chr17_", sample_names[j], ".pdf"))
	par(mar = c(6.1, 6, 4.1, 3))
	zz = split.screen(figs=matrix(c(0,1,.15,1, 0,1,0.0775,.4), nrow=2, ncol=4, byrow=TRUE))
	screen(zz[1])
	start = 0
	data(CytoBand)
	end = max(as.numeric(CytoBand[CytoBand[,1]==17,4]))
	plot(1, 1, type="n", xlim=c(start,end), ylim=c(-4,4), xlab="", ylab="", main="", frame.plot=FALSE, axes=FALSE)
	index = CN[,"Chromosome"]==17
	z0 = CN[index,c("Chromosome", "Position", "Log2Ratio"),drop=FALSE]
	z1 = winsorize(data=z0, tau=2.5, k=50)
	z2 = pcf(data=z1, kmin=15, gamma=50)
	tmp = z2[,c("chrom","start.pos","end.pos","mean")]
	colnames(tmp) = c("Chromosome", "Start", "End", "Log2Ratio")
	points(z1[,"pos"], z1[,"Log2Ratio"], type="p", pch=".", cex=3.25, col="grey75")
	for (i in 1:nrow(tmp)) {
		points(c(tmp[i,"Start"], tmp[i,"End"]), rep(tmp[i,"Log2Ratio"],2), type="l", col="red", lwd=4)
	}
	for (i in 1:(nrow(tmp)-1)) {
		points(c(tmp[i,"End"], tmp[i+1,"Start"]), c(tmp[i,"Log2Ratio"],tmp[i+1,"Log2Ratio"]), type="l", col="red", lwd=1)
	}
	abline(h=0, lwd=1)
	axis(2, at = c(-4,-3,-2,-1,0,1,2,3,4), labels=c(-4,-3,-2,-1,0,1,2,3,4), cex.axis = 1.25, las = 1, lwd=1.5, lwd.ticks=1.35)
	mtext(side = 2, text = expression("Log"[2]~"Ratio"), line = 4, cex = 1.5)
	box(lwd=2)
	screen(zz[2])
	assembly = read.csv(file=url_cytogenetic_data, header=FALSE, sep="\t", stringsAsFactors=FALSE)
	plotIdeogram(chrom=17, TRUE, cyto.data = assembly, cex = .75, unit = "bp")
	close.screen(all.screens=TRUE)
	dev.off()
}

## GRAIL
index = i_cn["ERBB2",]==2 & grepl("VB", colnames(i_cn))
sample_names = tracker$GRAIL_ID[index]
for (j in 1:length(sample_names)) {
	load(paste0("../res/rebuttal/uncollapsed_bam/cnvkit/totalcopy/", sample_names[j], "-T.RData"))
	pdf(file=paste0("../res/rebuttal/Chr17_", sample_names[j], ".pdf"))
	par(mar = c(6.1, 6, 4.1, 3))
	zz = split.screen(figs=matrix(c(0,1,.15,1, 0,1,0.0775,.4), nrow=2, ncol=4, byrow=TRUE))
	screen(zz[1])
	start = 0
	data(CytoBand)
	end = max(as.numeric(CytoBand[CytoBand[,1]==17,4]))
	plot(1, 1, type="n", xlim=c(start,end), ylim=c(-4,4), xlab="", ylab="", main="", frame.plot=FALSE, axes=FALSE)
	index = CN[,"Chromosome"]==17
	z0 = CN[index,c("Chromosome", "Position", "Log2Ratio"),drop=FALSE]
	z1 = winsorize(data=z0, tau=4.5)
	z2 = pcf(data=z1, kmin=10, gamma=50)
	tmp = z2[,c("chrom","start.pos","end.pos","mean")]
	colnames(tmp) = c("Chromosome", "Start", "End", "Log2Ratio")
	points(z1[,"pos"], CN[CN$Chromosome==17,"Log2Ratio"], type="p", pch=".", cex=3.25, col="grey75")
	for (i in 1:nrow(tmp)) {
		points(c(tmp[i,"Start"], tmp[i,"End"]), rep(tmp[i,"Log2Ratio"],2), type="l", col="red", lwd=4)
	}
	for (i in 1:(nrow(tmp)-1)) {
		points(c(tmp[i,"End"], tmp[i+1,"Start"]), c(tmp[i,"Log2Ratio"],tmp[i+1,"Log2Ratio"]), type="l", col="red", lwd=1)
	}
	abline(h=0, lwd=1)
	axis(2, at = c(-4,-3,-2,-1,0,1,2,3,4), labels=c(-4,-3,-2,-1,0,1,2,3,4), cex.axis = 1.25, las = 1, lwd=1.5, lwd.ticks=1.35)
	mtext(side = 2, text = expression("Log"[2]~"Ratio"), line = 4, cex = 1.5)
	box(lwd=2)
	screen(zz[2])
	assembly = read.csv(file=url_cytogenetic_data, header=FALSE, sep="\t", stringsAsFactors=FALSE)
	plotIdeogram(chrom=17, TRUE, cyto.data = assembly, cex = .75, unit = "bp")
	close.screen(all.screens=TRUE)
	dev.off()
}
 
 
#==================================================
# Lung MET amplified cases
#==================================================
tracker = read_tsv(file=url_master_key, col_types = cols(.default = col_character())) %>%
 		  type_convert() %>%
 		  dplyr::select(PATIENT_ID, GRAIL_ID, DMP_ID, TUMOR_ID, NORMAL_ID, GRAIL_alpha, GRAIL_psi, IMPACT_alpha, IMPACT_psi)
i_cn = read.csv(file="~/share/data/common/cbioportal_repos/msk-impact/msk-impact/data_CNA.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(i_cn) = gsub(pattern=".", replacement="-", x=colnames(i_cn), fixed=TRUE)
rownames(i_cn) = i_cn$Hugo_Symbol
i_cn = i_cn[,tracker$DMP_ID,drop=FALSE]
colnames(i_cn) = tracker$GRAIL_ID

## MSK-IMPACT
index = i_cn["MET",]==2 & grepl("VL", colnames(i_cn))
sample_names = tracker$TUMOR_ID[index]
for (j in 1:length(sample_names)) {
	load(paste0("../res/rebuttal/msk_impact/cnvkit/totalcopy/", sample_names[j], ".RData"))
	pdf(file=paste0("../res/rebuttal/Chr7_", sample_names[j], ".pdf"))
	par(mar = c(6.1, 6, 4.1, 3))
	zz = split.screen(figs=matrix(c(0,1,.15,1, 0,1,0.0775,.4), nrow=2, ncol=4, byrow=TRUE))
	screen(zz[1])
	start = 0
	data(CytoBand)
	end = max(as.numeric(CytoBand[CytoBand[,1]==7,4]))
	plot(1, 1, type="n", xlim=c(start,end), ylim=c(-4,4), xlab="", ylab="", main="", frame.plot=FALSE, axes=FALSE)
	index = CN[,"Chromosome"]==7
	z0 = CN[index,c("Chromosome", "Position", "Log2Ratio"),drop=FALSE]
	z1 = winsorize(data=z0, tau=2.5, k=50)
	z2 = pcf(data=z1, kmin=15, gamma=50)
	tmp = z2[,c("chrom","start.pos","end.pos","mean")]
	colnames(tmp) = c("Chromosome", "Start", "End", "Log2Ratio")
	points(z1[,"pos"], z1[,"Log2Ratio"], type="p", pch=".", cex=3.25, col="grey75")
	for (i in 1:nrow(tmp)) {
		points(c(tmp[i,"Start"], tmp[i,"End"]), rep(tmp[i,"Log2Ratio"],2), type="l", col="red", lwd=4)
	}
	for (i in 1:(nrow(tmp)-1)) {
		points(c(tmp[i,"End"], tmp[i+1,"Start"]), c(tmp[i,"Log2Ratio"],tmp[i+1,"Log2Ratio"]), type="l", col="red", lwd=1)
	}
	abline(h=0, lwd=1)
	axis(2, at = c(-4,-3,-2,-1,0,1,2,3,4), labels=c(-4,-3,-2,-1,0,1,2,3,4), cex.axis = 1.25, las = 1, lwd=1.5, lwd.ticks=1.35)
	mtext(side = 2, text = expression("Log"[2]~"Ratio"), line = 4, cex = 1.5)
	box(lwd=2)
	screen(zz[2])
	assembly = read.csv(file=url_cytogenetic_data, header=FALSE, sep="\t", stringsAsFactors=FALSE)
	plotIdeogram(chrom=7, TRUE, cyto.data = assembly, cex = .75, unit = "bp")
	close.screen(all.screens=TRUE)
	dev.off()
}

## GRAIL
index = i_cn["MET",]==2 & grepl("VL", colnames(i_cn))
sample_names = tracker$GRAIL_ID[index]
for (j in 1:length(sample_names)) {
	load(paste0("../res/rebuttal/uncollapsed_bam/cnvkit/totalcopy/", sample_names[j], "-T.RData"))
	pdf(file=paste0("../res/rebuttal/Chr7_", sample_names[j], ".pdf"))
	par(mar = c(6.1, 6, 4.1, 3))
	zz = split.screen(figs=matrix(c(0,1,.15,1, 0,1,0.0775,.4), nrow=2, ncol=4, byrow=TRUE))
	screen(zz[1])
	start = 0
	data(CytoBand)
	end = max(as.numeric(CytoBand[CytoBand[,1]==7,4]))
	plot(1, 1, type="n", xlim=c(start,end), ylim=c(-4,4), xlab="", ylab="", main="", frame.plot=FALSE, axes=FALSE)
	index = CN[,"Chromosome"]==7
	z0 = CN[index,c("Chromosome", "Position", "Log2Ratio"),drop=FALSE]
	z1 = winsorize(data=z0, tau=6.5, k=100)
	z2 = pcf(data=z1, kmin=10, gamma=50)
	tmp = z2[,c("chrom","start.pos","end.pos","mean")]
	colnames(tmp) = c("Chromosome", "Start", "End", "Log2Ratio")
	points(z1[,"pos"], z1[,"Log2Ratio"], type="p", pch=".", cex=3.25, col="grey75")
	for (i in 1:nrow(tmp)) {
		points(c(tmp[i,"Start"], tmp[i,"End"]), rep(tmp[i,"Log2Ratio"],2), type="l", col="red", lwd=4)
	}
	for (i in 1:(nrow(tmp)-1)) {
		points(c(tmp[i,"End"], tmp[i+1,"Start"]), c(tmp[i,"Log2Ratio"],tmp[i+1,"Log2Ratio"]), type="l", col="red", lwd=1)
	}
	abline(h=0, lwd=1)
	axis(2, at = c(-4,-3,-2,-1,0,1,2,3,4), labels=c(-4,-3,-2,-1,0,1,2,3,4), cex.axis = 1.25, las = 1, lwd=1.5, lwd.ticks=1.35)
	mtext(side = 2, text = expression("Log"[2]~"Ratio"), line = 4, cex = 1.5)
	box(lwd=2)
	screen(zz[2])
	assembly = read.csv(file=url_cytogenetic_data, header=FALSE, sep="\t", stringsAsFactors=FALSE)
	plotIdeogram(chrom=7, TRUE, cyto.data = assembly, cex = .75, unit = "bp")
	close.screen(all.screens=TRUE)
	dev.off()
}

# # i_bygene = foreach (i=1:nrow(tracker)) %dopar% {
# #  	cat(tracker$GRAIL_ID[i], "\n")
# #  	impact_path = paste0("../res/rebuttal/MSK-IMPACT/facets/cncf/", key_file$TUMOR_ID[i], "_", key_file$NORMAL_ID[i], ".Rdata")
# # 	impact_data = new.env()
# # 	load(impact_path, envir=impact_data)
# # 	
# # 	impact_seg = impact_data$fit$cncf %>%
# # 				 select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
# # 	impact_seg = prune_(x=impact_seg) %>%
# # 				 bind_cols(cn = absolute_(rho=key_file$IMPACT_alpha[i],
# # 										  psi=key_file$IMPACT_psi[i],
# # 										  x=impact_seg$log2)) %>%
# # 				 mutate(n = cumsum(n))
# # 				
# #  	Chromosome = impact_seg[,"chrom"]
# #  	Start = impact_seg[,"start"]
# #  	End = impact_seg[,"end"]
# #  	Calls = impact_seg[,"cn"]
# #  	res = data.frame(Chromosome, Start, End, Calls)
# #  	annot = read.csv(file="~/share/reference/IMPACT410_genes_for_copynumber.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
# #  			select(hgnc_symbol, chr, start_position, end_position) %>%
# #  			rename(Hugo_Symbol = hgnc_symbol,
# #  				   Chromosome = chr,
# #  				   Start = start_position,
# #  				   End = end_position) %>%
# #  			mutate(Chromosome = ifelse(Chromosome %in% "X", 23, Chromosome)) %>%
# #  			arrange(as.numeric(Chromosome), as.numeric(Start))
# #  				   
# #  	annot_by_gene <- annot %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Hugo_Symbol = Hugo_Symbol)
# #  	res_by_segment <- res %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Calls = Calls)
# #  	fo <- findOverlaps(res_by_segment, annot_by_gene)
# #  
# #  	df <- data.frame(Hugo_Symbol=mcols(annot_by_gene)[subjectHits(fo),], Calls=mcols(res_by_segment)[queryHits(fo),])
# #  	Hugo_Symbol = which(duplicated(df$Hugo_Symbol))
# #  	for (j in 1:length(Hugo_Symbol)) {
# #  		index = which(as.character(df$Hugo_Symbol)==as.character(df$Hugo_Symbol[Hugo_Symbol[j]]))
# #  		df[index,2] = mean(df[index,2], na.rm=TRUE)
# #  	}
# #  	df = df %>% filter(!duplicated(Hugo_Symbol))
# #  	df[,2] = round(df[,2])
# #  	res = rep(0, nrow(annot))
# #  	names(res) = annot[,"Hugo_Symbol"]
# #  	res[as.character(df$Hugo_Symbol)] = df$Calls
# #  	return(invisible(res))
# # }
# # i_bygene = do.call(cbind, i_bygene)
# # colnames(i_bygene) = tracker$GRAIL_ID
# # 
# #  
# # g_bygene = foreach (i=1:nrow(tracker)) %dopar% {
# #  	cat(tracker$GRAIL_ID[i], "\n")
# # 	grail_path = paste0("../res/rebuttal/GRAIL/facets/cncf/", key_file$GRAIL_ID[i], "_", key_file$GRAIL_ID[i], "-N.Rdata")
# # 	grail_data = new.env()
# # 	load(grail_path, envir=grail_data)
# # 	
# # 	grail_cn = grail_data$out2$jointseg %>%
# # 			   select(chrom, pos = maploc, log2 = cnlr)
# # 	grail_seg = grail_data$fit$cncf %>%
# # 				select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
# # 	
# # 	fixed_cn = fix_6(grail_cn, grail_seg)
# # 	grail_cn = fixed_cn[[1]]
# # 	grail_seg = fixed_cn[[2]]
# # 	
# # 	grail_seg = prune_(x=grail_seg) %>%
# # 				bind_cols(cn = absolute_(rho=key_file$GRAIL_alpha[i],
# # 										 psi=key_file$GRAIL_psi[i],
# # 										 x=grail_seg$log2)) %>%
# # 				mutate(n = cumsum(n))
# #  	
# #  	Chromosome = grail_seg[,"chrom"]
# #  	Start = grail_seg[,"start"]
# #  	End = grail_seg[,"end"]
# #  	Calls = grail_seg[,"cn"]
# #  	res = data.frame(Chromosome, Start, End, Calls)
# #  	annot = read.csv(file="~/share/reference/IMPACT410_genes_for_copynumber.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
# #  			select(hgnc_symbol, chr, start_position, end_position) %>%
# #  			rename(Hugo_Symbol = hgnc_symbol,
# #  				   Chromosome = chr,
# #  				   Start = start_position,
# #  				   End = end_position) %>%
# #  			mutate(Chromosome = ifelse(Chromosome %in% "X", 23, Chromosome)) %>%
# #  			arrange(as.numeric(Chromosome), as.numeric(Start))
# #  				   
# #  	annot_by_gene <- annot %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Hugo_Symbol = Hugo_Symbol)
# #  	res_by_segment <- res %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Calls = Calls)
# #  	fo <- findOverlaps(res_by_segment, annot_by_gene)
# #  
# #  	df <- data.frame(Hugo_Symbol=mcols(annot_by_gene)[subjectHits(fo),], Calls=mcols(res_by_segment)[queryHits(fo),])
# #  	Hugo_Symbol = which(duplicated(df$Hugo_Symbol))
# #  	for (j in 1:length(Hugo_Symbol)) {
# #  		index = which(as.character(df$Hugo_Symbol)==as.character(df$Hugo_Symbol[Hugo_Symbol[j]]))
# #  		df[index,2] = mean(df[index,2], na.rm=TRUE)
# #  	}
# #  	df = df %>% filter(!duplicated(Hugo_Symbol))
# #  	df[,2] = round(df[,2])
# #  	res = rep(0, nrow(annot))
# #  	names(res) = annot[,"Hugo_Symbol"]
# #  	res[as.character(df$Hugo_Symbol)] = df$Calls
# #  	return(invisible(res))
# # }
# # g_bygene = do.call(cbind, g_bygene)
# # colnames(g_bygene) = tracker$GRAIL_ID
# # 
# # index = c("CRLF2", "HLA-A", "HLA-B", "HLA-C", "AR", "HIST2H3D", "HIST2H3C", "HIST3H3",
# # 		  "HIST1H3A", "HIST1H3B", "HIST1H3C", "HIST1H1C", "HIST1H2BD", "HIST1H3D",
# # 		  "HIST1H3E", "HIST1H3F", "HIST1H3G", "HIST1H3H", "HIST1H3I", "HIST1H3J")
# # i_bygene = i_bygene[!(rownames(i_bygene) %in% index),,drop=FALSE]
# # g_bygene = g_bygene[!(rownames(g_bygene) %in% index),,drop=FALSE]
# # 
# # 
# # annot = read.csv(file="~/share/reference/IMPACT410_genes_for_copynumber.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
# #  			select(hgnc_symbol, chr, start_position, end_position) %>%
# #  			rename(Hugo_Symbol = hgnc_symbol,
# #  				   Chromosome = chr,
# #  				   Start = start_position,
# #  				   End = end_position) %>%
# #  			mutate(Chromosome = ifelse(Chromosome %in% "X", 23, Chromosome)) %>%
# #  			arrange(as.numeric(Chromosome), as.numeric(Start))
# # annot = annot[!(annot$Hugo_Symbol %in% index),,drop=FALSE]
# # 
# # n = 23
# # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# # ColSideColors = rep(NA, nrow(annot))
# # ColSideColors = col_vector[as.numeric(annot[,"Chromosome"])]
# # names(ColSideColors) = as.character(annot[,"Hugo_Symbol"])
# # ColSideColors = ColSideColors[rownames(i_bygene)]
# # RowSideColors = rep(NA, ncol(i_bygene))
# # RowSideColors[grep("VB", colnames(i_bygene), fixed=TRUE)] = "salmon"
# # RowSideColors[grep("VL", colnames(i_bygene), fixed=TRUE)] = "#FDAE61"
# # RowSideColors[grep("VP", colnames(i_bygene), fixed=TRUE)] = "#ABDDA4"
# # 
# # index = is.na(ColSideColors)
# # ColSideColors = ColSideColors[!index]
# # i_bygene = i_bygene[!index,,drop=FALSE]
# # g_bygene = g_bygene[!index,,drop=FALSE]
# # 
# # for (i in 1:nrow(tracker)) {
# # 	psi = round(tracker$IMPACT_psi[i])
# #  	psi = ifelse(psi==1, 2, psi)
# #  	x = i_bygene[,i]
# #  	x2 = rep(0, length(x))
# #  	if (psi==2) {
# # 	 	x2[x<(psi-1)] = -1
# # 		x2[x>(psi+3)] = 1
# # 	} else if (psi==3) {
# # 		x2[x<(psi-2)] = -1
# # 		x2[x>(psi+4)] = 1
# # 	} else {
# # 		x2[x<(psi-3)] = -1
# # 		x2[x>(psi+5)] = 1
# # 	}
# # 	i_bygene[,i] = x2
# # 	
# # 	psi = round(tracker$GRAIL_psi[i])
# #  	psi = ifelse(psi==1, 2, psi)
# #  	x = g_bygene[,i]
# #  	x2 = rep(0, length(x))
# #  	if (psi==2) {
# # 	 	x2[x<(psi-1)] = -1
# # 		x2[x>(psi+3)] = 1
# # 	} else if (psi==3) {
# # 		x2[x<(psi-2)] = -1
# # 		x2[x>(psi+4)] = 1
# # 	} else {
# # 		x2[x<(psi-3)] = -1
# # 		x2[x>(psi+4)] = 1
# # 	}
# # 	g_bygene[,i] = x2
# # }
# # 
# # pdf(file="../res/rebuttal/Heatmap_CN_tumor_abs_copy_all_genes.pdf", width=8)
# # hmi = heatmap(t(i_bygene), Colv=NA, Rowv=as.dendrogram(hclust(dist(t(i_bygene)))),
# # 			  col = c("steelblue", "white", "red"),
# # 			  scale="none",
# #  			  RowSideColors=RowSideColors,
# #  			  ColSideColors=ColSideColors,
# #  			  labRow=NA,
# #  			  labCol=NA)
# # dev.off()
# # 
# # pdf(file="../res/rebuttal/Heatmap_CN_cfDNA_abs_copy_all_genes.pdf", width=8)
# # hmg = heatmap(t(g_bygene), Colv=NA, Rowv=as.dendrogram(hclust(dist(t(i_bygene)))),
# # 			  col = c("steelblue", "white", "red"),
# # 			  scale = "none",
# #  			  RowSideColors=RowSideColors,
# #  			  ColSideColors=ColSideColors,
# #  			  labRow=NA,
# #  			  labCol=NA)
# # dev.off()
# # 
# # 
# # 
# # 

# # #==================================================
# # # Updated ploidy
# # #==================================================
# # if (TRUE) { ploidy = foreach (i=1:nrow(key_file)) %dopar% {
# # 	print(key_file$GRAIL_ID[i])
# # 	impact_path = paste0("../res/rebuttal/MSK-IMPACT/facets/cncf/", key_file$TUMOR_ID[i], "_", key_file$NORMAL_ID[i], ".Rdata")
# # 	impact_data = new.env()
# # 	load(impact_path, envir=impact_data)
# # 	impact_cn = impact_data$out2$jointseg %>%
# # 			    select(chrom, pos = maploc, log2 = cnlr)
# # 	impact_seg = impact_data$fit$cncf %>%
# # 				 select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
# # 	impact_seg = prune_(x=impact_seg) %>%
# # 				 bind_cols(cn = absolute_(rho=key_file$IMPACT_alpha[i],
# # 										  psi=key_file$IMPACT_psi[i],
# # 										  x=impact_seg$log2)) %>%
# # 				mutate(n = cumsum(n))
# # 	impact_psi = (t(impact_seg[,"end"] - impact_seg[,"start"])%*%impact_seg[,"cn"])/sum(impact_seg[,"end"] - impact_seg[,"start"])
# # 	grail_path = paste0("../res/rebuttal/GRAIL/facets/cncf/", key_file$GRAIL_ID[i], "_", key_file$GRAIL_ID[i], "-N.Rdata")
# # 	grail_data = new.env()
# # 	load(grail_path, envir=grail_data)
# # 	grail_cn = grail_data$out2$jointseg %>%
# # 			   select(chrom, pos = maploc, log2 = cnlr)
# # 	grail_seg = grail_data$fit$cncf %>%
# # 				select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
# # 	
# # 	fixed_cn = fix_6(grail_cn, grail_seg)
# # 	grail_cn = fixed_cn[[1]]
# # 	grail_seg = fixed_cn[[2]]
# # 	
# # 	grail_seg = prune_(x=grail_seg) %>%
# # 				bind_cols(cn = absolute_(rho=key_file$GRAIL_alpha[i],
# # 										 psi=key_file$GRAIL_psi[i],
# # 										 x=grail_seg$log2)) %>%
# # 				mutate(n = cumsum(n))
# # 	grail_psi = (t(grail_seg[,"end"] - grail_seg[,"start"])%*%grail_seg[,"cn"])/sum(grail_seg[,"end"] - grail_seg[,"start"])
# # 	return(c(grail_psi, impact_psi))
# # } }
# # ploidy = do.call(rbind, ploidy)
# # 
# # tmp = cbind(tracker, pa)
# # tmp[,"GRAIL_psi"] = ploidy[,1]
# # tmp[,"IMPACT_psi"] = ploidy[,2]
# # 
# # tmp.0 = tmp %>%
# # 		mutate(GRAIL_alpha = ifelse(GRAIL_alpha==1, 0, GRAIL_alpha)) %>%
# # 		mutate(IMPACT_alpha = ifelse(IMPACT_alpha==1, 0, IMPACT_alpha)) %>%
# # 		mutate(Cat = "Both estimate available") %>%
# # 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha == 0, "No estimate in both", Cat)) %>%
# # 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha != 0, "No estimate in cfDNA", Cat)) %>%
# # 		mutate(Cat = ifelse(GRAIL_alpha != 0 & IMPACT_alpha == 0, "No estimate in Biopsy", Cat)) %>%
# # 		mutate(Purity = "Purity")
# # 		
# # plot.0 = ggplot(tmp.0, aes(y = GRAIL_alpha, x = IMPACT_alpha, shape = Tissue, fill = Cat)) +
# # 		 geom_abline(slope = 1, color = "goldenrod3", linetype = 1) +
# # 		 geom_point(alpha = .8, size = 2.5) +
# # 		 scale_fill_manual(values = c("Both estimate available"="salmon", "No estimate in both"="#FDAE61", "No estimate in cfDNA"="#ABDDA4", "No estimate in Biopsy"="steelblue")) +
# # 		 scale_shape_manual(values = c(24, 21, 22)) +
# # 		 theme_bw(base_size=15) +
# # 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.26, 0.7), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# # 		 labs(y="\ncfDNA\n", x="\nBiopsy\n") +
# # 		 coord_cartesian(xlim=c(0,1), ylim = c(0, 1)) +
# # 		 facet_wrap(~Purity) +
# # 		 guides(shape=guide_legend(title=c("Tissue"), override.aes=list(fill="black"))) +
# # 		 guides(fill=guide_legend(title=c("Purity")))
# # 		 
# # 		 
# # pdf(file="../res/rebuttal/Comparison_Purity.pdf", width=6, height=6)
# # print(plot.0)
# # dev.off()
# # 
# # tmp.0 = tmp %>%
# # 		mutate(GRAIL_alpha = ifelse(GRAIL_alpha==1, 0, GRAIL_alpha)) %>%
# # 		mutate(IMPACT_alpha = ifelse(IMPACT_alpha==1, 0, IMPACT_alpha)) %>%
# # 		mutate(Cat = "Both estimate available") %>%
# # 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha == 0, "No estimate in both", Cat)) %>%
# # 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha != 0, "No estimate in cfDNA", Cat)) %>%
# # 		mutate(Cat = ifelse(GRAIL_alpha != 0 & IMPACT_alpha == 0, "No estimate in Biopsy", Cat)) %>%
# # 		mutate(Ploidy = "Ploidy")
# # 		
# # plot.0 = ggplot(tmp.0, aes(y = GRAIL_psi, x = IMPACT_psi, shape = Tissue, fill = Cat)) +
# # 		 geom_abline(slope = 1, color = "goldenrod3", linetype = 1) +
# # 		 geom_point(alpha = .8, size = 2.5) +
# #  		 scale_fill_manual(values = c("Both estimate available"="salmon", "No estimate in both"="#FDAE61", "No estimate in cfDNA"="#ABDDA4", "No estimate in Biopsy"="steelblue")) +
# # 		 scale_shape_manual(values = c(24, 21, 22)) +
# # 		 theme_bw(base_size=15) +
# # 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.26, 0.7), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# # 		 labs(y="\ncfDNA\n", x="\nBiopsy\n") +
# # 		 coord_cartesian(xlim=c(1.5,4), ylim = c(1.5,4)) +
# # 		 facet_wrap(~Ploidy) +
# # 		 guides(shape=guide_legend(title=c("Tissue"), override.aes=list(fill="black"))) +
# # 		 guides(fill=guide_legend(title=c("Purity")))
# # 		 
# # pdf(file="../res/rebuttal/Comparison_Ploidy.pdf", width=6, height=6)
# # print(plot.0)
# # dev.off()
# # 
# # 
# # tmp.0 = tmp %>%
# # 		mutate(Tissue = "All samples")
# # 
# # plot.0 = ggplot(tmp.0, aes(x = 100*pa)) +
# # 		 geom_histogram(color="black", fill="#2B83BA") +
# # 		 coord_cartesian(xlim=c(0,100)) +
# # 		 theme_bw(base_size=15) +
# # 		 facet_wrap(~Tissue) +
# # 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# # 		 labs(y="\nFrequency\n", x="\n Agreement (%)\n")
# # 		 
# # pdf(file="../res/rebuttal/Percent_Agreement_Combined.pdf", width=6.5)
# # print(plot.0)
# # dev.off()
# # 
# # plot.0 = ggplot(tmp, aes(x = 100*pa)) +
# # 		 geom_histogram(color="black", fill="#2B83BA") +
# # 		 coord_cartesian(xlim=c(0,100)) +
# # 		 theme_bw(base_size=15) +
# # 		 facet_wrap(~Tissue) +
# # 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# # 		 labs(y="\nFrequency\n", x="\n Agreement (%)\n")
# # 		 
# # pdf(file="../res/rebuttal/Percent_Agreement_Tissue.pdf", width=15)
# # print(plot.0)
# # dev.off()
# # 
# # tmp.0 = tmp %>%
# # 		mutate(GRAIL_alpha = ifelse(GRAIL_alpha==1, 0, GRAIL_alpha)) %>%
# # 		mutate(IMPACT_alpha = ifelse(IMPACT_alpha==1, 0, IMPACT_alpha)) %>%
# # 		mutate(Cat = "Both estimate available") %>%
# # 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha == 0, "No estimate in both", Cat)) %>%
# # 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha != 0, "No estimate in cfDNA", Cat)) %>%
# # 		mutate(Cat = ifelse(GRAIL_alpha != 0 & IMPACT_alpha == 0, "No estimate in Biopsy", Cat)) %>%
# # 		mutate(Ploidy = "Agreement versus cfDNA ploidy")
# # 		
# # 
# # plot.0 = ggplot(tmp.0, aes(y = GRAIL_psi, x = 100*pa, shape = Tissue, fill = Cat)) +
# # 		 geom_point(alpha = .8, size = 2.5) +
# # 		 scale_fill_manual(values = c("Both estimate available"="salmon", "No estimate in both"="#FDAE61", "No estimate in cfDNA"="#ABDDA4", "No estimate in Biopsy"="steelblue")) +
# # 		 scale_shape_manual(values = c(24, 21, 22)) +
# # 		 theme_bw(base_size=15) +
# # 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.26, 0.7), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# # 		 labs(y="\nPloidy in cfDNA\n", x="\nAgreement (%)\n") +
# # 		 coord_cartesian(xlim=c(0,100), ylim = c(1.5,4)) +
# # 		 facet_wrap(~Ploidy) +
# # 		 guides(shape=guide_legend(title=c("Tissue"), override.aes=list(fill="black"))) +
# # 		 guides(fill=guide_legend(title=c("Purity")))
# # 		 
# # pdf(file="../res/rebuttal/Agreement_Ploidy_cfDNA.pdf", width=6, height=6)
# # print(plot.0)
# # dev.off()
# # 
# # tmp.0 = tmp %>%
# # 		mutate(GRAIL_alpha = ifelse(GRAIL_alpha==1, 0, GRAIL_alpha)) %>%
# # 		mutate(IMPACT_alpha = ifelse(IMPACT_alpha==1, 0, IMPACT_alpha)) %>%
# # 		mutate(Cat = "Both estimate available") %>%
# # 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha == 0, "No estimate in both", Cat)) %>%
# # 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha != 0, "No estimate in cfDNA", Cat)) %>%
# # 		mutate(Cat = ifelse(GRAIL_alpha != 0 & IMPACT_alpha == 0, "No estimate in Biopsy", Cat)) %>%
# # 		mutate(Ploidy = "Agreement versus tumor ploidy")
# # 		
# # 
# # plot.0 = ggplot(tmp.0, aes(y = IMPACT_psi, x = 100*pa, shape = Tissue, fill = Cat)) +
# # 		 geom_point(alpha = .8, size = 2.5) +
# # 		 scale_fill_manual(values = c("Both estimate available"="salmon", "No estimate in both"="#FDAE61", "No estimate in cfDNA"="#ABDDA4", "No estimate in Biopsy"="steelblue")) +
# # 		 scale_shape_manual(values = c(24, 21, 22)) +
# # 		 theme_bw(base_size=15) +
# # 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.8, 0.7), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# # 		 labs(y="\nPloidy in tumor\n", x="\nAgreement (%)\n") +
# # 		 coord_cartesian(xlim=c(0,100), ylim = c(1.5,4)) +
# # 		 facet_wrap(~Ploidy) +
# # 		 guides(shape=guide_legend(title=c("Tissue"), override.aes=list(fill="black"))) +
# # 		 guides(fill=guide_legend(title=c("Purity")))
# # 		 
# # pdf(file="../res/rebuttal/Agreement_Ploidy_biopsy.pdf", width=6, height=6)
# # print(plot.0)
# # dev.off()
# # 
# # tmp.0 = tmp %>%
# # 		mutate(GRAIL_alpha = ifelse(GRAIL_alpha==1, 0, GRAIL_alpha)) %>%
# # 		mutate(IMPACT_alpha = ifelse(IMPACT_alpha==1, 0, IMPACT_alpha)) %>%
# # 		mutate(Cat = "Both estimate available") %>%
# # 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha == 0, "No estimate in both", Cat)) %>%
# # 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha != 0, "No estimate in cfDNA", Cat)) %>%
# # 		mutate(Cat = ifelse(GRAIL_alpha != 0 & IMPACT_alpha == 0, "No estimate in Biopsy", Cat)) %>%
# # 		mutate(Ploidy = "Agreement versus ploidy difference")
# # 		
# # plot.0 = ggplot(tmp.0, aes(y = abs(IMPACT_psi - GRAIL_psi), x = 100*pa, shape = Tissue, fill = Cat)) +
# # 		 geom_point(alpha = .8, size = 2.5) +
# # 		 scale_fill_manual(values = c("Both estimate available"="salmon", "No estimate in both"="#FDAE61", "No estimate in cfDNA"="#ABDDA4", "No estimate in Biopsy"="steelblue")) +
# # 		 scale_shape_manual(values = c(24, 21, 22)) +
# # 		 theme_bw(base_size=15) +
# # 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.8, 0.7), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# # 		 labs(y="\n| cfDNA - Tumor |\n", x="\nAgreement (%)\n") +
# # 		 coord_cartesian(xlim=c(0,100), ylim = c(0,2)) +
# # 		 facet_wrap(~Ploidy) +
# # 		 guides(shape=guide_legend(title=c("Tissue"), override.aes=list(fill="black"))) +
# # 		 guides(fill=guide_legend(title=c("Purity")))
# # 		 
# # pdf(file="../res/rebuttal/Agreement_Ploidy_biopsy_cfDNA.pdf", width=6, height=6)
# # print(plot.0)
# # dev.off()
# # 
# # tmp.0 = tmp %>%
# # 		mutate(GRAIL_alpha = ifelse(GRAIL_alpha==1, 0, GRAIL_alpha)) %>%
# # 		mutate(IMPACT_alpha = ifelse(IMPACT_alpha==1, 0, IMPACT_alpha)) %>%
# # 		mutate(Cat = "Both estimate available") %>%
# # 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha == 0, "No estimate in both", Cat)) %>%
# # 		mutate(Cat = ifelse(GRAIL_alpha == 0 & IMPACT_alpha != 0, "No estimate in cfDNA", Cat)) %>%
# # 		mutate(Cat = ifelse(GRAIL_alpha != 0 & IMPACT_alpha == 0, "No estimate in Biopsy", Cat)) %>%
# # 		mutate(Ploidy = "Agreement versus purity difference")
# # 		
# # 
# # plot.0 = ggplot(tmp.0, aes(y = abs(IMPACT_alpha - GRAIL_alpha), x = 100*pa, shape = Tissue, fill = Cat)) +
# # 		 geom_point(alpha = .8, size = 2.5) +
# # 		 scale_fill_manual(values = c("Both estimate available"="salmon", "No estimate in both"="#FDAE61", "No estimate in cfDNA"="#ABDDA4", "No estimate in Biopsy"="steelblue")) +
# # 		 scale_shape_manual(values = c(24, 21, 22)) +
# # 		 theme_bw(base_size=15) +
# # 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.8, 0.7), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
# # 		 labs(y="\n| cfDNA - Tumor |\n", x="\nAgreement (%)\n") +
# # 		 coord_cartesian(xlim=c(0,100), ylim = c(0,1)) +
# # 		 facet_wrap(~Ploidy) +
# # 		 guides(shape=guide_legend(title=c("Tissue"), override.aes=list(fill="black"))) +
# # 		 guides(fill=guide_legend(title=c("Purity")))
# # 		 
# # pdf(file="../res/rebuttal/Agreement_Purity_biopsy_cfDNA.pdf", width=6, height=6)
# # print(plot.0)
# # dev.off()
# 
# 
# # for (i in 1:nrow(tracker)) {
# # 	x = g_cn_jrf[,i]
# # 	y = ifelse(i_cn[,i]==1, 1, 0)
# # 	z = rep(0, length(x))
# # 	if (length(unique(y))==2) {
# # 		pred = prediction(x,y)
# # 		perf = performance(pred, measure = "acc")
# # 		ind = which.max( slot(perf, "y.values")[[1]] )
# # 		co = slot(perf, "x.values")[[1]][ind]
# # 	} else {
# # 		co = Inf
# # 	}
# # 	z[x>co] = 1
# # 	
# # 	y = ifelse(i_cn[,i]==-1, 0, 1)
# #     if (length(unique(y))==2) {
# #     	pred = prediction(x,y)
# #     	perf = performance(pred, measure = "acc")
# #         ind = which.max( slot(perf, "y.values")[[1]] )
# #         co = slot(perf, "x.values")[[1]][ind]
# #     } else {
# #     	co = -Inf
# #     }
# #     z[x<co] = -1
# #     g_cn_jrf[,i] = z
# # }
# # 

# # i_cn = i_cn[rownames(g_cn),,drop=FALSE]
# # i_cn[i_cn>0] = 1
# # i_cn[i_cn<0] = -1
# # i_cn[is.na(i_cn)] = 0
# # i_cn = as.matrix(i_cn)
# # 
# # for (i in 1:nrow(tracker)) {
# # 	x = g_cn_cmo[,i]
# # 	y = ifelse(i_cn[,i]==1, 1, 0)
# # 	z = rep(0, length(x))
# # 	if (length(unique(y))==2) {
# # 		pred = prediction(x,y)
# # 		perf = performance(pred, measure = "acc")
# # 		ind = which.max( slot(perf, "y.values")[[1]] )
# # 		co = slot(perf, "x.values")[[1]][ind]
# # 	} else {
# # 		co = Inf
# # 	}
# # 	z[x>co] = 1
# # 	
# # 	y = ifelse(i_cn[,i]==-1, 0, 1)
# #     if (length(unique(y))==2) {
# #     	pred = prediction(x,y)
# #     	perf = performance(pred, measure = "acc")
# #         ind = which.max( slot(perf, "y.values")[[1]] )
# #         co = slot(perf, "x.values")[[1]][ind]
# #     } else {
# #     	co = -Inf
# #     }
# #     z[x<co] = -1
# #     g_cn_cmo[,i] = z
# # }
# # 
# # g_cn = g_cn_jrf
# # g_cn[g_cn_cmo==-1] = -1
# # g_cn[g_cn_cmo==1] = 1
# # 
# # index = which.max(apply(g_cn, 2, function(x) { sum(x==-1) }))
# # g_cn[g_cn[,index]==-1,index] = 0
# # 
# # i_cn = read.csv(file="~/share/data/common/cbioportal_repos/msk-impact/msk-impact/data_CNA.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
# #            filter(Hugo_Symbol %in% rownames(g_cn))
# # colnames(i_cn) = gsub(pattern=".", replacement="-", x=colnames(i_cn), fixed=TRUE)
# # rownames(i_cn) = i_cn$Hugo_Symbol
# # i_cn = i_cn[,tracker$DMP_ID,drop=FALSE]
# # colnames(i_cn) = tracker$GRAIL_ID
# # i_cn = i_cn[rownames(g_cn),,drop=FALSE]
# # i_cn[i_cn>0] = 1
# # i_cn[i_cn<0] = -1
# # i_cn[is.na(i_cn)] = 0
# # i_cn = as.matrix(i_cn)
# # 
# # annot = read.csv(file="~/share/reference/IMPACT410_genes_for_copynumber.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
# # 				 select(hgnc_symbol, chr, start_position, end_position) %>%
# # 				 rename(Hugo_Symbol = hgnc_symbol,
# # 				   		Chromosome = chr,
# # 				   		Start = start_position,
# # 				   		End = end_position) %>%
# # 				 mutate(Chromosome = ifelse(Chromosome %in% "X", 23, Chromosome))
# # n = 23
# # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# # ColSideColors = rep(NA, nrow(annot))
# # ColSideColors = col_vector[as.numeric(annot[,"Chromosome"])]
# # names(ColSideColors) = as.character(annot[,"Hugo_Symbol"])
# # ColSideColors = ColSideColors[rownames(i_cn)]
# # RowSideColors = rep(NA, ncol(i_cn))
# # RowSideColors[grep("VB", colnames(i_cn), fixed=TRUE)] = "salmon"
# # RowSideColors[grep("VL", colnames(i_cn), fixed=TRUE)] = "#FDAE61"
# # RowSideColors[grep("VP", colnames(i_cn), fixed=TRUE)] = "#ABDDA4"
# # 
# # index = is.na(ColSideColors)
# # ColSideColors = ColSideColors[!index]
# # i_cn = i_cn[!index,,drop=FALSE]
# # g_cn = g_cn[!index,,drop=FALSE]
# # 
# # pdf(file="../res/rebuttal/Heatmap_CN_IMPACT_All_Genes.pdf", width=8)
# # hmi = heatmap(t(i_cn), Colv=NA, Rowv=as.dendrogram(hclust(dist(t(i_cn)))), col=c("blue", "grey99", "red"), scale="none",
# # 			  RowSideColors=RowSideColors,
# # 			  ColSideColors=ColSideColors,
# # 			  labRow=NA,
# # 			  labCol=NA)
# # dev.off()
# # 
# # pdf(file="../res/rebuttal/Heatmap_CN_GRAIL_All_Genes.pdf", width=8)
# # hmg = heatmap(t(g_cn), Colv=NA, Rowv=as.dendrogram(hclust(dist(t(i_cn)))), col=c("blue", "grey99", "red"), scale="none",
# # 			  RowSideColors=RowSideColors,
# # 			  ColSideColors=ColSideColors,
# # 			  labRow=NA,
# # 			  labCol=NA)
# # dev.off()
# # 
# # index = apply(i_cn, 1, function(x) {sum(x==0)==length(x)})
# # ColSideColors = ColSideColors[!index]
# # i_cn = i_cn[!index,,drop=FALSE]
# # g_cn = g_cn[!index,,drop=FALSE]
# # 
# # pdf(file="../res/rebuttal/Heatmap_CN_IMPACT_HD_and_Amp_Only.pdf", width=8)
# # hmi = heatmap(t(i_cn), Colv=NA, Rowv=as.dendrogram(hclust(dist(t(i_cn)))), col=c("blue", "grey99", "red"), scale="none",
# # 			  RowSideColors=RowSideColors,
# # 			  ColSideColors=ColSideColors,
# # 			  labRow=NA,
# # 			  labCol=NA)
# # dev.off()
# # 
# # pdf(file="../res/rebuttal/Heatmap_CN_GRAIL_HD_and_Amp_Only.pdf", width=8)
# # hmg = heatmap(t(g_cn), Colv=NA, Rowv=as.dendrogram(hclust(dist(t(i_cn)))), col=c("blue", "grey99", "red"), scale="none",
# # 			  RowSideColors=RowSideColors,
# # 			  ColSideColors=ColSideColors,
# # 			  labRow=NA,
# # 			  labCol=NA)
# # dev.off()
