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
# update msk_impact tumor (alpha, psi)
#==================================================
key_file = read_tsv(file=url_master_key, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   dplyr::select(PATIENT_ID, GRAIL_ID, DMP_ID, TUMOR_ID, NORMAL_ID, GRAIL_alpha, GRAIL_psi, IMPACT_alpha, IMPACT_psi)

if (FALSE) { foreach (i=1:nrow(key_file)) %dopar% {
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

if (FALSE) { foreach (i=1:nrow(key_file)) %dopar% {
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
if (FALSE) { foreach (i=1:nrow(key_file)) %dopar% {
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

plot.0 = ggplot(tmp.0 %>% mutate(facet="All samples (n = 124)"), aes(x = ctdna_fraction_cat, y = correlation_coefficient)) + 
		 geom_boxplot(alpha=1, outlier.size=2.5, outlier.shape=21, fill="salmon") + 
		 facet_wrap(~facet) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15)) +
		 labs(x="\nctDNA fraction (%)", y="Correlation coefficient\n")
		 
pdf(file="../res/rebuttal/ctDNA_Fraction_versus_Correlation_Coefficient_Log2_BoxPlot.pdf", width=7, height=5)
print(plot.0)
dev.off()

#==================================================
# ROC curves
#==================================================
tracker = read_tsv(file=url_master_key, col_types = cols(.default = col_character())) %>%
		  type_convert() %>%
 		  dplyr::select(PATIENT_ID, GRAIL_ID, DMP_ID, TUMOR_ID, NORMAL_ID, GRAIL_alpha, GRAIL_psi, IMPACT_alpha, IMPACT_psi)
		   
ctDNA_fraction = read_csv(file=url_ctdna_frac, col_types = cols(.default = col_character())) %>%
		   		 type_convert() %>%
		   		 rename(GRAIL_ID = ID)
		   		 
tracker = left_join(tracker, ctDNA_fraction, by="GRAIL_ID") %>%
 		  filter(!is.na(ctdna_frac)) %>%
		  filter(ctdna_frac>.1)
 
sse = foreach (i=1:nrow(tracker)) %dopar% {
  	cat(tracker$GRAIL_ID[i], "\n")
  	if (!is.na(tracker$ctdna_frac[i]) & tracker$ctdna_frac[i]>.1) {
 	  	grail_path = paste0("../res/rebuttal/uncollapsed_bam/cnvkit/totalcopy/", tracker$GRAIL_ID[i], "-T.RData")
 	 	grail_data = new.env()
 	 	load(grail_path, envir=grail_data)
  	 	grail_seg = grail_data$tmp %>%
  				    dplyr::select(chrom=Chromosome, start = Start, end = End, log2 = Log2Ratio, n=N) %>%
  				    filter(chrom<23)
  		grail_seg = prune_(x=grail_seg)
  		sse = ploidy_(x=grail_seg$log2, y=tracker$ctdna_frac[i], w=grail_seg$n)
  	} else {
  		sse = rep(NA, 100)
  	}
  	return(sse)
}
min_sse = unlist(lapply(sse, min))
index = rep(NA, length(sse))
index[!is.na(min_sse)] = unlist(lapply(sse, which.min))
ploidy = seq(from=1.5, to=3.5, length=100)[index]
tracker$GRAIL_psi = ploidy
tracker = tracker %>%
	   	  filter(!is.na(GRAIL_psi) & !is.na(ctdna_frac))
 		   
i_cn = foreach (i=1:nrow(tracker)) %dopar% {
 	cat(tracker$GRAIL_ID[i], "\n")
  	path = paste0("../res/rebuttal/msk_impact/cnvkit/totalcopy/", tracker$TUMOR_ID[i], ".RData")
 	data = new.env()
 	load(path, envir=data)
  	seg = data$tmp %>%
  		  dplyr::select(chrom=Chromosome, start = Start, end = End, log2 = Log2Ratio, n=N) %>%
  		  filter(chrom<23)
  	seg = prune_(x=seg)
  	a_cn = round(absolute_(rho=tracker$IMPACT_alpha[i], psi=tracker$GRAIL_psi[i], gamma=.85, seg$log2))
 	c_cn = rep(0, length(a_cn))
  	if (round(key_file$IMPACT_psi[i])==2) {
  		c_cn[a_cn<0.5] = -1
   		c_cn[a_cn>6] = 1
  	} else if (round(key_file$IMPACT_psi[i])==3) {
  		c_cn[a_cn<1] = -1
   		c_cn[a_cn>7] = 1
  	} else if (round(key_file$IMPACT_psi[i])==4) {
  		c_cn[a_cn<1] = -1
   		c_cn[a_cn>8] = 1
 	}
  	
 	sample_names = c("MSK-VB-0066", "MSK-VP-0005", "MSK-VP-0025", "MSK-VP-0047")
 
  	if (tracker$GRAIL_ID[i] %in% sample_names) {
  		c_cn = rep(0, length(a_cn))
  		c_cn[a_cn<0] = -1
   		c_cn[a_cn>7] = 1
   	}
   	
 	seg = cbind(seg, a_cn, c_cn)
  	
  	Chromosome = seg[,"chrom"]
  	Start = seg[,"start"]
  	End = seg[,"end"]
  	Calls = seg[,"c_cn"]
  	res = data.frame(Chromosome, Start, End, Calls)
  	annot = read.csv(file="~/share/reference/IMPACT410_genes_for_copynumber.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
  			dplyr::select(hgnc_symbol, chr, start_position, end_position) %>%
  			rename(Hugo_Symbol = hgnc_symbol,
  				   Chromosome = chr,
  				   Start = start_position,
  				   End = end_position) %>%
  			mutate(Chromosome = ifelse(Chromosome %in% "X", 23, Chromosome)) %>%
  			arrange(as.numeric(Chromosome), as.numeric(Start))
  				   
  	annot_by_gene <- annot %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Hugo_Symbol = Hugo_Symbol)
  	res_by_segment <- res %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Calls = Calls)
  	fo <- findOverlaps(res_by_segment, annot_by_gene)
  
  	df <- data.frame(Hugo_Symbol=mcols(annot_by_gene)[subjectHits(fo),], Calls=mcols(res_by_segment)[queryHits(fo),])
  	Hugo_Symbol = which(duplicated(df$Hugo_Symbol))
  	for (j in 1:length(Hugo_Symbol)) {
  		index = which(as.character(df$Hugo_Symbol)==as.character(df$Hugo_Symbol[Hugo_Symbol[j]]))
  		df[index,2] = max(df[index,2], na.rm=TRUE)
  	}
  	df = df %>% filter(!duplicated(Hugo_Symbol))
  	df[,2] = round(df[,2])
  	res = rep(0, nrow(annot))
  	names(res) = annot[,"Hugo_Symbol"]
  	res[as.character(df$Hugo_Symbol)] = df$Calls
  	return(invisible(res))
}
i_cn = do.call(cbind, i_cn)
colnames(i_cn) = tracker$GRAIL_ID
i_cn["ERBB2","MSK-VB-0044"] = 1
 
g_cn = foreach (i=1:nrow(tracker)) %dopar% {
	cat(tracker$GRAIL_ID[i], "\n")
 	path = paste0("../res/rebuttal/uncollapsed_bam/cnvkit/totalcopy/", tracker$GRAIL_ID[i], "-T.RData")
 	data = new.env()
 	load(path, envir=data)
  	seg = data$tmp %>%
  		  dplyr::select(chrom=Chromosome, start = Start, end = End, log2 = Log2Ratio, n=N) %>%
  		  filter(chrom<23)
  	seg = prune_(x=seg)
  	a_cn = round(absolute_(rho=tracker$ctdna_frac[i], psi=tracker$GRAIL_psi[i], gamma=.85, seg$log2))
 	c_cn = rep(0, length(a_cn))
  	if (round(tracker$GRAIL_psi[i])==2) {
  		c_cn[a_cn<0.5] = -1
  		c_cn[a_cn>6] = 1
  	} else if (round(tracker$GRAIL_psi[i])==3) {
  		c_cn[a_cn<1] = -1
  		c_cn[a_cn>7] = 1
  	} else if (round(tracker$IMPACT_psi[i])==4) {
  		c_cn[a_cn<1] = -1
   		c_cn[a_cn>8] = 1
  	}
  	
  	sample_names = c("MSK-VB-0024", "MSK-VP-0028", "MSK-VB-0032", "MSK-VB-0046", "MSK-VP-0031")
 
  	if (tracker$GRAIL_ID[i] %in% sample_names) {
  		c_cn = rep(0, length(a_cn))
  		c_cn[a_cn<0] = -1
   		c_cn[a_cn>7] = 1
   	}
   	
   	
  	
  	seg = cbind(seg, a_cn, c_cn)
  	
  	Chromosome = seg[,"chrom"]
  	Start = seg[,"start"]
  	End = seg[,"end"]
  	Calls = seg[,"a_cn"]
  	res = data.frame(Chromosome, Start, End, Calls)
  	annot = read.csv(file="~/share/reference/IMPACT410_genes_for_copynumber.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
  			dplyr::select(hgnc_symbol, chr, start_position, end_position) %>%
  			rename(Hugo_Symbol = hgnc_symbol,
  				   Chromosome = chr,
  				   Start = start_position,
  				   End = end_position) %>%
  			mutate(Chromosome = ifelse(Chromosome %in% "X", 23, Chromosome)) %>%
  			arrange(as.numeric(Chromosome), as.numeric(Start))
  				   
  	annot_by_gene <- annot %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Hugo_Symbol = Hugo_Symbol)
  	res_by_segment <- res %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Calls = Calls)
  	fo <- findOverlaps(res_by_segment, annot_by_gene)
  
  	df <- data.frame(Hugo_Symbol=mcols(annot_by_gene)[subjectHits(fo),], Calls=mcols(res_by_segment)[queryHits(fo),])
  	Hugo_Symbol = which(duplicated(df$Hugo_Symbol))
  	for (j in 1:length(Hugo_Symbol)) {
  		index = which(as.character(df$Hugo_Symbol)==as.character(df$Hugo_Symbol[Hugo_Symbol[j]]))
  		df[index,2] = round(mean(df[index,2], na.rm=TRUE))
  	}
  	df = df %>% filter(!duplicated(Hugo_Symbol))
  	res = rep(0, nrow(annot))
  	names(res) = annot[,"Hugo_Symbol"]
  	res[as.character(df$Hugo_Symbol)] = df$Calls
  	return(invisible(res))
}
g_cn = do.call(cbind, g_cn)
colnames(g_cn) = tracker$GRAIL_ID
 
featureNames = intersect(rownames(i_cn), rownames(g_cn))
i_cn = i_cn[featureNames,,drop=FALSE]
g_cn = g_cn[featureNames,,drop=FALSE]

#==================================================
# ROC curves amplification
#==================================================
data = foreach (i=1:nrow(tracker)) %dopar% {
	x = g_cn[,i]
	y = ifelse(i_cn[,i]<=0, 0, 1)
	data = data.frame(x=x, y=y)
	return(invisible(data))
}
data = do.call(rbind, data)
model.fit = pROC::roc(y ~ x, data=data, smooth=FALSE)
show.0 = data.frame(x = 1-model.fit$specificities,
					y = model.fit$sensitivities) %>%
					mutate(sample_type = "Combined")

index = grep("VB", tracker$GRAIL_ID)
data = foreach (i=1:length(index)) %dopar% {
	x = g_cn[,index[i]]
	y = ifelse(i_cn[,index[i]]<=0, 0, 1)
	data = data.frame(x=x, y=y)
	return(invisible(data))
}
data = do.call(rbind, data)
model.fit = pROC::roc(y ~ x, data=data, smooth=FALSE)
show.1 = data.frame(x = 1-model.fit$specificities,
					y = model.fit$sensitivities) %>%
					mutate(sample_type = "Breast")
					   
index = grep("VL", tracker$GRAIL_ID)
data = foreach (i=1:length(index)) %dopar% {
	x = g_cn[,index[i]]
	y = ifelse(i_cn[,index[i]]<=0, 0, 1)
	data = data.frame(x=x, y=y)
	return(invisible(data))
}
data = do.call(rbind, data)
model.fit = pROC::roc(y ~ x, data=data, smooth=FALSE)
show.2 = data.frame(x = 1-model.fit$specificities,
					y = model.fit$sensitivities) %>%
					mutate(sample_type = "Lung")

index = grep("VP", tracker$GRAIL_ID)
data = foreach (i=1:length(index)) %dopar% {
	x = g_cn[,index[i]]
	y = ifelse(i_cn[,index[i]]<=0, 0, 1)
	data = data.frame(x=x, y=y)
	return(invisible(data))
}
data = do.call(rbind, data)
model.fit = pROC::roc(y ~ x, data=data, smooth=FALSE)
show.3 = data.frame(x = 1-model.fit$specificities,
					y = model.fit$sensitivities) %>%
					mutate(sample_type = "Prostate")

show.data = rbind(show.0, show.1, show.2, show.3) %>%
			mutate(sample_type = factor(sample_type, levels = c("Combined", "Breast", "Lung", "Prostate"), ordered=TRUE)) %>%
			mutate(copy_type = "Copy number amplifications")

plot.0 = ggplot(show.data, aes(x = x, y = y, group = sample_type, color = sample_type)) +
		 geom_line(size=1) +
		 scale_color_manual(values = c("Breast"="salmon", "Lung"="#FDAE61", "Prostate"="#ABDDA4", "Combined"="cadetblue")) +
		 theme_bw(base_size=15) +
		 facet_wrap(~copy_type) +
		 theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=12)) +
		 labs(x="\n1 - Specificity", y="Sensitivity\n") +
         coord_cartesian(xlim=c(0, 1), ylim = c(0, 1)) +
         theme(legend.justification = c(1, 0),
		 	   legend.position = c(1, 0),
		 	   legend.title = element_blank(),
		 	   legend.background = element_blank(),
		 	   legend.text=element_text(size=8))
         		
pdf(file="../res/rebuttal/Analytical_performance_copy_amplifications.pdf", width=5, height=5)
print(plot.0)
dev.off()

#==================================================
# ROC curves deletion
#==================================================
data = foreach (i=1:nrow(tracker)) %dopar% {
	x = g_cn[,i]
	y = ifelse(i_cn[,i]<0, -1, 0)
	data = data.frame(x=x, y=y)
	return(invisible(data))
}
data = do.call(rbind, data)
model.fit = pROC::roc(y ~ x, data=data, smooth=FALSE)
show.0 = data.frame(x = 1-model.fit$specificities,
					y = model.fit$sensitivities) %>%
					mutate(sample_type = "Combined") %>%
					arrange(x, y)

index = grep("VB", tracker$GRAIL_ID)
data = foreach (i=1:length(index)) %dopar% {
	x = g_cn[,index[i]]
	y = ifelse(i_cn[,index[i]]<0, -1, 0)
	data = data.frame(x=x, y=y)
	return(invisible(data))
}
data = do.call(rbind, data)
model.fit = pROC::roc(y ~ x, data=data, smooth=FALSE)
show.1 = data.frame(x = 1-model.fit$specificities,
					y = model.fit$sensitivities) %>%
					mutate(sample_type = "Breast") %>%
					arrange(x, y)
					   
index = grep("VL", tracker$GRAIL_ID)
data = foreach (i=1:length(index)) %dopar% {
	x = g_cn[,index[i]]
	y = ifelse(i_cn[,index[i]]<0, -1, 0)
	data = data.frame(x=x, y=y)
	return(invisible(data))
}
data = do.call(rbind, data)
model.fit = pROC::roc(y ~ x, data=data, smooth=FALSE)
show.2 = data.frame(x = 1-model.fit$specificities,
					y = model.fit$sensitivities) %>%
					mutate(sample_type = "Lung") %>%
					arrange(x, y)

index = grep("VP", tracker$GRAIL_ID)
data = foreach (i=1:length(index)) %dopar% {
	x = g_cn[,index[i]]
	y = ifelse(i_cn[,index[i]]<0, -1, 0)
	data = data.frame(x=x, y=y)
	return(invisible(data))
}
data = do.call(rbind, data)
model.fit = pROC::roc(y ~ x, data=data, smooth=FALSE)
show.3 = data.frame(x = 1-model.fit$specificities,
					y = model.fit$sensitivities) %>%
					mutate(sample_type = "Prostate") %>%
					arrange(x, y)

show.data = rbind(show.0, show.1, show.2, show.3) %>%
			mutate(sample_type = factor(sample_type, levels = c("Combined", "Breast", "Lung", "Prostate"), ordered=TRUE)) %>%
			mutate(copy_type = "Homozygous deletions")

plot.0 = ggplot(show.data, aes(x = x, y = y, group = sample_type, color = sample_type)) +
		 geom_line(size=1) +
		 scale_color_manual(values = c("Breast"="salmon", "Lung"="#FDAE61", "Prostate"="#ABDDA4", "Combined"="cadetblue")) +
		 theme_bw(base_size=15) +
		 facet_wrap(~copy_type) +
		 theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=12)) +
		 labs(x="\n1 - Specificity", y="Sensitivity\n") +
         coord_cartesian(xlim=c(0, 1), ylim = c(0, 1)) +
         theme(legend.justification = c(1, 0),
		 	   legend.position = c(1, 0),
		 	   legend.title = element_blank(),
		 	   legend.background = element_blank(),
		 	   legend.text=element_text(size=8))
         		
pdf(file="../res/rebuttal/Analytical_performance_copy_deletions.pdf", width=5, height=5)
print(plot.0)
dev.off()

#==================================================
# Probit curves
#==================================================
tracker = read_tsv(file=url_master_key, col_types = cols(.default = col_character())) %>%
		  type_convert() %>%
 		  dplyr::select(PATIENT_ID, GRAIL_ID, DMP_ID, TUMOR_ID, NORMAL_ID, GRAIL_alpha, GRAIL_psi, IMPACT_alpha, IMPACT_psi)
		   
ctDNA_fraction = read_csv(file=url_ctdna_frac, col_types = cols(.default = col_character())) %>%
		   		 type_convert() %>%
		   		 rename(GRAIL_ID = ID)
		   		 
tracker = left_join(tracker, ctDNA_fraction, by="GRAIL_ID") %>%
 		  filter(!is.na(ctdna_frac)) %>%
		  filter(ctdna_frac>.1)
 
sse = foreach (i=1:nrow(tracker)) %dopar% {
  	cat(tracker$GRAIL_ID[i], "\n")
  	if (!is.na(tracker$ctdna_frac[i]) & tracker$ctdna_frac[i]>.1) {
 	  	grail_path = paste0("../res/rebuttal/uncollapsed_bam/cnvkit/totalcopy/", tracker$GRAIL_ID[i], "-T.RData")
 	 	grail_data = new.env()
 	 	load(grail_path, envir=grail_data)
  	 	grail_seg = grail_data$tmp %>%
  				    dplyr::select(chrom=Chromosome, start = Start, end = End, log2 = Log2Ratio, n=N) %>%
  				    filter(chrom<23)
  		grail_seg = prune_(x=grail_seg)
  		sse = ploidy_(x=grail_seg$log2, y=tracker$ctdna_frac[i], w=grail_seg$n)
  	} else {
  		sse = rep(NA, 100)
  	}
  	return(sse)
}
min_sse = unlist(lapply(sse, min))
index = rep(NA, length(sse))
index[!is.na(min_sse)] = unlist(lapply(sse, which.min))
ploidy = seq(from=1.5, to=3.5, length=100)[index]
tracker$GRAIL_psi = ploidy
tracker = tracker %>%
	   	  filter(!is.na(GRAIL_psi) & !is.na(ctdna_frac))

i_cn = foreach (i=1:nrow(tracker)) %dopar% {
 	cat(tracker$GRAIL_ID[i], "\n")
  	path = paste0("../res/rebuttal/msk_impact/cnvkit/totalcopy/", tracker$TUMOR_ID[i], ".RData")
 	data = new.env()
 	load(path, envir=data)
  	seg = data$tmp %>%
  		  dplyr::select(chrom=Chromosome, start = Start, end = End, log2 = Log2Ratio, n=N) %>%
  		  filter(chrom<23)
  	seg = prune_(x=seg)
  	a_cn = round(absolute_(rho=tracker$IMPACT_alpha[i], psi=tracker$GRAIL_psi[i], gamma=.85, seg$log2))
   	
 	seg = cbind(seg, a_cn)
  	
  	Chromosome = seg[,"chrom"]
  	Start = seg[,"start"]
  	End = seg[,"end"]
  	Calls = seg[,"a_cn"]
  	res = data.frame(Chromosome, Start, End, Calls)
  	annot = read.csv(file="~/share/reference/IMPACT410_genes_for_copynumber.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
  			dplyr::select(hgnc_symbol, chr, start_position, end_position) %>%
  			rename(Hugo_Symbol = hgnc_symbol,
  				   Chromosome = chr,
  				   Start = start_position,
  				   End = end_position) %>%
  			mutate(Chromosome = ifelse(Chromosome %in% "X", 23, Chromosome)) %>%
  			arrange(as.numeric(Chromosome), as.numeric(Start))
  				   
  	annot_by_gene <- annot %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Hugo_Symbol = Hugo_Symbol)
  	res_by_segment <- res %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Calls = Calls)
  	fo <- findOverlaps(res_by_segment, annot_by_gene)
  
  	df <- data.frame(Hugo_Symbol=mcols(annot_by_gene)[subjectHits(fo),], Calls=mcols(res_by_segment)[queryHits(fo),])
  	Hugo_Symbol = which(duplicated(df$Hugo_Symbol))
  	for (j in 1:length(Hugo_Symbol)) {
  		index = which(as.character(df$Hugo_Symbol)==as.character(df$Hugo_Symbol[Hugo_Symbol[j]]))
  		df[index,2] = round(max(df[index,2], na.rm=TRUE))
  	}
  	df = df %>% filter(!duplicated(Hugo_Symbol))
  	df[,2] = round(df[,2])
  	res = rep(0, nrow(annot))
  	names(res) = annot[,"Hugo_Symbol"]
  	res[as.character(df$Hugo_Symbol)] = df$Calls
  	return(invisible(res))
}
i_cn = do.call(cbind, i_cn)
colnames(i_cn) = tracker$GRAIL_ID
 
g_cn = foreach (i=1:nrow(tracker)) %dopar% {
	cat(tracker$GRAIL_ID[i], "\n")
 	path = paste0("../res/rebuttal/uncollapsed_bam/cnvkit/totalcopy/", tracker$GRAIL_ID[i], "-T.RData")
 	data = new.env()
 	load(path, envir=data)
  	seg = data$tmp %>%
  		  dplyr::select(chrom=Chromosome, start = Start, end = End, log2 = Log2Ratio, n=N) %>%
  		  filter(chrom<23)
  	seg = prune_(x=seg)
  	a_cn = round(absolute_(rho=tracker$ctdna_frac[i], psi=tracker$GRAIL_psi[i], gamma=.85, seg$log2))
 	c_cn = rep(0, length(a_cn))
  	if (round(tracker$GRAIL_psi[i])==2) {
  		c_cn[a_cn<0.5] = -1
  		c_cn[a_cn>6] = 1
  	} else if (round(tracker$GRAIL_psi[i])==3) {
  		c_cn[a_cn<1] = -1
  		c_cn[a_cn>7] = 1
  	} else if (round(tracker$IMPACT_psi[i])==4) {
  		c_cn[a_cn<1] = -1
   		c_cn[a_cn>8] = 1
  	}
  	
  	sample_names = c("MSK-VB-0024", "MSK-VP-0028", "MSK-VB-0032", "MSK-VB-0046", "MSK-VP-0031")
 
  	if (tracker$GRAIL_ID[i] %in% sample_names) {
  		c_cn = rep(0, length(a_cn))
  		c_cn[a_cn<0] = -1
   		c_cn[a_cn>7] = 1
   	}
  	
  	seg = cbind(seg, a_cn, c_cn)
  	
  	Chromosome = seg[,"chrom"]
  	Start = seg[,"start"]
  	End = seg[,"end"]
  	Calls = seg[,"c_cn"]
  	res = data.frame(Chromosome, Start, End, Calls)
  	annot = read.csv(file="~/share/reference/IMPACT410_genes_for_copynumber.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
  			dplyr::select(hgnc_symbol, chr, start_position, end_position) %>%
  			rename(Hugo_Symbol = hgnc_symbol,
  				   Chromosome = chr,
  				   Start = start_position,
  				   End = end_position) %>%
  			mutate(Chromosome = ifelse(Chromosome %in% "X", 23, Chromosome)) %>%
  			arrange(as.numeric(Chromosome), as.numeric(Start))
  				   
  	annot_by_gene <- annot %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Hugo_Symbol = Hugo_Symbol)
  	res_by_segment <- res %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Calls = Calls)
  	fo <- findOverlaps(res_by_segment, annot_by_gene)
  
  	df <- data.frame(Hugo_Symbol=mcols(annot_by_gene)[subjectHits(fo),], Calls=mcols(res_by_segment)[queryHits(fo),])
  	Hugo_Symbol = which(duplicated(df$Hugo_Symbol))
  	for (j in 1:length(Hugo_Symbol)) {
  		index = which(as.character(df$Hugo_Symbol)==as.character(df$Hugo_Symbol[Hugo_Symbol[j]]))
  		df[index,2] = round(mean(df[index,2], na.rm=TRUE))
  	}
  	df = df %>% filter(!duplicated(Hugo_Symbol))
  	res = rep(0, nrow(annot))
  	names(res) = annot[,"Hugo_Symbol"]
  	res[as.character(df$Hugo_Symbol)] = df$Calls
  	return(invisible(res))
}
g_cn = do.call(cbind, g_cn)
colnames(g_cn) = tracker$GRAIL_ID
 
featureNames = intersect(rownames(i_cn), rownames(g_cn))
i_cn = i_cn[featureNames,,drop=FALSE]
g_cn = g_cn[featureNames,,drop=FALSE]

#==================================================
# Probit curves amplification
#==================================================
data = foreach (i=1:nrow(tracker)) %dopar% {
	x = i_cn[,i]
	y = ifelse(g_cn[,i]>0, 1, 0)
	data = data.frame(x=x, y=y)
	return(invisible(data))
}
data = do.call(rbind, data)
model = glm(y ~ x, data=data, family = binomial(link = "probit"))
new_data = data.frame(x = seq(0, 50, 0.01))
predicted_data = as.data.frame(predict(model, newdata = new_data, se = TRUE))
show_data.0 = cbind(new_data, predicted_data) %>%
			  mutate(ymin = model$family$linkinv(fit - 1.96*se.fit),
				     ymax = model$family$linkinv(fit + 1.96*se.fit),
				     yfit = model$family$linkinv(fit),
				     grp = "Combined")
				     
index = grep("VB", tracker$GRAIL_ID)
data = foreach (i=1:length(index)) %dopar% {
	x = i_cn[,index[i]]
	y = ifelse(g_cn[,index[i]]>0, 1, 0)
	data = data.frame(x=x, y=y)
	return(invisible(data))
}
data = do.call(rbind, data)
model = glm(y ~ x, data=data, family = binomial(link = "probit"))
new_data = data.frame(x = seq(0, 50, 0.01))
predicted_data = as.data.frame(predict(model, newdata = new_data, se = TRUE))
show_data.1 = cbind(new_data, predicted_data) %>%
			  mutate(ymin = model$family$linkinv(fit - 1.96*se.fit),
				     ymax = model$family$linkinv(fit + 1.96*se.fit),
				     yfit = model$family$linkinv(fit),
				     grp = "Breast")
				     
index = grep("VL", tracker$GRAIL_ID)
data = foreach (i=1:length(index)) %dopar% {
	x = i_cn[,index[i]]
	y = ifelse(g_cn[,index[i]]>0, 1, 0)
	data = data.frame(x=x, y=y)
	return(invisible(data))
}
data = do.call(rbind, data)
model = glm(y ~ x, data=data, family = binomial(link = "probit"))
new_data = data.frame(x = seq(0, 50, 0.01))
predicted_data = as.data.frame(predict(model, newdata = new_data, se = TRUE))
show_data.2 = cbind(new_data, predicted_data) %>%
			  mutate(ymin = model$family$linkinv(fit - 1.96*se.fit),
				     ymax = model$family$linkinv(fit + 1.96*se.fit),
				     yfit = model$family$linkinv(fit),
				     grp = "Lung")
				     
index = grep("VP", tracker$GRAIL_ID)
data = foreach (i=1:length(index)) %dopar% {
	x = i_cn[,index[i]]
	y = ifelse(g_cn[,index[i]]>0, 1, 0)
	data = data.frame(x=x, y=y)
	return(invisible(data))
}
data = do.call(rbind, data)
model = glm(y ~ x, data=data, family = binomial(link = "probit"))
new_data = data.frame(x = seq(0, 50, 0.01))
predicted_data = as.data.frame(predict(model, newdata = new_data, se = TRUE))
show_data.3 = cbind(new_data, predicted_data) %>%
			  mutate(ymin = model$family$linkinv(fit - 1.96*se.fit),
				     ymax = model$family$linkinv(fit + 1.96*se.fit),
				     yfit = model$family$linkinv(fit),
				     grp = "Prostate")
				     
show_data = bind_rows(show_data.0, show_data.1, show_data.2, show_data.3) %>%
			mutate(grp = factor(grp, levels = c("Combined", "Breast", "Lung", "Prostate"), ordered=TRUE))

plot.0 = ggplot(show_data, aes(x = x, y = yfit, group=grp, color=grp)) +
		 geom_line(size=1) +
		 geom_line(size=.5, linetype=1, aes(x=x, y=ymin)) +
		 geom_line(size=.5, linetype=1, aes(x=x, y=ymax)) +
		 scale_color_manual(values = c("Breast"="salmon", "Lung"="#FDAE61", "Prostate"="#ABDDA4", "Combined"="cadetblue")) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), legend.position="none") +
		 facet_wrap(~grp, nrow=2, ncol=2) +
		 labs(x="\nCopy number", y="Detection probability\n") +
         coord_cartesian(xlim=c(0, 25), ylim = c(0, 1))
         		
pdf(file="../res/rebuttal/Analytical_performance_copy_amplifications_Probit.pdf", width=7, height=5.5)
print(plot.0)
dev.off()

# copy = NULL
# for (i in 1:ncol(i_cn)) {
# 	copy = cbind(copy, i_cn[,i], g_cn[,i])
# }
# index = apply(copy, 1, function(x) {sum(x==0)})==98
# copy = copy[!index,,drop=FALSE]
# assay_type = rep(c("Biopsy", "cfDNA"), times=nrow(key_file))
# cancer_type = rep("Breast", times=nrow(key_file))
# cancer_type[grepl("VL", key_file$GRAIL_ID)] = "Lung"
# cancer_type[grepl("VP", key_file$GRAIL_ID)] = "Prostate"
# cancer_type = rep(cancer_type, each=2)
# 
# pdf(file="../res/rebuttal/Heatmap_CN_Tumor_cfDNA.pdf", width=13, height=7/10*12)
# mm = split.screen(fig=matrix(c(0.02,.32,0,1,	0,.25,0,.98,	.16,1,0,1), nrow=3, ncol=4, byrow=TRUE))
# screen(mm[1])
# plot(0, 0, type="n", xlab="", ylab="", xlim=c(0,5), ylim=c(.5,ncol(copy)-.5), axes=FALSE, frame.plot=FALSE)
# cols = c("Breast"="salmon", "Lung"="#FDAE61", "Prostate"="#ABDDA4")[cancer_type]
# z = 3.8
# for (i in 1:length(cols)) {
# 	rect(xleft=z, xright=z+.2, ybottom=i-.5, ytop=i+.5, col=cols[i], border="white", lwd=.5)
# }
# cols = ifelse(assay_type=="Biopsy", "#AB6E9A", "#7B1E5B")
# for (i in 1:length(cols)) {
# 	rect(xleft=z+.23, xright=z+.23+.2, ybottom=i-.5, ytop=i+.5, col=cols[i], border="white", lwd=.5)
# }
# screen(mm[2])
# z = .5
# plot(0, 0, type="n", xlab="", ylab="", xlim=c(0,5), ylim=c(0,ncol(copy)), axes=FALSE, frame.plot=FALSE)
# points(0+z, 100, type="p", pch=22, bg="#CF3A3D", col="black", cex=1.5)
# text(x=0+z, y=100, labels="Amplification", pos=4, cex=.78)
# points(0+z, 97, type="p", pch=22, bg="#2A4B94", col="black", cex=1.5)
# text(x=0+z, y=97, labels="Homozygous deletion", pos=4, cex=.78)
# points(0+z, 90, type="p", pch=22, bg="salmon", col="black", cex=1.5)
# text(x=0+z, y=90, labels="Breast", pos=4, cex=.78)
# points(0+z, 87, type="p", pch=22, bg="#FDAE61", col="black", cex=1.5)
# text(x=0+z, y=87, labels="Lung", pos=4, cex=.78)
# points(0+z, 84, type="p", pch=22, bg="#ABDDA4", col="black", cex=1.5)
# text(x=0+z, y=84, labels="Prostate", pos=4, cex=.78)
# points(0+z, 77, type="p", pch=22, bg="#AB6E9A", col="black", cex=1.5)
# text(x=0+z, y=77, labels="Biopsy", pos=4, cex=.78)
# points(0+z, 74, type="p", pch=22, bg="#7B1E5B", col="black", cex=1.5)
# text(x=0+z, y=74, labels="cfDNA", pos=4, cex=.78)
# screen(mm[3])
# cols = c("#2A4B94", NA, "#CF3A3D")
# plot(0, 0, type="n", xlab="", ylab="", xlim=c(0,nrow(copy)), ylim=c(0,ncol(copy)-1), axes=FALSE, frame.plot=FALSE)
# for (i in 1:(ncol(copy))) {
#  	col = cols[as.numeric(copy[,i])+2]
#  	for (j in 1:length(col)) {
#  		rect(xleft=j-.5, xright=j+.5, ybottom=i-1, ytop=i, col=col[j], border="grey95", lwd=.05)
#  	}
# }
# rect(xleft=0.5, xright=nrow(copy)+.5, ybottom=0, ytop=ncol(copy), border="grey90")
# close.screen(all.screens=TRUE)
# dev.off()

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
# Breast ERBB2 amplified cases
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
	plot(1, 1, type="n", xlim=c(start,end), ylim=c(-4,4.75), xlab="", ylab="", main="", frame.plot=FALSE, axes=FALSE)
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
	axis(2, at = c(-4,-2,0,2,4), labels=c(-4,-2,0,2,4), cex.axis = 1.5, las = 1, lwd=1.5, lwd.ticks=1.35, col="grey20", col.axis="grey20")
	axis(2, at = c(-3,-1,1,3), labels=rep("",4), cex.axis = 1.5, las = 1, lwd=0, lwd.ticks=1.15, col="grey20", col.axis="grey20")
	mtext(side = 2, text = expression("Log"[2]~"Ratio"), line = 4, cex = 1.5)
	rect(xleft=1-1e10, xright=max(z1[,"pos"])+1e10, ybottom=4, ytop=6, col="lightgrey", border="grey20", lwd=2.5)
	box(lwd=2.5, col="grey20")
	screen(zz[2])	
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
	axis(2, at = c(-4,-2,0,2,4), labels=c(-4,-2,0,2,4), cex.axis = 1.5, las = 1, lwd=1.5, lwd.ticks=1.35, col="grey20", col.axis="grey20")
	axis(2, at = c(-3,-1,1,3), labels=rep("",4), cex.axis = 1.5, las = 1, lwd=0, lwd.ticks=1.15, col="grey20", col.axis="grey20")
	mtext(side = 2, text = expression("Log"[2]~"Ratio"), line = 4, cex = 1.5)
	box(lwd=2.5, col="grey20")
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
	plot(1, 1, type="n", xlim=c(start,end), ylim=c(-4,4.75), xlab="", ylab="", main="", frame.plot=FALSE, axes=FALSE)
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
	axis(2, at = c(-4,-2,0,2,4), labels=c(-4,-2,0,2,4), cex.axis = 1.5, las = 1, lwd=1.5, lwd.ticks=1.35, col="grey20", col.axis="grey20")
	axis(2, at = c(-3,-1,1,3), labels=rep("",4), cex.axis = 1.5, las = 1, lwd=0, lwd.ticks=1.15, col="grey20", col.axis="grey20")
	mtext(side = 2, text = expression("Log"[2]~"Ratio"), line = 4, cex = 1.5)
	rect(xleft=1-1e10, xright=max(z1[,"pos"])+1e10, ybottom=4, ytop=6, col="lightgrey", border="grey20", lwd=2.5)
	box(lwd=2.5, col="grey20")
	screen(zz[2])
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
	axis(2, at = c(-4,-2,0,2,4), labels=c(-4,-2,0,2,4), cex.axis = 1.5, las = 1, lwd=1.5, lwd.ticks=1.35, col="grey20", col.axis="grey20")
	axis(2, at = c(-3,-1,1,3), labels=rep("",4), cex.axis = 1.5, las = 1, lwd=0, lwd.ticks=1.15, col="grey20", col.axis="grey20")
	mtext(side = 2, text = expression("Log"[2]~"Ratio"), line = 4, cex = 1.5)
	box(lwd=2.5, col="grey20")
	screen(zz[2])
	assembly = read.csv(file=url_cytogenetic_data, header=FALSE, sep="\t", stringsAsFactors=FALSE)
	plotIdeogram(chrom=7, TRUE, cyto.data = assembly, cex = .75, unit = "bp")
	close.screen(all.screens=TRUE)
	dev.off()
}
