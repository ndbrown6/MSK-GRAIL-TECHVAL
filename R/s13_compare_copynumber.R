#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/rebuttal")) {
	dir.create("../res/rebuttal")
}

'prune_segments_cn' <- function(x, n=10)
{
		cnm = matrix(NA, nrow=nrow(x), ncol=nrow(x))
		for (j in 1:nrow(x)) {
			cnm[,j] = abs(2^x[j,"cnlr.median"] - 2^x[,"cnlr.median"])
		}
		cnt = hclust(as.dist(cnm), "average")
		cnc = cutree(tree=cnt, k=n)
		for (j in unique(cnc)) {
			indx = which(cnc==j)
			if (length(indx)>2) {
				mcl = mean(x[indx,"cnlr.median"])
				scl = sd(x[indx,"cnlr.median"])
				ind = which(x[indx,"cnlr.median"]<(mcl+1.96*scl) & x[indx,"cnlr.median"]>(mcl-1.96*scl))
				x[indx[ind],"cnlr.median"] = mean(x[indx[ind],"cnlr.median"])
			} else {
				x[indx,"cnlr.median"] = mean(x[indx,"cnlr.median"])
			}
		}
		return(x)
}

'plot_log2_ratio' <- function(x, y, xlab=FALSE, n=10, facet)
{
    segmented = x$jointseg
    chr = segmented$chrom
    len = table(chr)
    tmp = cumsum(len)
    start = c(1, tmp[-length(len)] + 1)
    end = tmp
    mid = start + len/2
    
    purity = y$purity
    ploidy = y$ploidy
    ki = ky = c(1, 3, 7, 13, 21)
    for (i in 1:length(ki)) {
   		ky[i] = (log2(((purity)*ki[i] + (1-purity)*2)/((purity)*ploidy + (1-purity)*2)))
   	}
   	
    
    
    tmp = bind_cols(segmented %>%
    	  				select(chrom, cnlr) %>%
    	  				rename(Chromosome = chrom, Log2Ratio = cnlr) %>%
    	  				mutate(Residual = as.factor(Chromosome %% 2)) %>%
    	  				mutate(Facet = facet),
    	  			Index = 1:nrow(segmented))
    	  
    	  
    tmp[,"Log2Ratio"] = winsorize(segmented[,c("chrom", "maploc", "cnlr")], method="pcf", gamma=100, verbose = FALSE)[,3]
    
    if (xlab) {
	    plot.0 = ggplot(tmp, aes(x=Index, y=Log2Ratio, fill=Residual, color=Residual)) +
  		 	 geom_point(alpha=1, size=.3, shape=20) +
			 geom_hline(yintercept=ky, color="goldenrod3", linetype=3, size=.25) +
  		 	 geom_vline(xintercept=mid, size=.25) +
  		 	 scale_fill_manual(guide = FALSE, values = c("grey75", "cadetblue")) +
  		 	 scale_color_manual(guide = FALSE, values = c("grey75", "cadetblue")) +
  		 	 facet_wrap(~Facet) +
		 	 theme_bw(base_size=9) +
		 	 labs(x="", y=expression(Log[2]~"Ratio")) +
		 	 coord_cartesian(ylim = c(-4,4)) +
 		 	 theme(axis.text.y = element_text(size=10),
 		 	 	   axis.text.x = element_blank(),
 		 	 	   axis.ticks.x = element_blank(),
 		 	 	   panel.grid.major = element_blank(),
 		 	 	   panel.grid.minor = element_blank(),
 		 	 	   legend.justification = c(1, 0),
		 	   	   legend.position = c(1, 0),
		 	   	   legend.title = element_blank(),
		 	   	   legend.background = element_blank(),
		 	   	   legend.text=element_text(size=8))
	} else {
		plot.0 = ggplot(tmp, aes(x=Index, y=Log2Ratio, fill=Residual, color=Residual)) +
  		 	 geom_point(alpha=1, size=.3, shape=20) +
			 geom_hline(yintercept=ky, color="goldenrod3", linetype=3, size=.25) +
  		 	 scale_fill_manual(guide = FALSE, values = c("grey75", "cadetblue")) +
  		 	 scale_color_manual(guide = FALSE, values = c("grey75", "cadetblue")) +
  		 	 facet_wrap(~Facet) +
		 	 theme_bw(base_size=9) +
		 	 labs(x="", y=expression(Log[2]~"Ratio")) +
		 	 coord_cartesian(ylim = c(-4,4)) +
 		 	 theme(axis.text.y = element_text(size=9),
 		 	 	   axis.text.x = element_blank(),
 		 	 	   axis.ticks.x = element_blank(),
 		 	 	   panel.grid.major = element_blank(),
 		 	 	   panel.grid.minor = element_blank(),
 		 	 	   legend.justification = c(1, 0),
		 	   	   legend.position = c(1, 0),
		 	   	   legend.title = element_blank(),
		 	   	   legend.background = element_blank(),
		 	   	   legend.text=element_text(size=8))
	}
	print(plot.0)
}


#==================================================
# Examples of GRAIL cfDNA
#==================================================
load("../res/rebuttal/GRAIL/facets/cncf/MSK-VB-0069_MSK-VB-0069-N.Rdata")
fit$ploidy = 2
fit$purity = 1
pdf(file="../res/rebuttal/MSK-VB-0069_MSK-VB-0069-N.pdf", width=7, height=2)
plot_log2_ratio(x=out2, y=fit, xlab=FALSE, n=4, facet="MSK-VB-0069 | cfDNA")
dev.off()

load("../res/rebuttal/GRAIL/facets/cncf/MSK-VL-0056_MSK-VL-0056-N.Rdata")
fit$ploidy = 2
fit$purity = 1
pdf(file="../res/rebuttal/MSK-VL-0056_MSK-VL-0056-N.pdf", width=7, height=2)
plot_log2_ratio(x=out2, y=fit, xlab=FALSE, n=5, facet="MSK-VL-0056 | cfDNA")
dev.off()

load("../res/rebuttal/GRAIL/facets/cncf/MSK-VP-0004_MSK-VP-0004-N.Rdata")
fit$ploidy = 2
fit$purity = 1
pdf(file="../res/rebuttal/MSK-VP-0004_MSK-VP-0004-N.pdf", width=7, height=2)
plot_log2_ratio(x=out2, y=fit, xlab=FALSE, n=7, facet="MSK-VP-0004 | cfDNA")
dev.off()

load("../res/rebuttal/GRAIL/facets/cncf/MSK-VB-0023_MSK-VB-0023-N.Rdata")
fit$ploidy = 2
fit$purity = 1
pdf(file="../res/rebuttal/MSK-VB-0023_MSK-VB-0023-N.pdf", width=7, height=2)
plot_log2_ratio(x=out2, y=fit, xlab=FALSE, n=7, facet="MSK-VB-0023 | cfDNA")
dev.off()

#==================================================
# Examples of MSK-IMPACT tumor biopsies
#==================================================
load("../res/rebuttal/MSK-IMPACT/facets/cncf/UE782797-T_TG260267-N.Rdata")
fit$ploidy = 2
fit$purity = 1
pdf(file="../res/rebuttal/UE782797-T_TG260267-N.pdf", width=7, height=2)
plot_log2_ratio(x=out2, y=fit, xlab=FALSE, n=6, facet="MSK-VB-0069 | Tumor biopsy")
dev.off()

load("../res/rebuttal/MSK-IMPACT/facets/cncf/GK457095-T_GK457095-N.Rdata")
fit$ploidy = 2
fit$purity = 1
pdf(file="../res/rebuttal/GK457095-T_GK457095-N.pdf", width=7, height=2)
plot_log2_ratio(x=out2, y=fit, xlab=FALSE, n=9, facet="MSK-VL-0056 | Tumor biopsy")
dev.off()

load("../res/rebuttal/MSK-IMPACT/facets/cncf/GV175163-T_GV175163-N.Rdata")
fit$ploidy = 2
fit$purity = 1
pdf(file="../res/rebuttal/GV175163-T_GV175163-N.pdf", width=7, height=2)
plot_log2_ratio(x=out2, y=fit, xlab=FALSE, n=9, facet="MSK-VP-0004 | Tumor biopsy")
dev.off()

load("../res/rebuttal/MSK-IMPACT/facets/cncf/DW801941-T_DW801941-N.Rdata")
fit$ploidy = 2
fit$purity = 1
pdf(file="../res/rebuttal/DW801941-T_DW801941-N.pdf", width=7, height=2)
plot_log2_ratio(x=out2, y=fit, xlab=FALSE, n=9, facet="MSK-VB-0023 | Tumor biopsy")
dev.off()


# ==================================================
#  Comparison of copy number aberrations
# ==================================================
# 'absolute_copies' <- function(rho, psi, gamma=1, x) {
# 	return(invisible(((((2^(x/gamma))*(rho*psi+(1-rho)*2)) - ((1-rho)*2))/rho)))
# }
# 
# tracker = read_tsv(file="../res/rebuttal/GRAIL/dev/data/master_sample_key.tsv", col_types = cols(.default = col_character()))  %>%
#  		  type_convert()
# 
# i_cn = foreach (i=1:nrow(tracker)) %dopar% {
# 	cat(i, "\n")
# 	load(paste0("../res/rebuttal/MSK-IMPACT/facets/cncf/", tracker$TUMOR_ID[i], "_", tracker$NORMAL_ID[i], ".Rdata"))
# 	fit$purity = ifelse(is.na(fit$purity), 0, fit$purity)
# 	fit$ploidy = ifelse(is.na(fit$ploidy), 2, fit$ploidy)
# 	fit$cncf$tcn.em = round(absolute_copies(fit$purity, fit$ploidy, gamma=1, x=fit$cncf$cnlr.median))
# 	fit$cncf$called = rep(0, nrow(fit$cncf))
# 
# 	if (fit$purity>=.5) {
# 		fit$cncf$called = ifelse(fit$cncf$cnlr.median>(.55), 2, 0)
#  		fit$cncf$called = ifelse(fit$cncf$cnlr.median<(-.65), -2, fit$cncf$called)
# 	} else if (fit$purity>=.2 & fit$purity<.5) {
# 		fit$cncf$called = ifelse(fit$cncf$cnlr.median>(.45), 2, 0)
#  		fit$cncf$called = ifelse(fit$cncf$cnlr.median<(-.50), -2, fit$cncf$called)
#     } else {
#      	fit$cncf$called = ifelse(fit$cncf$cnlr.median>(.35), 2, 0)
#      	fit$cncf$called = ifelse(fit$cncf$cnlr.median<(-.30), -2, fit$cncf$called)
#  	}
# 
# 	fit$cncf$called = ifelse(fit$cncf$cnlr.median>(median(fit$cncf$cnlr.median) + 2*sd(fit$cncf$cnlr.median)), 2, fit$cncf$called)
# 	fit$cncf$called = ifelse(fit$cncf$cnlr.median<(median(fit$cncf$cnlr.median) - 2*sd(fit$cncf$cnlr.median)), -2, fit$cncf$called)
# 
# 	Chromosome = fit$cncf[,"chrom"]
# 	Start = fit$cncf[,"start"]
# 	End = fit$cncf[,"end"]
# 	Calls = fit$cncf$called
# 	res = data.frame(Chromosome, Start, End, Calls)
# 	annot = read.csv(file="~/share/reference/IMPACT410_genes_for_copynumber.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
# 			select(hgnc_symbol, chr, start_position, end_position) %>%
# 			rename(Hugo_Symbol = hgnc_symbol,
# 				   Chromosome = chr,
# 				   Start = start_position,
# 				   End = end_position) %>%
# 			mutate(Chromosome = ifelse(Chromosome %in% "X", 23, Chromosome)) %>%
# 			arrange(as.numeric(Chromosome), as.numeric(Start))
# 				   
# 	annot_by_gene <- annot %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Hugo_Symbol = Hugo_Symbol)
# 	res_by_segment <- res %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Calls = Calls)
# 	fo <- findOverlaps(res_by_segment, annot_by_gene)
# 
# 	df <- data.frame(Hugo_Symbol=mcols(annot_by_gene)[subjectHits(fo),], Calls=mcols(res_by_segment)[queryHits(fo),])
# 	Hugo_Symbol = which(duplicated(df$Hugo_Symbol))
# 	for (j in 1:length(Hugo_Symbol)) {
# 		index = which(as.character(df$Hugo_Symbol)==as.character(df$Hugo_Symbol[Hugo_Symbol[j]]))
# 		df[index,2] = round(mean(df[index,2], na.rm=TRUE))
# 	}
# 	df = df %>% filter(!duplicated(Hugo_Symbol))
# 	res = rep(0, nrow(annot))
# 	names(res) = annot[,"Hugo_Symbol"]
# 	res[as.character(df$Hugo_Symbol)] = df$Calls
# 	res[res<0] = -1
# 	res[res>0] = 1
# 	return(invisible(res))
# }
# i_cn = do.call(cbind, i_cn)
# colnames(i_cn) = tracker$GRAIL_ID
# 
# 
# g_cn = foreach (i=1:nrow(tracker)) %dopar% {
# 	cat(i, "\n")
# 	load(paste0("../res/rebuttal/GRAIL/facets/cncf/", tracker$GRAIL_ID[i], "_", tracker$GRAIL_ID[i], "-N.Rdata"))
# 	
# 	Chromosome = out2$jointseg[,"chrom"]
# 	Start = out2$jointseg[,"maploc"]
# 	End = out2$jointseg[,"maploc"]
# 	Calls = out2$jointseg[,"cnlr"]
# 	res = data.frame(Chromosome, Start, End, Calls)
# 	annot = read.csv(file="~/share/reference/IMPACT410_genes_for_copynumber.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
# 			select(hgnc_symbol, chr, start_position, end_position) %>%
# 			rename(Hugo_Symbol = hgnc_symbol,
# 				   Chromosome = chr,
# 				   Start = start_position,
# 				   End = end_position) %>%
# 			mutate(Chromosome = ifelse(Chromosome %in% "X", 23, Chromosome)) %>%
# 			arrange(as.numeric(Chromosome), as.numeric(Start))
# 				   
# 	annot_by_gene <- annot %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Hugo_Symbol = Hugo_Symbol)
# 	res_by_segment <- res %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Calls = Calls)
# 	fo <- findOverlaps(res_by_segment, annot_by_gene)
# 
# 	df <- data.frame(Hugo_Symbol=mcols(annot_by_gene)[subjectHits(fo),], Calls=mcols(res_by_segment)[queryHits(fo),])
# 	Hugo_Symbol = which(duplicated(df$Hugo_Symbol))
# 	for (j in 1:length(Hugo_Symbol)) {
# 		index = which(as.character(df$Hugo_Symbol)==as.character(df$Hugo_Symbol[Hugo_Symbol[j]]))
# 		df[index,2] = mean(df[index,2], na.rm=TRUE)
# 	}
# 	df = df %>% filter(!duplicated(Hugo_Symbol))
# 	res = rep(0, nrow(annot))
# 	names(res) = annot[,"Hugo_Symbol"]
# 	res[as.character(df$Hugo_Symbol)] = df$Calls
# 	return(invisible(res))
# }
# g_cn = do.call(cbind, g_cn)
# colnames(g_cn) = tracker$GRAIL_ID
# 
# g_cn_jrf = g_cn_cmo = g_cn
# 
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
