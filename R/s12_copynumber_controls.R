#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/figureS12")) {
	dir.create("../res/figureS12")
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
# histogram of log2 ratios grail cfdna controls
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

pdf(file="../res/figureS12/Histogram_Log2_Ratio_Healthy_M.pdf", width=6.5, height=7)
par(mar=c(6.1, 6.5, 4.1, 1.1))
plot(0, 0, type="n", xlab="", ylab="", frame.plot=FALSE, axes=FALSE, ylim=c(0,10.625), xlim=c(-1.5,1.5))
for (i in 1:length(fileNames)) {
	if (res[[i]]$z=="M") {
		points(res[[i]]$x, res[[i]]$y, type="l", col=transparent_rgb("salmon", 205), lwd=2.5)
	}
}
abline(v=0, lty=3, col="goldenrod3", lwd=1.5)
abline(v=-.85, lty=2, col="goldenrod3", lwd=1.5)
axis(1, at=seq(-1.5,1.5, by=.5), labels=seq(-1.5,1.5, by=.5), cex.axis = 1.75, padj = 0.25, lwd=1.85, lwd.ticks=1.75)
axis(2, at=seq(0,10, by=2), labels=seq(0,10, by=2), cex.axis = 1.75, las = 1, lwd=1.85, lwd.ticks=1.75)
mtext(side = 1, text = expression(Log[2]~"Ratio"), line = 4, cex = 1.75)
mtext(side = 2, text = "Density", line = 4, cex = 1.75)
title(main=paste0("Male controls (n=", sum(unlist(lapply(res, function(x) {x$z} ))=="M"), ")"), cex.main=1.0)
dev.off()

pdf(file="../res/figureS12/Histogram_Log2_Ratio_Healthy_F.pdf", width=6.5, height=7)
par(mar=c(6.1, 6.5, 4.1, 1.1))
plot(0, 0, type="n", xlab="", ylab="", frame.plot=FALSE, axes=FALSE, ylim=c(0,10.625), xlim=c(-1.5,1.5))
for (i in 1:length(fileNames)) {
	if (res[[i]]$z=="F") {
		points(res[[i]]$x, res[[i]]$y, type="l", col=transparent_rgb("steelblue", 205), lwd=2.5)
	}
}
abline(v=0, lty=3, col="goldenrod3", lwd=1.5)
abline(v=-.85, lty=2, col="goldenrod3", lwd=1.5)
axis(1, at=seq(-1.5,1.5, by=.5), labels=seq(-1.5,1.5, by=.5), cex.axis = 1.75, padj = 0.25, lwd=1.85, lwd.ticks=1.75)
axis(2, at=seq(0,10, by=2), labels=seq(0,10, by=2), cex.axis = 1.75, las = 1, lwd=1.85, lwd.ticks=1.75)
mtext(side = 1, text = expression(Log[2]~"Ratio"), line = 4, cex = 1.75)
mtext(side = 2, text = "Density", line = 4, cex = 1.75)
title(main=paste0("Female controls (n=", sum(unlist(lapply(res, function(x) {x$z} ))=="F"), ")"), cex.main=1.0)
dev.off()
