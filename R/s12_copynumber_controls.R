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

cols = colorRampPalette(c("#fee0d2", "#fc9272", "#de2d26"))(sum(unlist(lapply(res, function(x) {x[[3]]}))=="M"))
pdf(file="../res/figureS12/Histogram_Log2_Ratio_Healthy_M.pdf", width=5.5, height=7)
par(mar=c(6.1, 6.5, 4.1, 1.1))
plot(0, 0, type="n", xlab="", ylab="", frame.plot=FALSE, axes=FALSE, ylim=c(0,11), xlim=c(-1.5,1.5))
ii = 1
for (i in 1:length(fileNames)) {
	if (res[[i]]$z=="M") {
		points(res[[i]]$x, res[[i]]$y, type="l", col=cols[ii], lwd=2.5)
		ii = ii + 1
	}
}
points(c(0,0), c(-1,10), type="l", lty=3, col="goldenrod3", lwd=1.5)
points(-c(0.85,.85), c(-1,10), type="l", lty=3, col="goldenrod3", lwd=1.5)
axis(1, at=seq(-1.5,1.5, by=.5), labels=seq(-1.5,1.5, by=.5), cex.axis = 1.45, padj = 0.25, lwd=1.85, lwd.ticks=1.75)
axis(2, at=seq(0,10, by=2), labels=seq(0,10, by=2), cex.axis = 1.45, las = 1, lwd=1.85, lwd.ticks=1.75)
axis(2, at=seq(1,9, by=2), labels=rep("", 5), lwd=0, lwd.ticks=1, tcl=-.35)
mtext(side = 1, text = expression(Log[2]~"Ratio"), line = 4, cex = 1.55)
mtext(side = 2, text = "Density", line = 4, cex = 1.55)
dev.off()

cols = colorRampPalette(c("#efedf5", "#bcbddc", "#756bb1"))(sum(unlist(lapply(res, function(x) {x[[3]]}))=="F"))
pdf(file="../res/figureS12/Histogram_Log2_Ratio_Healthy_F.pdf", width=5.5, height=7)
par(mar=c(6.1, 6.5, 4.1, 1.1))
plot(0, 0, type="n", xlab="", ylab="", frame.plot=FALSE, axes=FALSE, ylim=c(0,11), xlim=c(-1.5,1.5))
ii = 1
for (i in 1:length(fileNames)) {
	if (res[[i]]$z=="F") {
		points(res[[i]]$x, res[[i]]$y, type="l", col=cols[ii], lwd=2.5)
		ii = ii + 1
	}
}
points(c(0,0), c(-1,10), type="l", lty=3, col="goldenrod3", lwd=1.5)
axis(1, at=seq(-1.5,1.5, by=.5), labels=seq(-1.5,1.5, by=.5), cex.axis = 1.45, padj = 0.25, lwd=1.85, lwd.ticks=1.75)
axis(2, at=seq(0,10, by=2), labels=seq(0,10, by=2), cex.axis = 1.45, las = 1, lwd=1.85, lwd.ticks=1.75)
axis(2, at=seq(1,9, by=2), labels=rep("", 5), lwd=0, lwd.ticks=1, tcl=-.35)
mtext(side = 1, text = expression(Log[2]~"Ratio"), line = 4, cex = 1.55)
mtext(side = 2, text = "Density", line = 4, cex = 1.55)
dev.off()

#==================================================
# genome-wide log2 ratios grail cfdna controls
#==================================================
fileNames = dir(path="../res/rebuttal/uncollapsed_bam/cnvkit/totalcopy", pattern="W", full.names=TRUE)
res = foreach (i=1:length(fileNames)) %dopar% {
	load(fileNames[i])
	sex = ifelse(mean(CN[CN[,"Chromosome"]==23,"Log2Ratio"])<(-.5), "M", "F")
	cn = smooth.spline(CN$Log2Ratio, spar=0.3)
	return(list(sex, cn))
}

pdf(file="../res/figureS12/Log2_Ratio_Healthy_M.pdf", width=12, height=5)
par(mar=c(6.1, 6.5, 4.1, 1.1))
plot(0, 0, type="n", xlab="", ylab="", frame.plot=FALSE, axes=FALSE, ylim=c(-1.25,1.25), xlim=c(0,19316))
cols = colorRampPalette(c("#fee0d2", "#fc9272", "#de2d26"))(sum(unlist(lapply(res, function(x) {x[[1]]}))=="M"))
for (i in 1:length(res)) {
	if (res[[i]][[1]]=="M") {
		points(res[[i]][[2]]$x, res[[i]][[2]]$y, type="l", col=cols[i], lwd=2)
	}
}
points(c(-550,19316), c(0,0), type="l", col="black")
points(c(-550,19316), c(-.85,-.85), type="l", col="black", lty=3)
axis(2, at=seq(-1,1, by=.25), labels=seq(-1,1, by=.25), las=1, cex.axis = 1, lwd=1.15, lwd.ticks=0.95, line=-.5)
mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 4, cex = 1.25)

load(fileNames[1])
ix = NULL
iy = NULL
for (i in 1:23) {
	ix = c(ix, min(which(CN[,"Chromosome"]==i)))
	if (i==23) {
		ix = c(ix, max(which(CN[,"Chromosome"]==i)))
	}
	iy = c(iy, mean(range(which(CN[,"Chromosome"]==i))))
}
axis(side=1, at=ix, labels=rep("", length(ix)), tcl=.5)
axis(side=1, at=iy, labels=c(1:22, "X"), tcl=-.5, lwd=0, lwd.ticks=1, tcl=-.25)
dev.off()

pdf(file="../res/figureS12/Log2_Ratio_Healthy_F.pdf", width=12, height=5)
par(mar=c(6.1, 6.5, 4.1, 1.1))
plot(0, 0, type="n", xlab="", ylab="", frame.plot=FALSE, axes=FALSE, ylim=c(-1.25,1.25), xlim=c(0,19316))
cols = colorRampPalette(c("#efedf5", "#bcbddc", "#756bb1"))(sum(unlist(lapply(res, function(x) {x[[1]]}))=="M"))
for (i in 1:length(res)) {
	if (res[[i]][[1]]=="F") {
		points(res[[i]][[2]]$x, res[[i]][[2]]$y, type="l", col=cols[i], lwd=2)
	}
}
points(c(-1000,19316), c(0,0), type="l", col="black")
axis(2, at=seq(-1,1, by=.25), labels=seq(-1,1, by=.25), las=1, cex.axis = 1, lwd=1.15, lwd.ticks=0.95, line=-.5)
mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 4, cex = 1.25)
dev.off()
