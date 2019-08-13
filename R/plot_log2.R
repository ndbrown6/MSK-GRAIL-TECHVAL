'plot_log2' <- function(x,
						 y,
						 axis = TRUE,
						 ylim = c(-2.5, 2.5))
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
	abline(h=0, col="black", lty=1, lwd=1)
	if (axis) {
		axis(side=1, at=c(start, end[length(end)]), labels=rep("", length(start)+1), tcl=.5)
		axis(1, at = .5*(start+end), labels=c(1:22, "X"), tcl=-.5, lwd=0, lwd.ticks=1, tcl=-.25)
	}
	return(invisible(0))

}
