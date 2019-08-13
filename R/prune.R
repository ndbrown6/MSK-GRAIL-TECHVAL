'prune' <- function(x,
					 n = 10)
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
	return(invisible(x))
}
