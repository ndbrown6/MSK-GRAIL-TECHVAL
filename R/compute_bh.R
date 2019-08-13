'compute_bh' <- function(p_vals)
{
	base_rate = 0.05
	xtmp = sort(p_vals,decreasing=FALSE)
	xthreshold = 1:length(xtmp)/length(xtmp) * base_rate
	if (sum(xthreshold>xtmp)>0) {
		kk = max(which(xthreshold>xtmp))
		kthresh = xthreshold[kk]
	} else {
		kthresh = base_rate/length(xtmp)
	}
	return(invisible(kthresh))
}
