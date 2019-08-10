'order_samples' <- function(x)
{
	res = rep(1, length(x))
	index = grep("VL", x, fixed=TRUE)
	res[index] = 2
	index = grep("VP", x, fixed=TRUE)
	res[index] = 3
	return(invisible(res))
}
