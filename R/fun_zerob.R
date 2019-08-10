'fun_zerob' <- function(x,
			y,
			n = 100,
			seed = 0)
{
	set.seed(seed)
	init = data.frame(y, x)
	y0 = NULL
	for (i in 1:n) {
		index = sample(x=(1:nrow(init)), size=nrow(init), replace=TRUE)
		data = init[index,,drop=FALSE]
		z = try(zeroinfl(y ~ x, dist = "poisson", data = data), silent=TRUE)
		if ("try-error" %in% is(z)) {
			y0[[i]] = rep(NA, times=100)
		} else {
			x0 = data.frame(x=seq(20, 90, l=100))
			y0[[i]] = predict(object=z, newdata=x0, type = "count")
		}
	}
	y0 = do.call(cbind, y0)
	return(invisible(y0))
}
