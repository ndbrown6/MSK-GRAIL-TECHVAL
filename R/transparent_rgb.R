'transparent_rgb' <- function (col = "black",
			       alpha = 85)
{
	tmp = c(col2rgb(col), alpha, 255)
	names(tmp) = c("red", "green", "blue", "alpha", "maxColorValue")
	out = do.call("rgb", as.list(tmp))
	return(invisible(out))
}
