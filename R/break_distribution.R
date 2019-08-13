'break_distribution' <- function(x,
								 y)
{
	a = as.numeric(str_split(x," ",simplify=TRUE)[2:101])
	b = as.numeric(str_split(y," ",simplify=TRUE)[2:101])
	c = 1:100
	return(invisible(rbind(c,a,b))
}
