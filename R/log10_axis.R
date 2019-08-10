'log10_axis' <- function(side,
			 at,
			 ...)
{
    minor = NULL
    for (j in 2:length(at)) {
    	minor = c(minor, seq(from=at[j-1], to=at[j], length=10))
	}
	axis(side=side, at=minor, labels=NA, tcl=par("tcl")*0.65, ...)
}
