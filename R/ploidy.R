'ploidy' <- function(x,
					  y,
					  w)
{
	psi = seq(from=1.5, to=3.5, length=100)
	sse = vector(mode="numeric", length=length(psi))
	for (i in 1:length(psi)) {
		z = absolute_(rho=y, psi=psi[i], gamma=.85, x=x)
		sse[i] = sum(((w/sum(w)) * (z-round(z)))^2)
	}
	return(invisible(sse))
}
