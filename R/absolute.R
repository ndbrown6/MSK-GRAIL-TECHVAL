'absolute' <- function(rho,
		       psi,
		       gamma = 1,
		       x)
{
	rho = ifelse(is.na(rho), 1, rho)
	psi = ifelse(is.na(psi), 2, psi)
	return(invisible(((((2^(x/gamma))*(rho*psi+(1-rho)*2)) - ((1-rho)*2))/rho)))
}
