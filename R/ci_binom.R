'ci_binom' <- function(n_succuss,
					   n_trial)
{
	test_res = binom.test(n_succuss, n_trial)
	ci = paste(round(100 * test_res$conf.int, 0), collapse = ",")
	return(invisible(ci))
}
