'jensen_shannon' <- function(proportion_vals)
{
  return(invisible(compute_entropy(apply(proportion_vals,2,mean))-mean(apply(proportion_vals,1,compute_entropy))))
}
