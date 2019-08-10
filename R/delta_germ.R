'delta_germ' <- function(zz)
{
  return(invisible(mean(pmin((zz - 0)^2, (zz - 0.5)^2, (zz - 1.0)^2), na.rm = T)))
}
