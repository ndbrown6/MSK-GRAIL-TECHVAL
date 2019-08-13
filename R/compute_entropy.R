'compute_entropy' <- function(probs)
{
  return(invisible(-sum(probs*log2(probs),na.rm=TRUE)))
}
