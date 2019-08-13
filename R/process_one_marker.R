'process_one_marker' <- function(marker_mat)
{
  x_ndx = which(apply(marker_mat[2:3,],2,sum)>0)
  row_test = apply(marker_mat[2:3,],1,sum)
  x_sum = sum(row_test>=20)
  t_sum = sum(marker_mat[2:3,])
  valid = FALSE
  geno = NA
  xsq_pval = NA
  delta_entropy = NA
  js_metric = NA
  if (length(x_ndx)>1 & x_sum>1 & t_sum>=20) {
    valid = TRUE
    marker_vals = marker_mat[2:3,x_ndx]
    proportion_vals = sweep(marker_vals,1,apply(marker_vals,1,sum),FUN="/")
    ss_geno = sort(proportion_vals[1,],decreasing=TRUE)
    geno = TRUE
    if (ss_geno[1]>0.7) {
    	geno = FALSE
    }
    if (ss_geno[1]*0.66>ss_geno[2]) {
    	geno = FALSE
    }
    with_genotype = FALSE
    first_ratio = ss_geno[1]
    second_ratio = ss_geno[2]
    if (first_ratio>0.7) {
    	with_genotype=TRUE
    } else if (first_ratio>0.3 & second_ratio<first_ratio*0.66) {
    	with_genotype=TRUE
    } else if (first_ratio>0.3 & second_ratio>=first_ratio*0.66) {
    	with_genotype=TRUE
    }
    xsq_test = chisq.test(marker_vals)
    xsq_pval = xsq_test$p.value
    delta_entropy = diff(apply(proportion_vals,1,compute_entropy))
    js_metric = sqrt(JS(proportion_vals))
  }
  if (length(x_ndx)==1 & x_sum>1 & t_sum>=20) {
    valid = TRUE
    xsq_pval = 1
    delta_entropy = 0
    geno = FALSE
    js_metric = 0
  }
  return(invisible(list(valid = valid,
  						geno = geno,
  						xsq_pval = xsq_pval,
  						delta_entropy = delta_entropy,
  						js_metric = js_metric,
  						normal_depth = row_test[1],
  						tumor_depth = row_test[2])))
}
