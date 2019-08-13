'parse_msi_distribution' <- function(data_dir,
									 msi_dis_file)
{
	dd = readLines(sprintf("%s/%s",data_dir,msi_dis_file))
	dd_var_ndx = which(!grepl("N:|T:",dd))
	tidy_dd = data_frame(first=dd[dd_var_ndx],second=dd[dd_var_ndx+1],third=dd[dd_var_ndx+2])
	tidy_dd_f = tidy_dd %>%
				separate(first,c("chr","pos","flank_left","ms","flank_right"),sep=" ")
	tidy_dd_x = tidy_dd_f %>%
				mutate(mat=mapply(break_distribution,second,third,SIMPLIFY=FALSE,USE.NAMES=FALSE)) %>%
				select(-second,-third)
	return(invisible(tidy_dd_x))
}
