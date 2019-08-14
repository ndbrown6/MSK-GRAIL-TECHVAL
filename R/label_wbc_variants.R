'label_wbc_variants' <- function(wbc_stacked)
{
	n_samples = wbc_stacked %>%
				distinct(patient_id) %>%
				count()
	
	recurrence = wbc_stacked %>%
				 group_by(chrom, pos, ref, alt) %>%
				 summarize(n_recur = n(),
				 		   f_recur = n_recur / n_samples$n)
	
	wbc_stacked = wbc_stacked %>%
				  left_join(recurrence) %>%
				  mutate(filter = case_when(
				  		qualnobaq < 60 ~ "low_qual",
				  		dpnobaq < 500 ~ "low_depth",
				  		adnobaq / dpnobaq > 0.3 ~ "germline",
				  		f_recur > 0.05 ~ "recurrent",
				  		TRUE ~ "PASS"))
	return(invisible(wbc_stacked))
}
