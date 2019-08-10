'label_bio_source' <- function(small_vars_plasma)
{
	n_samples <- small_vars_plasma %>%
    distinct(patient_id) %>%
    count()
    
    recurrence <- small_vars_plasma %>%
    group_by(loc) %>%
    summarize(n_recur = length(ref_orig),
              g_recur = sum(grepl("GDNA", filter) , na.rm = TRUE),
              q_germ = delta_germ(adgdna / (dpgdna + 0.5)),
              m_recur = mean(adgdna / dpgdna, na.rm=TRUE),
              f_recur = n_recur / n_samples$n,
              gf_recur = g_recur / n_samples$n)
  
	small_vars_plasma <- small_vars_plasma %>%
  					   mutate(gratio = (adnobaq + 2) * (dpgdna + 4) / ((adgdna + 2) * (dpnobaq + 4)),
  					   gzero = adgdna / sqrt(adgdna + 2)) %>%
  					   left_join(recurrence)
  
	variants <- small_vars_plasma %>%
  				mutate(category = case_when(
  			  			.$isedge == TRUE ~ "edge",
  			  			.$dpnobaq < 200 ~ "low_depth",
  			  			grepl("QUAL_LT", .$filter) ~ "low_qual",
  			  			(.$q_germ < 0.005 & .$adgdna / (.$dpgdna + 0.5) > germline_alpha) ~ "germline",
  			  			(.$adgdna / (.$dpgdna + 0.5) > germline_alpha) ~ "germlineish",
  			  			.$gf_recur > 0.05 ~ "artifact",
  			  			!is_nonsyn ~ "SSV",
  			  			grepl("GDNA",.$filter) & .$adgdna > 0 ~ "blood",
  			  			grepl("GDNA", .$filter) ~ "bloodish",
  			  			.$gzero>2 & .$gratio<4 ~ "bloodier",
  			  			grepl("PASS", .$filter) ~ "somatic",
  			  			TRUE ~ "other"))
 
	variants <- variants %>%
				mutate(bio_source = case_when(
  			  					category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
  			  					category %in% c("germline", "germlineish") ~ "germline",
  			  					category %in% c("blood", "bloodier") ~ "WBC_matched",
  			  					MSK == 1 ~ "biopsy_matched",
  			  					category == "somatic" ~ "VUSo",
  			  					TRUE ~ "other"))
  
	return(invisible(variants))
}
