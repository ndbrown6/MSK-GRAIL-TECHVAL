'update_filter' <- function(filter,
							qualnobaq,
							pgtkxgdna,
							is_edge,
							min_p,
							min_qual = 60,
							sep = ";")
{
	high_qual = which(qualnobaq >= min_qual)
	low_qual = which(qualnobaq < min_qual)
	gdna = which(pgtkxgdna < min_p)
	edge = which(is_edge)
	pass = high_qual %>%
    	   setdiff(gdna) %>%
    	   setdiff(edge)
    low_qual_string = sprintf("QUAL_LT_%d", min_qual)
    gdna_strings = sprintf("PGTKXGDNA_LT_%0.2f", min_p)
    edge_string = "IS_EDGE"
  
	ids = data_frame(i = seq_along(filter))
	filters = bind_rows(
				data_frame(i = pass, filter = "PASS"),
				data_frame(i = low_qual, filter = low_qual_string),
				data_frame(i = gdna, filter = gdna_strings[gdna]),
				data_frame(i = edge, filter = edge_string))
	filters = filters %>%
			  group_by(i) %>%
			  summarize(filter = str_c(filter, collapse = sep)) %>%
			  ungroup() %>%
			  full_join(ids) %>%
			  arrange(i)
	return(invisible(filters$filter[filters$i]))
}
