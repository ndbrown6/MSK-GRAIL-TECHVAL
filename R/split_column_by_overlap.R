'split_column_by_overlap' <- function (query,
				       subject,
				       column = "ENTREZID",
				       ...)
{
    olaps = findOverlaps(query, subject, ...)
    f1 = factor(subjectHits(olaps), levels=seq_len(subjectLength(olaps)))
    return(invisible(splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)))
}
