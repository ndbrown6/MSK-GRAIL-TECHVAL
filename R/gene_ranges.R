'gene_ranges' <-  function(db,
						   column="ENTREZID")
{
    g = genes(db, columns=column)
    col = mcols(g)[[column]]
    genes = granges(g)[rep(seq_along(g), elementNROWS(col))]
    mcols(genes)[[column]] = as.character(unlist(col))
    return(invisible(genes))
}
