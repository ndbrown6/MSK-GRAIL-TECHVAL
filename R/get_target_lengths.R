'get_target_lengths' <- function()
{
	bed = read.csv(file=common_bed, header=FALSE, sep="\t", stringsAsFactors=FALSE)
	colnames(bed) = c("chrom", "start", "end")
	grbed = makeGRangesFromDataFrame(bed)
	gene_symbols = gene_ranges(Homo.sapiens, column="SYMBOL")
	res = unlist(split_column_by_overlap(gene_symbols, grbed, "SYMBOL"))
	gene_symbols = rep(NA, nrow(bed))
	gene_symbols[as.numeric(names(res))] = res
	bed = cbind(bed, gene_symbols)
	gene_symbols = as.character(unique(bed$gene_symbols))
	gene_symbols = gene_symbols[!is.na(gene_symbols)]
	target_lengths = unlist(foreach (i=1:length(gene_symbols)) %dopar% {
		index = which(as.character(bed[,"gene_symbols"])==gene_symbols[i])
		return(invisible(sum(bed[index,"end"] - bed[index,"start"])))
	})
	names(target_lengths) = gene_symbols
	return(invisible(target_lengths))
}
