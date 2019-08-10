'in_common_bed' <- function (chromosome,
							 position)
{
	bed = read.csv(file=common_bed, header=FALSE, sep="\t", stringsAsFactors=FALSE)
	index = rep(FALSE, length(chromosome))
	for (i in 1:length(chromosome)) {
		indx = which(as.character(bed[,1])==as.character(chromosome[i]) & as.numeric(bed[,2])<=as.numeric(position[i]) & as.numeric(bed[,3])>=as.numeric(position[i]))
		if (length(indx)!=0) {
			index[i] = TRUE
		}
	}
	return(invisible(index))
}
