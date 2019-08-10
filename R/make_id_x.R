'make_id_x' <- function (patient_id,
			 chrom,
			 pos,
			 ref,
			 alt)
{
	return(invisible(paste0(patient_id, ":", chrom, ":", pos, ":", ref, ">", alt)))
}
