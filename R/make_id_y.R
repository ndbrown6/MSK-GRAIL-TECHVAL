'make_id_y' <- function (patient_id,
						 hgvsp)
{
	return(invisible(paste0(patient_id, ":", hgvsp)))
}
