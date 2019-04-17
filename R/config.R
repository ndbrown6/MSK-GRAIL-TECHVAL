#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("GenomicAlignments"))
suppressPackageStartupMessages(library("GenomicFeatures"))
suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene"))
suppressPackageStartupMessages(library("Homo.sapiens"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("annotate"))
suppressPackageStartupMessages(library("gdata"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("doMC"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("pander"))
suppressPackageStartupMessages(library("broom"))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("ggridges"))
suppressPackageStartupMessages(library("plotrix"))
suppressPackageStartupMessages(library("deconstructSigs"))
suppressPackageStartupMessages(library("gplots"))
suppressPackageStartupMessages(library("fuzzyjoin"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("googlesheets"))
suppressPackageStartupMessages(library("Publish"))
suppressPackageStartupMessages(library("Palimpsest"))
suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg19"))
suppressPackageStartupMessages(library("vioplot"))
suppressPackageStartupMessages(library("corrplot"))
suppressPackageStartupMessages(library("awtools"))
suppressPackageStartupMessages(library("beeswarm"))
suppressPackageStartupMessages(library("clinfun"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("mblm"))
suppressPackageStartupMessages(library("pscl"))
suppressPackageStartupMessages(library("boot"))
suppressPackageStartupMessages(library("snow"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("magrittr"))

#suppressPackageStartupMessages(library("eulerr"))
#suppressPackageStartupMessages(library("maftools"))
#suppressPackageStartupMessages(library("ggrepel"))

registerDoMC(32)

FLAG <- TRUE

snv_file <- list(
	scored = str_c(
		"../",
		"modified_v11/",
		"Variants_Calls/",
		"Joined_cfDNA_IMPACT_variants/",
		"scored_merged_snvs_20171115.tsv"
	)
)

indel_file <- list(
	scored = str_c(
		"../",
		"modified_v11/",
		"Variants_Calls/",
		"Joined_cfDNA_IMPACT_variants/",
		"scored_merged_indels_20171115.tsv"
	)
)

hypermutators <- list(
	patient_id = c(
		"MSK-VB-0023",
		"MSK-VB-0044",
		"MSK-VB-0046",
		"MSK-VB-0050",
		"MSK-VB-0057",
		"MSK-VL-0035",
		"MSK-VL-0054",
		"MSK-VP-0031",
		"MSK-VP-0054"
	)
)

msi_hypermutators <- list(
	patient_id = c(
		"MSK-VP-0041"
	)
)

common_bed <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"Bed_files/",
	"pan_v2_wo_decoy_wo_iSNP_wo_CNV_IMPACT_common_regions.merged.bed"
)

common_bed_annotated <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"Bed_files/",
	"pan_v2_wo_decoy_wo_iSNP_wo_CNV_IMPACT_common_regions.merged_annotated.bed"
)

common_bed_annotated_w_introns <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"Bed_files/",
	"pan_v2_wo_decoy_wo_iSNP_wo_CNV_IMPACT_common_regions.merged_annotated_w_introns.bed"
)

patient_tracker <- str_c(
	"../",
	"modified_v11/",
	"Tracker_files/",
	"s3_sample_tracker_TechVal_Merlin.csv"
)

impact_tracker <- str_c(
	"../",
	"modified_v11/",
    "Tracker_files/",
    "sample_tracker_IMPACT_BAM_20170227.csv"
)

wbc_variants <- list(
	scored = str_c(
		"../",
		"modified_v11/",
		"Variants_Calls/",
		"Stacked_Scored_WBC/",
		"wbc_scored_annotated.tsv"
	),
	annotations = str_c(
		"../",
		"res/",
		"etc/",
		"wbc_annotations.maf"
	)
)

clinical_file <- str_c(
	"../",
	"modified_v11/",
	"Variants_Calls/",
	"Joined_cfDNA_IMPACT_variants/",
	"clinical_no_dates_1_201703241549.tsv"
)

clinical_file_updated <- str_c(
	"../",
	"modified_v11/",
	"Variants_Calls/",
	"Joined_cfDNA_IMPACT_variants/",
	"clinical_no_dates_1_201703241549_updated.tsv"
)

somatic_vars_breast <- str_c(
	"../",
	"res/",
	"etc/",
	"normal_silent_common_intron_target.tables_merged_breast.tsv"
)

somatic_vars_lung <- str_c(
	"../",
	"res/",
	"etc/",
	"normal_silent_common_intron_target.tables_merged_lung.tsv"
)

somatic_vars_prostate <- str_c(
	"../",
	"res/",
	"etc/",
	"normal_silent_common_intron_target.tables_merged_prostate.tsv"
)

cosmic_file <- str_c(
	"../",
	"res/",
	"etc/",
	"cosmic_db_v84.RData"
)

gnomad_file <- str_c(
	"../",
	"res/",
	"etc/",
	"gnomad_db.r2.0.1.RData"
)

hotspot_file <- str_c(
	"../",
	"res/",
	"etc/",
	"cancer_hotspots_v2.RData"
)

all_vars_and_clinical <- str_c(
	"../",
	"res/",
	"etc/",
	"all_vars_and_clinical.RData"
)

chip_genes <- c(
	"DNMT3A",
	"TET2",
	"ASXL1",
	"PPM1D",
	"TP53",
	"JAK2",
	"RUNX1",
	"SF3B1",
	"SRSF2",
	"IDH1",
	"IDH2",
	"U2AF1",
	"CBL",
	"ATM",
	"CHEK2"
)

somatic_snvs_grail <- list(
	scored = str_c(
    "../",
    "modified_v11/",
    "Variants_Calls/",
    "Joined_cfDNA_IMPACT_variants/",
    "scored_merged_snvs_20171115.tsv"
  )
)

somatic_indels_grail <- list(
  scored = str_c(
    "../",
    "modified_v11/",
    "Variants_Calls/",
    "Joined_cfDNA_IMPACT_variants/",
    "scored_merged_indels_20171115.tsv"
  )
)

msk_anno_joined <- str_c("../",
						 "modified_v11/",
						 "Variants_Calls/",
						 "Joined_cfDNA_IMPACT_variants/",
						 "scripts2/",
						 "20180412_MSK_TechVal_Grail_small_variants_bio_source_label.filter_vus_biopsy_only.all_cases_ccf.tsv")
						 
						 
url_sample.tracker <- patient_tracker
url_msk.snv <- snv_file$scored
url_msk.indel <- indel_file$scored
url_techval.repeats <- "../annotated.tsv"
url_techval.repeats <- "../modified_v11/Variants_Calls/Joined_cfDNA_IMPACT_variants/annotated.tsv"
url_target.bed <- common_bed
url_ddpcr <- "https://docs.google.com/spreadsheets/d/1s9PKGYZ1v3vgUiC25N6D5jXX1ai6c0vbisZsKIgwGlg/edit?usp=sharing"
url_original <- url_msk.snv
url_retest <- "../modified_v11/Variants_Calls/Stacked_annotated_retest/TechVal_retest_annotated_stack.tsv"
url_cell.line <- "../res/etc/2018-01-10/scored/annotated.tsv"
url_HD753.truth <- "../res/etc/2018-01-10/hd753_manifest.tsv"

germline_alpha <- 0.15

bam_read_count_cutoff <- 2

ad_cutoff <- 5

ref_genome <- BSgenome.Hsapiens.UCSC.hg19

psa_file <- "../res/etc/PSA_MSK_VP_0041.txt"

recist_file <- "../res/etc/RECIST_MSK_VP_0041.txt"

'in_common_bed' <- function (chromosome, position)
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


'make_id_x' <- function (patient_id, chrom, pos, ref, alt)
{
	return(invisible(paste0(patient_id, ":", chrom, ":", pos, ":", ref, ">", alt)))
}


'make_id_y' <- function (patient_id, hgvsp)
{
	return(invisible(paste0(patient_id, ":", hgvsp)))
}


'transparentRgb' <- function (col = "black", alpha = 85)
{
	tmp = c(col2rgb(col), alpha, 255)
	names(tmp) = c("red", "green", "blue", "alpha", "maxColorValue")
    out = do.call("rgb", as.list(tmp))
    return(invisible(out))
}

'geneRanges' <-  function(db, column="ENTREZID")
{
    g = genes(db, columns=column)
    col = mcols(g)[[column]]
    genes = granges(g)[rep(seq_along(g), elementNROWS(col))]
    mcols(genes)[[column]] = as.character(unlist(col))
    return(invisible(genes))
}


'splitColumnByOverlap' <- function(query, subject, column="ENTREZID", ...)
{
    olaps = findOverlaps(query, subject, ...)
    f1 = factor(subjectHits(olaps), levels=seq_len(subjectLength(olaps)))
    return(invisible(splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)))
}


'get_target_lengths' <- function()
{
	bed = read.csv(file=common_bed, header=FALSE, sep="\t", stringsAsFactors=FALSE)
	colnames(bed) = c("chrom", "start", "end")
	grbed = makeGRangesFromDataFrame(bed)
	gene_symbols = geneRanges(Homo.sapiens, column="SYMBOL")
	res = unlist(splitColumnByOverlap(gene_symbols, grbed, "SYMBOL"))
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

'label_bio_source' <- function(small_vars_plasma) {
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
              gf_recur = g_recur / n_samples$n
    )
  
  small_vars_plasma <- small_vars_plasma %>%
    mutate(gratio = (adnobaq + 2) * (dpgdna + 4) / ((adgdna + 2) * (dpnobaq + 4)),
           gzero = adgdna / sqrt(adgdna + 2)) %>%
    left_join(recurrence)
  
  variants <- small_vars_plasma %>%
    mutate(
      category = case_when(
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
        TRUE ~ "other"
      )
    )
  
  variants <- variants %>%
    mutate(bio_source = case_when(
      category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
      category %in% c("germline", "germlineish") ~ "germline",
      category %in% c("blood", "bloodier") ~ "WBC_matched",
      MSK == 1 ~ "biopsy_matched",
      category == "somatic" ~ "VUSo",
      TRUE ~ "other")
    )
  
  return(variants)
}

'label_bio_source_hypermutator' <- function(small_vars_plasma, small_vars_hypermutator) {
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
  
  small_vars_hypermutator <- small_vars_hypermutator %>%
  							 mutate(gratio = (adnobaq + 2) * (dpgdna + 4) / ((adgdna + 2) * (dpnobaq + 4)), gzero = adgdna / sqrt(adgdna + 2)) %>%
  							 left_join(recurrence)
  
  variants <- small_vars_hypermutator %>%
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
  
  return(variants)
}

'label_wbc_variants' <- function(wbc_stacked) {
  n_samples <- wbc_stacked %>%
  			   distinct(patient_id) %>%
  			   count()
  
  recurrence <- wbc_stacked %>%
  				group_by(chrom, pos, ref, alt) %>%
  				summarize(n_recur = n(), f_recur = n_recur / n_samples$n)
  
  wbc_stacked <- wbc_stacked %>%
  				 left_join(recurrence) %>%
  				 mutate(filter = case_when(
  				 				qualnobaq < 60 ~ "low_qual",
  				 				dpnobaq < 500 ~ "low_depth",
  				 				adnobaq / dpnobaq > germline_alpha ~ "germline",
  				 				f_recur > 0.05 ~ "recurrent",
  				 				TRUE ~ "PASS"))
  
  return(wbc_stacked)
}

'delta_germ' <- function(zz) {
  mean(pmin((zz - 0)^2, (zz - 0.5)^2, (zz - 1.0)^2), na.rm = T)
}

'update_filter' <- function (filter, qualnobaq, pgtkxgdna, is_edge, min_p,
                          min_qual = 60, sep = ";")
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
	filters$filter[filters$i]
}

'plot96_mutation_spectrum' <- function (vcf, sample.col = "sample", mutcat3.col = "mutcat3",
										ymax = NULL, averageProp = FALSE, plot.file = NULL)
{
    bases <- c("A", "C", "G", "T")
    ctxt16 <- paste(rep(bases, each = 4), rep(bases, 4), sep = ".")
    mt <- c("CA", "CG", "CT", "TA", "TC", "TG")
    types96 <- paste(rep(mt, each = 16), rep(ctxt16, 6), sep = "_")
    types96 <- sapply(types96, function(z) {
        sub("\\.", substr(z, 1, 1), z)
    })
    context <- substr(types96, 4, 6)
    nsamp <- length(unique(vcf[, sample.col]))
    if (averageProp & nsamp > 1) {
        tmp <- makeMutypeMatFromVcf(vcf, sample.col = "CHCID", 
            mutcat.col = "mutcat3", mutypes = types96)
        freq <- apply(tmp, 1, mean)
    }
    else {
        freq <- sapply(types96, function(z) {
            mean(vcf[, mutcat3.col] == z, na.rm = T)
        })
    }
    if (!is.null(plot.file)) {
        pdf(plot.file, width = 24, height = 5)
    }
    col96 <- c(rep("skyblue3", 16), rep("black", 16), rep("red", 
        16), rep("grey", 16), rep("green", 16), rep("pink", 16))
    labs <- c(rep("C>A", 16), rep("C>G", 16), rep("C>T", 16), 
        rep("T>A", 16), rep("T>C", 16), rep("T>G", 16))
    if (is.null(ymax)) {
        ymax <- 100*ceiling(max(freq) * 100)/100
        ymax <- ifelse(ymax>10, 30, 10)
    }
    bp <- barplot(freq*100, col = col96, border = col96, las = 2, 
        width = 1, space = .35, yaxt = "n", xaxt = "n", ylim = c(0, 
            ymax * 1.2))
    title(ylab = "Fraction of mutations (%)", mgp = c(1, 1, 0), 
        cex.lab = 1.6)
    axis(1, at = bp, labels = context, pos = 0, las = 2, cex.axis = 1.5, 
        tick = F, cex.axis = 1, lwd=-1)
    if (ymax==40) {
	    axis(2, at = c(0,10,20,30,40), labels=c(0,10,20,30,40), pos = 0, las = 1, cex.axis = 1.5)
	} else if (ymax==30) {
	    axis(2, at = c(0,5,10,15,20,25,30), labels=c(0,5,10,15,20,25,30), pos = 0, las = 1, cex.axis = 1.5)
	} else if (ymax==20) {
		axis(2, at = c(0,5,10,15,20), labels=c(0,5,10,15,20), pos = 0, las = 1, cex.axis = 1.5)
	} else if (ymax==10) {
		axis(2, at = c(0,2,4,6,8,10), labels=c(0,2,4,6,8,10), pos = 0, las = 1, cex.axis = 1.5)
	}
    for (i in seq(1, 81, by = 16)) {
        rect(bp[i], par()$usr[4], bp[i + 15], par()$usr[4] - 
            0.05 * diff(par()$usr[3:4]), col = col96[i], border = col96[i])
        text((bp[i] + bp[i + 15])/2, par()$usr[4] + 0.09 * diff(par()$usr[3:4]), 
            labels = labs[i], xpd = TRUE, cex = 2)
    }
    if (!is.null(plot.file)) {
        dev.off()
    }
}

'corrplot2' <- function (corr, corr2, method = c("circle", "square", "ellipse", "number",
    "shade", "color", "pie"), type = c("full", "lower", "upper"), 
    add = FALSE, col = NULL, bg = "white", title = "", is.corr = TRUE, 
    diag = TRUE, outline = FALSE, mar = c(0, 0, 0, 0), addgrid.col = NULL, 
    addCoef.col = NULL, addCoefasPercent = FALSE, order = c("original", 
    "AOE", "FPC", "hclust", "alphabet"), hclust.method = c("complete", 
    "ward", "ward.D", "ward.D2", "single", "average", "mcquitty", 
    "median", "centroid"), addrect = NULL, rect.col = "black", 
    rect.lwd = 2, tl.pos = NULL, tl.cex = 1, tl.col = "red", 
    tl.offset = 0.4, tl.srt = 90, cl.pos = NULL, cl.lim = NULL, cl.lim2 = NULL, 
    cl.length = NULL, cl.cex = 0.8, cl.ratio = 0.15, cl.align.text = "c", 
    cl.offset = 0.5, number.cex = 1, number.font = 2, number.digits = NULL, 
    addshade = c("negative", "positive", "all"), shade.lwd = 1, 
    shade.col = "white", p.mat = NULL, sig.level = 0.05, insig = c("pch", 
    "p-value", "blank", "n", "label_sig"), pch = 4, pch.col = "black", 
    pch.cex = 3, plotCI = c("n", "square", "circle", "rect"), 
    lowCI.mat = NULL, uppCI.mat = NULL, na.label = "?", na.label.col = "black", 
    win.asp = 1, ...) 
{
    method <- match.arg(method)
    type <- match.arg(type)
    order <- match.arg(order)
    hclust.method <- match.arg(hclust.method)
    addshade <- match.arg(addshade)
    insig <- match.arg(insig)
    plotCI <- match.arg(plotCI)
    if (win.asp != 1 && !(method %in% c("circle", "square"))) {
        stop("Parameter 'win.asp' is supported only for circle and square methods.")
    }
    asp_rescale_factor <- min(1, win.asp)/max(1, win.asp)
    stopifnot(asp_rescale_factor >= 0 && asp_rescale_factor <= 
        1)
    if (!is.matrix(corr) && !is.data.frame(corr)) {
        stop("Need a matrix or data frame!")
    }
    if (is.null(addgrid.col)) {
        addgrid.col <- switch(method, color = NA, shade = NA, 
            "grey")
    }
    if (any(corr < cl.lim[1]) || any(corr > cl.lim[2])) {
        stop("color limits should cover matrix")
    }
    if (is.null(cl.lim)) {
        if (is.corr) {
            cl.lim <- c(-1, 1)
        }
        else {
            corr_tmp <- corr
            diag(corr_tmp) <- ifelse(diag, diag(corr_tmp), NA)
            cl.lim <- c(min(corr_tmp, na.rm = TRUE), max(corr_tmp, 
                na.rm = TRUE))
        }
    }
    if (is.null(cl.lim2)) {
       corr2_tmp <- corr2
       diag(corr2_tmp) <- ifelse(diag, diag(corr2_tmp), NA)
       cl.lim2 <- c(min(corr2_tmp, na.rm = TRUE), max(corr2_tmp, na.rm = TRUE))
    }
    intercept = intercept2 <- 0
    zoom = zoom2 <- 1
    if (!is.corr) {
        c_max <- max(corr, na.rm = TRUE); c_max2 <- max(corr2, na.rm = TRUE)
        c_min <- min(corr, na.rm = TRUE); c_min2 <- min(corr2, na.rm = TRUE)
        
        if (c_max <= 0) {
            intercept <- -cl.lim[2]
            zoom <- 1/(diff(cl.lim))
        }
        else if (c_min >= 0) {
            intercept <- -cl.lim[1]
            zoom <- 1/(diff(cl.lim))
        }
        else {
            stopifnot(c_max * c_min < 0)
            stopifnot(c_min < 0 && c_max > 0)
            intercept <- 0
            zoom <- 1/max(abs(cl.lim))
        }
        
        if (c_max2 <= 0) {
            intercept2 <- -cl.lim2[2]
            zoom2 <- 1/(diff(cl.lim2))
        }
        else if (c_min2 >= 0) {
            intercept2 <- -cl.lim2[1]
            zoom2 <- 1/(diff(cl.lim2))
        }
        else {
            stopifnot(c_max2 * c_min2 < 0)
            stopifnot(c_min2 < 0 && c_max2 > 0)
            intercept2 <- 0
            zoom2 <- 1/max(abs(cl.lim2))
        }
        
        if (zoom == Inf) {
            stopifnot(cl.lim[1] == 0 && cl.lim[2] == 0)
            zoom <- 0
        }
        corr <- (intercept + corr) * zoom
        
        if (zoom2 == Inf) {
            stopifnot(cl.lim2[1] == 0 && cl.lim2[2] == 0)
            zoom2 <- 0
        }
        corr2 <- (intercept2 + corr2) * zoom2

    }
    cl.lim2 <- (intercept2 + cl.lim2) * zoom2
    int <- intercept * zoom
    if (is.corr) {
        if (min(corr, na.rm = TRUE) < -1 - .Machine$double.eps^0.75 || 
            max(corr, na.rm = TRUE) > 1 + .Machine$double.eps^0.75) {
            stop("The matrix is not in [-1, 1]!")
        }
    }
    if (is.null(col)) {
        col <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
            "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
            "#4393C3", "#2166AC", "#053061"))(200)
    }
    n <- nrow(corr)
    m <- ncol(corr)
    min.nm <- min(n, m)
    ord <- seq_len(min.nm)
    if (order != "original") {
        ord <- corrMatOrder(corr, order = order, hclust.method = hclust.method)
        corr <- corr[ord, ord]
        corr2 <- corr2[ord, ord]
    }
    if (is.null(rownames(corr))) {
        rownames(corr) <- seq_len(n)
    }
    if (is.null(colnames(corr))) {
        colnames(corr) <- seq_len(m)
    }
    apply_mat_filter <- function(mat) {
        x <- matrix(1:n * m, nrow = n, ncol = m)
        switch(type, upper = mat[row(x) > col(x)] <- Inf, lower = mat[row(x) < 
            col(x)] <- Inf)
        if (!diag) {
            diag(mat) <- Inf
        }
        return(mat)
    }
    getPos.Dat <- function(mat) {
        tmp <- apply_mat_filter(mat)
        Dat <- tmp[is.finite(tmp)]
        ind <- which(is.finite(tmp), arr.ind = TRUE)
        Pos <- ind
        Pos[, 1] <- ind[, 2]
        Pos[, 2] <- -ind[, 1] + 1 + n
        return(list(Pos, Dat))
    }
    getPos.NAs <- function(mat) {
        tmp <- apply_mat_filter(mat)
        ind <- which(is.na(tmp), arr.ind = TRUE)
        Pos <- ind
        Pos[, 1] <- ind[, 2]
        Pos[, 2] <- -ind[, 1] + 1 + n
        return(Pos)
    }
    Pos <- getPos.Dat(corr)[[1]]
    if (any(is.na(corr)) && is.character(na.label)) {
        PosNA <- getPos.NAs(corr)
    }
    else {
        PosNA <- NULL
    }
    AllCoords <- rbind(Pos, PosNA)
    n2 <- max(AllCoords[, 2])
    n1 <- min(AllCoords[, 2])
    nn <- n2 - n1
    m2 <- max(AllCoords[, 1])
    m1 <- min(AllCoords[, 1])
    mm <- max(1, m2 - m1)
    expand_expression <- function(s) {
        ifelse(grepl("^[:=$]", s), parse(text = substring(s, 
            2)), s)
    }
    newrownames <- sapply(rownames(corr)[(n + 1 - n2):(n + 1 - 
        n1)], expand_expression)
    newcolnames <- sapply(colnames(corr)[m1:m2], expand_expression)
    DAT <- getPos.Dat(corr)[[2]]; DAT2 <- getPos.Dat(corr2)[[2]]
    len.DAT <- length(DAT)
    rm(expand_expression)
    assign.color <- function(dat = DAT, color = col) {
        newcorr <- (dat + 1)/2
        newcorr[newcorr <= 0] <- 0
        newcorr[newcorr >= 1] <- 1 - 1e-16
        color[floor(newcorr * length(color)) + 1]
    }
    col.fill <- assign.color(dat = DAT2, color = col)
    isFALSE <- function(x) identical(x, FALSE)
    isTRUE <- function(x) identical(x, TRUE)
    if (isFALSE(tl.pos)) {
        tl.pos <- "n"
    }
    if (is.null(tl.pos) || isTRUE(tl.pos)) {
        tl.pos <- switch(type, full = "lt", lower = "ld", upper = "td")
    }
    if (isFALSE(cl.pos)) {
        cl.pos <- "n"
    }
    if (is.null(cl.pos) || isTRUE(cl.pos)) {
        cl.pos <- switch(type, full = "r", lower = "b", upper = "r")
    }
    if (isFALSE(outline)) {
        col.border <- col.fill
    }
    else if (isTRUE(outline)) {
        col.border <- "black"
    }
    else if (is.character(outline)) {
        col.border <- outline
    }
    else {
        stop("Unsupported value type for parameter outline")
    }
    oldpar <- par(mar = mar, bg = "white")
    on.exit(par(oldpar), add = TRUE)
    if (!add) {
        plot.new()
        xlabwidth <- max(strwidth(newrownames, cex = tl.cex))
        ylabwidth <- max(strwidth(newcolnames, cex = tl.cex))
        laboffset <- strwidth("W", cex = tl.cex) * tl.offset
        for (i in 1:50) {
            xlim <- c(m1 - 0.5 - laboffset - xlabwidth * (grepl("l", 
                tl.pos) | grepl("d", tl.pos)), m2 + 0.5 + mm * 
                cl.ratio * (cl.pos == "r") + xlabwidth * abs(cos(tl.srt * 
                pi/180)) * grepl("d", tl.pos)) + c(-0.35, 0.15) + 
                c(-1, 0) * grepl("l", tl.pos)
            ylim <- c(n1 - 0.5 - nn * cl.ratio * (cl.pos == "b") - 
                laboffset, n2 + 0.5 + laboffset + ylabwidth * 
                abs(sin(tl.srt * pi/180)) * grepl("t", tl.pos)) + 
                c(-0.15, 0) + c(0, -1) * (type == "upper" && 
                tl.pos != "n") + c(0, 1) * grepl("d", tl.pos)
            plot.window(xlim, ylim, asp = 1, xaxs = "i", yaxs = "i")
            x.tmp <- max(strwidth(newrownames, cex = tl.cex))
            y.tmp <- max(strwidth(newcolnames, cex = tl.cex))
            laboffset.tmp <- strwidth("W", cex = tl.cex) * tl.offset
            if (max(x.tmp - xlabwidth, y.tmp - ylabwidth, laboffset.tmp - 
                laboffset) < 0.001) {
                break
            }
            xlabwidth <- x.tmp
            ylabwidth <- y.tmp
            laboffset <- laboffset.tmp
            if (i == 50) {
                warning(c("Not been able to calculate text margin, ", 
                  "please try again with a clean new empty window using ", 
                  "{plot.new(); dev.off()} or reduce tl.cex"))
            }
        }
        if (.Platform$OS.type == "windows") {
            grDevices::windows.options(width = 7, height = 7 * 
                diff(ylim)/diff(xlim))
        }
        plot.window(xlim = xlim, ylim = ylim, asp = win.asp, 
            xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    }
    laboffset <- strwidth("W", cex = tl.cex) * tl.offset
    symbols(Pos, add = TRUE, inches = FALSE, rectangles = matrix(1, 
        len.DAT, 2), bg = bg, fg = bg)
    if (method == "circle" && plotCI == "n") {
        symbols(Pos, add = TRUE, inches = FALSE, circles = asp_rescale_factor * 
            0.9 * abs(DAT)^0.5/2, fg = "black", bg = col.fill)
    }
    if (method == "ellipse" && plotCI == "n") {
        ell.dat <- function(rho, length = 99) {
            k <- seq(0, 2 * pi, length = length)
            x <- cos(k + acos(rho)/2)/2
            y <- cos(k - acos(rho)/2)/2
            cbind(rbind(x, y), c(NA, NA))
        }
        ELL.dat <- lapply(DAT, ell.dat)
        ELL.dat2 <- 0.85 * matrix(unlist(ELL.dat), ncol = 2, 
            byrow = TRUE)
        ELL.dat2 <- ELL.dat2 + Pos[rep(1:length(DAT), each = 100), 
            ]
        polygon(ELL.dat2, border = col.border, col = col.fill)
    }
    if (is.null(number.digits)) {
        number.digits <- switch(addCoefasPercent + 1, 2, 0)
    }
    stopifnot(number.digits%%1 == 0)
    stopifnot(number.digits >= 0)
    if (method == "number" && plotCI == "n") {
        text(Pos[, 1], Pos[, 2], font = number.font, col = col.fill, 
            labels = round((DAT - int) * ifelse(addCoefasPercent, 
                100, 1)/zoom, number.digits), cex = number.cex)
    }
    NA_LABEL_MAX_CHARS <- 2
    if (is.matrix(PosNA) && nrow(PosNA) > 0) {
        stopifnot(is.matrix(PosNA))
        if (na.label == "square") {
            symbols(PosNA, add = TRUE, inches = FALSE, squares = rep(1, 
                nrow(PosNA)), bg = na.label.col, fg = na.label.col)
        }
        else if (nchar(na.label) %in% 1:NA_LABEL_MAX_CHARS) {
            symbols(PosNA, add = TRUE, inches = FALSE, squares = rep(1, 
                nrow(PosNA)), fg = bg, bg = bg)
            text(PosNA[, 1], PosNA[, 2], font = number.font, 
                col = na.label.col, labels = na.label, cex = number.cex, 
                ...)
        }
        else {
            stop(paste("Maximum number of characters for NA label is:", 
                NA_LABEL_MAX_CHARS))
        }
    }
    if (method == "pie" && plotCI == "n") {
        symbols(Pos, add = TRUE, inches = FALSE, circles = rep(0.5, 
            len.DAT) * 0.85, fg = col.border)
        pie.dat <- function(theta, length = 100) {
            k <- seq(pi/2, pi/2 - theta, length = 0.5 * length * 
                abs(theta)/pi)
            x <- c(0, cos(k)/2, 0)
            y <- c(0, sin(k)/2, 0)
            cbind(rbind(x, y), c(NA, NA))
        }
        PIE.dat <- lapply(DAT * 2 * pi, pie.dat)
        len.pie <- unlist(lapply(PIE.dat, length))/2
        PIE.dat2 <- 0.85 * matrix(unlist(PIE.dat), ncol = 2, 
            byrow = TRUE)
        PIE.dat2 <- PIE.dat2 + Pos[rep(1:length(DAT), len.pie), 
            ]
        polygon(PIE.dat2, border = "black", col = col.fill)
    }
    if (method == "shade" && plotCI == "n") {
        symbols(Pos, add = TRUE, inches = FALSE, squares = rep(1, 
            len.DAT), bg = col.fill, fg = addgrid.col)
        shade.dat <- function(w) {
            x <- w[1]
            y <- w[2]
            rho <- w[3]
            x1 <- x - 0.5
            x2 <- x + 0.5
            y1 <- y - 0.5
            y2 <- y + 0.5
            dat <- NA
            if ((addshade == "positive" || addshade == "all") && 
                rho > 0) {
                dat <- cbind(c(x1, x1, x), c(y, y1, y1), c(x, 
                  x2, x2), c(y2, y2, y))
            }
            if ((addshade == "negative" || addshade == "all") && 
                rho < 0) {
                dat <- cbind(c(x1, x1, x), c(y, y2, y2), c(x, 
                  x2, x2), c(y1, y1, y))
            }
            return(t(dat))
        }
        pos_corr <- rbind(cbind(Pos, DAT))
        pos_corr2 <- split(pos_corr, 1:nrow(pos_corr))
        SHADE.dat <- matrix(na.omit(unlist(lapply(pos_corr2, 
            shade.dat))), byrow = TRUE, ncol = 4)
        segments(SHADE.dat[, 1], SHADE.dat[, 2], SHADE.dat[, 
            3], SHADE.dat[, 4], col = shade.col, lwd = shade.lwd)
    }
    if (method == "square" && plotCI == "n") {
        draw_method_square(Pos, DAT, asp_rescale_factor, col.border, 
            col.fill)
    }
    if (method == "color" && plotCI == "n") {
        draw_method_color(Pos, col.border, col.fill)
    }
    corrplot:::draw_grid(AllCoords, addgrid.col)
    if (plotCI != "n") {
        if (is.null(lowCI.mat) || is.null(uppCI.mat)) {
            stop("Need lowCI.mat and uppCI.mat!")
        }
        if (order != "original") {
            lowCI.mat <- lowCI.mat[ord, ord]
            uppCI.mat <- uppCI.mat[ord, ord]
        }
        pos.lowNew <- getPos.Dat(lowCI.mat)[[1]]
        lowNew <- getPos.Dat(lowCI.mat)[[2]]
        pos.uppNew <- getPos.Dat(uppCI.mat)[[1]]
        uppNew <- getPos.Dat(uppCI.mat)[[2]]
        if (!method %in% c("circle", "square")) {
            stop("Method shoud be circle or square if drawing confidence intervals.")
        }
        k1 <- (abs(uppNew) > abs(lowNew))
        bigabs <- uppNew
        bigabs[which(!k1)] <- lowNew[!k1]
        smallabs <- lowNew
        smallabs[which(!k1)] <- uppNew[!k1]
        sig <- sign(uppNew * lowNew)
        color_bigabs <- col[ceiling((bigabs + 1) * length(col)/2)]
        color_smallabs <- col[ceiling((smallabs + 1) * length(col)/2)]
        if (plotCI == "circle") {
            symbols(pos.uppNew[, 1], pos.uppNew[, 2], add = TRUE, 
                inches = FALSE, circles = 0.95 * abs(bigabs)^0.5/2, 
                bg = ifelse(sig > 0, col.fill, color_bigabs), 
                fg = ifelse(sig > 0, col.fill, color_bigabs))
            symbols(pos.lowNew[, 1], pos.lowNew[, 2], add = TRUE, 
                inches = FALSE, circles = 0.95 * abs(smallabs)^0.5/2, 
                bg = ifelse(sig > 0, bg, color_smallabs), fg = ifelse(sig > 
                  0, col.fill, color_smallabs))
        }
        if (plotCI == "square") {
            symbols(pos.uppNew[, 1], pos.uppNew[, 2], add = TRUE, 
                inches = FALSE, squares = abs(bigabs)^0.5, bg = ifelse(sig > 
                  0, col.fill, color_bigabs), fg = ifelse(sig > 
                  0, col.fill, color_bigabs))
            symbols(pos.lowNew[, 1], pos.lowNew[, 2], add = TRUE, 
                inches = FALSE, squares = abs(smallabs)^0.5, 
                bg = ifelse(sig > 0, bg, color_smallabs), fg = ifelse(sig > 
                  0, col.fill, color_smallabs))
        }
        if (plotCI == "rect") {
            rect.width <- 0.25
            rect(pos.uppNew[, 1] - rect.width, pos.uppNew[, 2] + 
                smallabs/2, pos.uppNew[, 1] + rect.width, pos.uppNew[, 
                2] + bigabs/2, col = col.fill, border = col.fill)
            segments(pos.lowNew[, 1] - rect.width, pos.lowNew[, 
                2] + DAT/2, pos.lowNew[, 1] + rect.width, pos.lowNew[, 
                2] + DAT/2, col = "black", lwd = 1)
            segments(pos.uppNew[, 1] - rect.width, pos.uppNew[, 
                2] + uppNew/2, pos.uppNew[, 1] + rect.width, 
                pos.uppNew[, 2] + uppNew/2, col = "black", lwd = 1)
            segments(pos.lowNew[, 1] - rect.width, pos.lowNew[, 
                2] + lowNew/2, pos.lowNew[, 1] + rect.width, 
                pos.lowNew[, 2] + lowNew/2, col = "black", lwd = 1)
            segments(pos.lowNew[, 1] - 0.5, pos.lowNew[, 2], 
                pos.lowNew[, 1] + 0.5, pos.lowNew[, 2], col = "grey70", 
                lty = 3)
        }
    }
    if (!is.null(p.mat) && insig != "n") {
        if (order != "original") {
            p.mat <- p.mat[ord, ord]
        }
        pos.pNew <- getPos.Dat(p.mat)[[1]]
        pNew <- getPos.Dat(p.mat)[[2]]
        if (insig == "label_sig") {
            if (!is.character(pch)) 
                pch <- "*"
            place_points <- function(sig.locs, point) {
                text(pos.pNew[, 1][sig.locs], pos.pNew[, 2][sig.locs], 
                  labels = point, col = pch.col, cex = pch.cex, 
                  lwd = 2)
            }
            if (length(sig.level) == 1) {
                place_points(sig.locs = which(pNew < sig.level), 
                  point = pch)
            }
            else {
                l <- length(sig.level)
                for (i in seq_along(sig.level)) {
                  iter <- l + 1 - i
                  pchTmp <- paste(rep(pch, i), collapse = "")
                  if (i == length(sig.level)) {
                    locs <- which(pNew < sig.level[iter])
                    if (length(locs)) {
                      place_points(sig.locs = locs, point = pchTmp)
                    }
                  }
                  else {
                    locs <- which(pNew < sig.level[iter] & pNew > 
                      sig.level[iter - 1])
                    if (length(locs)) {
                      place_points(sig.locs = locs, point = pchTmp)
                    }
                  }
                }
            }
        }
        else {
            ind.p <- which(pNew > sig.level)
            p_inSig <- length(ind.p) > 0
            if (insig == "pch" && p_inSig) {
                points(pos.pNew[, 1][ind.p], pos.pNew[, 2][ind.p], 
                  pch = pch, col = pch.col, cex = pch.cex, lwd = 2)
            }
            if (insig == "p-value" && p_inSig) {
                text(pos.pNew[, 1][ind.p], pos.pNew[, 2][ind.p], 
                  round(pNew[ind.p], 2), col = pch.col)
            }
            if (insig == "blank" && p_inSig) {
                symbols(pos.pNew[, 1][ind.p], pos.pNew[, 2][ind.p], 
                  inches = FALSE, squares = rep(1, length(pos.pNew[, 
                    1][ind.p])), fg = addgrid.col, bg = bg, add = TRUE)
            }
        }
    }
    if (cl.pos != "n") {
        colRange <- assign.color(dat = cl.lim2)
        ind1 <- which(col == colRange[1])
        ind2 <- which(col == colRange[2])
        colbar <- col[ind1:ind2]
        if (is.null(cl.length)) {
            cl.length <- ifelse(length(colbar) > 20, 11, length(colbar) + 
                1)
        }
        labels <- seq(cl.lim[1], cl.lim[2], length = cl.length)
        if (cl.pos == "r") {
            vertical <- TRUE
            xlim <- c(m2 + 0.5 + mm * 0.02, m2 + 0.5 + mm * cl.ratio)
            ylim <- c(n1 - 0.5, n2 + 0.5)
        }
        if (cl.pos == "b") {
            vertical <- FALSE
            xlim <- c(m1 - 0.5, m2 + 0.5)
            ylim <- c(n1 - 0.5 - nn * cl.ratio, n1 - 0.5 - nn * 
                0.02)
        }
        colorlegend(colbar = colbar, labels = round(labels, 2), 
            offset = cl.offset, ratio.colbar = 0.3, cex = cl.cex, 
            xlim = xlim, ylim = ylim, vertical = vertical, align = cl.align.text)
    }
    if (tl.pos != "n") {
        pos.xlabel <- cbind(m1:m2, n2 + 0.5 + laboffset)
        pos.ylabel <- cbind(m1 - 0.5, n2:n1)
        if (tl.pos == "td") {
            if (type != "upper") {
                stop("type should be \"upper\" if tl.pos is \"dt\".")
            }
            pos.ylabel <- cbind(m1:(m1 + nn) - 0.5, n2:n1)
        }
        if (tl.pos == "ld") {
            if (type != "lower") {
                stop("type should be \"lower\" if tl.pos is \"ld\".")
            }
            pos.xlabel <- cbind(m1:m2, n2:(n2 - mm) + 0.5 + laboffset)
        }
        if (tl.pos == "d") {
            pos.ylabel <- cbind(m1:(m1 + nn) - 0.5, n2:n1)
            pos.ylabel <- pos.ylabel[1:min(n, m), ]
            symbols(pos.ylabel[, 1] + 0.5, pos.ylabel[, 2], add = TRUE, 
                bg = bg, fg = addgrid.col, inches = FALSE, squares = rep(1, 
                  length(pos.ylabel[, 1])))
            text(pos.ylabel[, 1] + 0.5, pos.ylabel[, 2], newcolnames[1:min(n, 
                m)], col = tl.col, cex = tl.cex, ...)
        }
        else {
            text(pos.xlabel[, 1], pos.xlabel[, 2], newcolnames, 
                srt = tl.srt, adj = ifelse(tl.srt == 0, c(0.5, 
                  0), c(0, 0)), col = tl.col, cex = tl.cex, offset = tl.offset, 
                ...)
            text(pos.ylabel[, 1], pos.ylabel[, 2], newrownames, 
                col = tl.col, cex = tl.cex, pos = 2, offset = tl.offset, 
                ...)
        }
    }
    title(title, ...)
    if (!is.null(addCoef.col) && method != "number") {
        text(Pos[, 1], Pos[, 2], col = addCoef.col, labels = round((DAT - 
            int) * ifelse(addCoefasPercent, 100, 1)/zoom, number.digits), 
            cex = number.cex, font = number.font)
    }
    if (type == "full" && plotCI == "n" && !is.null(addgrid.col)) {
        rect(m1 - 0.5, n1 - 0.5, m2 + 0.5, n2 + 0.5, border = addgrid.col)
    }
    if (!is.null(addrect) && order == "hclust" && type == "full") {
        corrRect.hclust(corr, k = addrect, method = hclust.method, 
            col = rect.col, lwd = rect.lwd)
    }
    invisible(corr)
}

'order_samples' <- function(x) {
	res = rep(1, length(x))
	index = grep("VL", x, fixed=TRUE)
	res[index] = 2
	index = grep("VP", x, fixed=TRUE)
	res[index] = 3
	return(res)
}

'fun_lodmdl' <- function(df, mdl, grp, ...)
{
	model <- glm(call ~ expected_af, data=df, family = binomial(link = mdl))
	temp.data <- data.frame(expected_af = seq(0.01, max(df$expected_af), 0.01))
	predicted.data <- as.data.frame(predict(model, newdata = temp.data, se = TRUE))
	show.data <- cbind(temp.data, predicted.data) %>%
				 mutate(ymin = model$family$linkinv(fit - 1.96*se.fit),
				        ymax = model$family$linkinv(fit + 1.96*se.fit),
				        yfit = model$family$linkinv(fit),
				        group = grp)
  	return(show.data)
}

'fun_zerob' <- function(x, y, n=100, seed=0)
{
	set.seed(seed)
	init = data.frame(y, x)
	y0 = NULL
	for (i in 1:n) {
		index = sample(x=(1:nrow(init)), size=nrow(init), replace=TRUE)
		data = init[index,,drop=FALSE]
		z = try(zeroinfl(y ~ x, dist = "poisson", data = data), silent=TRUE)
		if ("try-error" %in% is(z)) {
			y0[[i]] = rep(NA, times=100)
		} else {
			x0 = data.frame(x=seq(20, 90, l=100))
			y0[[i]] = predict(object=z, newdata=x0, type = "count")
		}
	}
	y0 = do.call(cbind, y0)
	return(invisible(y0))
}

'box_plot' <- function (x, main = "", sub  = "", xlab = "", ylab = "", col, lwd  = 2, ...)
{
	if ("list" %in% is(x)) {
		nboxes = length(x)
		if (nboxes>10) {
			stop("no more than 10 boxes allowed")
		}
	} else if ("matrix" %in% is(x)) {
		nboxes = ncol(x)
		if (nboxes>10) {
			stop("no more than 10 boxes allowed")
		}
	}
	if (missing(col)) {
		col = vector(mode="character", length=nboxes)
		ucol = rainbow_hcl(nboxes, start=30, end=300)
		for (i in 1:nboxes) {
			col[i] = transparentRgb(col=ucol[i], alpha=126)
		}
	}
    boxplot(x, outline=FALSE, main="", xlab="", ylab="", axes=FALSE, frame=FALSE, lwd=1, ...)
    if (any(xlab=="")) {
    	axis(1, at=1:nboxes, labels=rep(xlab, nboxes), cex.axis=1.5, padj=0.25)
    } else {
    	axis(1, at=1:nboxes, labels=xlab, cex.axis=1.5, padj=0.25)
	}
    axis(2, at=NULL, cex.axis=1.5, las=1)
    mtext(side=2, text=ylab, line=3.5, cex=1.5)
    if (sub=="") {
    	title(main=main, cex.main=1.5)
    } else {
    	title(main=paste(main, "\n", sep=""), cex.main=2.0)
    	title(main=paste("\n", sub, sep=""), cex.main=1.5)
    }
    if ("list" %in% is(x)) {
		for (i in 1:nboxes) {
			points(jitter(rep(i, length(x[[i]])), amount=.25), x[[i]], pch=21, col = "black", bg=col[[i]], cex=2.0)
		}
    } else if ("matrix" %in% is(x)) {
		for (i in 1:nboxes) {
	    	points(jitter(rep(i, nrow(x)), amount=.25), x[,i], pch=21, col = "black", bg=col[[i]], cex=2.0)
	    }
	}
    box(lwd=lwd)
}
