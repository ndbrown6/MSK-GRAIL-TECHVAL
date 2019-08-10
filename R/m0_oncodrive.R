'm0_oncodrive' <- function(maf,
						   AACol = NULL,
						   minMut = 5,
						   pvalMethod = "zscore",
						   nBgGenes = 100,
						   bgEstimate = TRUE,
						   ignoreGenes = NULL) 
{
	gl = read.csv(system.file("extdata", "prot_len.txt.gz", package = "maftools"), sep = "\t", header=TRUE, stringsAsFactors = FALSE)
	pval.options = c("zscore", "poisson", "combined")
	if (!pvalMethod %in% pval.options) {
		stop("pvalMethod can only be either zscore, poisson or combined")
	}
	if (length(pvalMethod) > 1) {
		stop("pvalMethod can only be either zscore, poisson or combined")
	}
	syn.maf = maf@maf.silent
	numSamples = as.numeric(maf@summary[3, summary])
	if (bgEstimate) {
		if (nrow(syn.maf) == 0) {
			message("No syn mutations found! Skipping background estimation. Using predefined values. (Mean = 0.279; SD = 0.13)")
			bg.mean = 0.279
			bg.sd = 0.13
		} else {
			message("Estimating background scores from synonymous variants..")
			syn.bg.scores = maftools:::parse_prot(dat = syn.maf, AACol = AACol, gl, m = minMut, calBg = TRUE, nBg = nBgGenes)
			if (is.null(syn.bg.scores)) {
	            message("Not enough genes to build background. Using predefined values. (Mean = 0.279; SD = 0.13)")
	            bg.mean = 0.279
	            bg.sd = 0.13
	        } else {
	        	if (nrow(syn.bg.scores) < nBgGenes) {
	        		message("Not enough genes to build background. Using predefined values. (Mean = 0.279; SD = 0.13)")
	                bg.mean = 0.279
	                bg.sd = 0.13
	            } else {
	                bg.mean = mean(syn.bg.scores$clusterScores)
	                bg.sd = sd(syn.bg.scores$clusterScores)
	                message(paste("Estimated background mean: ", bg.mean))
	                message(paste("Estimated background SD: ", bg.sd))
	            }
	        }
	    }
	} else {
		message("Using predefined values for background. (Mean = 0.279; SD = 0.13)")
        bg.mean = 0.279
        bg.sd = 0.13
    }
    non.syn.maf = maf@data
    if (!is.null(ignoreGenes)) {
        ignoreGenes.count = nrow(non.syn.maf[Hugo_Symbol %in% ignoreGenes])
        message(paste("Removed", ignoreGenes.count, "variants belonging to", paste(ignoreGenes, collapse = ", ", sep = ",")))
        non.syn.maf = non.syn.maf[!Hugo_Symbol %in% ignoreGenes]
    }
    message("Estimating cluster scores from non-syn variants..")
    nonsyn.scores = maftools:::parse_prot(dat = non.syn.maf, AACol = AACol, gl = gl, m = minMut, calBg = FALSE, nBg = nBgGenes)
    if (pvalMethod == "combined") {
        message("Comapring with background model and estimating p-values..")
        nonsyn.scores$zscore = (nonsyn.scores$clusterScores - 
            bg.mean)/bg.sd
        nonsyn.scores$tPval = 1 - pnorm(nonsyn.scores$zscore)
        nonsyn.scores$tFdr = p.adjust(nonsyn.scores$tPval, method = "fdr")
        nonsyn.scores = merge(getGeneSummary(maf), nonsyn.scores, 
            by = "Hugo_Symbol")
        nonsyn.scores[, `:=`(fract_muts_in_clusters, muts_in_clusters/total)]
        counts.glm = glm(formula = total ~ protLen + clusters, 
            family = poisson(link = identity), data = nonsyn.scores)
        nonsyn.scores$Expected = counts.glm$fitted.values
        observed_mut_colIndex = which(colnames(nonsyn.scores) == 
            "total")
        expected_mut_colIndex = which(colnames(nonsyn.scores) == 
            "Expected")
        nonsyn.scores$poissonPval = apply(nonsyn.scores, 1, function(x) {
            poisson.test(as.numeric(x[observed_mut_colIndex]), 
                as.numeric(x[expected_mut_colIndex]))$p.value
        })
        nonsyn.scores$poissonFdr = p.adjust(nonsyn.scores$poissonPval, 
            method = "fdr")
        nonsyn.scores = nonsyn.scores[order(poissonFdr)]
        nonsyn.scores$fdr = apply(nonsyn.scores[, .(tFdr, poissonFdr)], 
            MARGIN = 1, FUN = min)
    } else if (pvalMethod == "zscore") {
        message("Comapring with background model and estimating p-values..")
        nonsyn.scores$zscore = (nonsyn.scores$clusterScores - 
            bg.mean)/bg.sd
        nonsyn.scores$pval = 1 - pnorm(nonsyn.scores$zscore)
        nonsyn.scores$fdr = p.adjust(nonsyn.scores$pval, method = "fdr")
        nonsyn.scores = merge(getGeneSummary(maf), nonsyn.scores, 
            by = "Hugo_Symbol")
        nonsyn.scores[, `:=`(fract_muts_in_clusters, muts_in_clusters/total)]
        nonsyn.scores = nonsyn.scores[order(fdr)]
    } else {
        nonsyn.scores = merge(getGeneSummary(maf), nonsyn.scores, 
            by = "Hugo_Symbol")
        nonsyn.scores[, `:=`(fract_muts_in_clusters, muts_in_clusters/total)]
        counts.glm = glm(formula = total ~ protLen + clusters, 
            family = poisson(link = identity), data = nonsyn.scores)
        nonsyn.scores$Expected = counts.glm$fitted.values
        observed_mut_colIndex = which(colnames(nonsyn.scores) == 
            "total")
        expected_mut_colIndex = which(colnames(nonsyn.scores) == 
            "Expected")
        nonsyn.scores$pval = apply(nonsyn.scores, 1, function(x) {
            poisson.test(as.numeric(x[observed_mut_colIndex]), 
                as.numeric(x[expected_mut_colIndex]))$p.value
        })
        nonsyn.scores$fdr = p.adjust(nonsyn.scores$pval, method = "fdr")
        nonsyn.scores = nonsyn.scores[order(fdr)]
    }
    message("Done !")
    return(invisible(nonsyn.scores))
}
