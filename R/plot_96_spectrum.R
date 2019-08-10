'plot_96_spectrum' <- function(vcf,
							   sample.col = "sample",
							   mutcat3.col = "mutcat3",
							   ymax = NULL,
							   averageProp = FALSE,
							   file = NULL)
{
    bases = c("A", "C", "G", "T")
    ctxt16 = paste(rep(bases, each = 4), rep(bases, 4), sep = ".")
    mt = c("CA", "CG", "CT", "TA", "TC", "TG")
    types96 = paste(rep(mt, each = 16), rep(ctxt16, 6), sep = "_")
    types96 = sapply(types96, function(z) { sub("\\.", substr(z, 1, 1), z) })
    context = substr(types96, 4, 6)
    nsamp = length(unique(vcf[, sample.col]))
    if (averageProp & nsamp > 1) {
        tmp = makeMutypeMatFromVcf(vcf, sample.col = "CHCID", mutcat.col = "mutcat3", mutypes = types96)
        freq = apply(tmp, 1, mean)
    } else {
        freq = sapply(types96, function(z) { mean(vcf[, mutcat3.col] == z, na.rm = T) })
	}
    if (!is.null(file)) {
        pdf(file, width = 24, height = 5)
    }
    col96 = c(rep("skyblue3", 16), rep("black", 16), rep("red", 16), rep("grey", 16), rep("green", 16), rep("pink", 16))
    labs = c(rep("C>A", 16), rep("C>G", 16), rep("C>T", 16), rep("T>A", 16), rep("T>C", 16), rep("T>G", 16))
    if (is.null(ymax)) {
        ymax = ceiling(max(freq)*100)
        freq[freq>.4] = .4
    }
    if (ymax >= 30) {
    	ymax = 40
    	by = 10
	} else if (ymax >= 20) {
		ymax = 30
		by = 5
	} else if (ymax >= 10) {
		ymax = 20
		by = 5
	} else if (ymax < 10) {
		ymax = 10
		by = 2
	}
    bp = barplot(freq*100, col = col96, border = col96, las = 2, width = 1, space = .35, yaxt = "n", xaxt = "n", ylim = c(0, ymax * 1.2))
    title(ylab = "Fraction of\nmutations (%)", mgp = c(1, 1, 0), cex.lab = 1.85)
    axis(1, at = bp, labels = context, pos = 0, las = 2, cex.axis = 1.5, tick = F, cex.axis = 1, lwd=-1)
    axis(2, at = seq(0, ymax, by=by), labels=seq(0, ymax, by=by), pos = 0, las = 1, cex.axis = 1.75, lwd=1.5, lwd.ticks=1.35, line=3.5)
    for (i in seq(1, 81, by = 16)) {
        rect(bp[i], par()$usr[4], bp[i + 15], par()$usr[4] - 0.05 * diff(par()$usr[3:4]), col = col96[i], border = col96[i])
        text((bp[i] + bp[i + 15])/2, par()$usr[4] + 0.09 * diff(par()$usr[3:4]), labels = labs[i], xpd = TRUE, cex = 2)
    }
    if (!is.null(file)) {
        dev.off()
    }
    return(invisible(0))
}
