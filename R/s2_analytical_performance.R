#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
source('config.R')

if (!dir.exists("../res/figureS2")) {
	dir.create("../res/figureS2")
}

#==================================================
# Scatter plot of analytical performance
#==================================================
cell_line = read_tsv(url_cell.line, col_types = NULL) %>%
		   	type_convert()
hd753_joint = cell_line %>%
			   filter(startsWith(patient_id, "HD753") | startsWith(patient_id, "NA12878"))
hd753_truth = read_tsv(url_HD753.truth, col_types = NULL) %>%
		   	  type_convert()
hd753_truth_join = hd753_truth %>%
				   mutate(hgvsp_chk = paste0(gene, ":", hgvsp)) %>%
				   dplyr::select(chrom, hgvsp_chk, snv, indel, stock_af)
dilution = hd753_joint %>%
		   dplyr::select(patient_id) %>%
		   distinct() %>%
		   mutate(dilution_factor = 
		   			ifelse(grepl("HD753_0_1", patient_id), 0.02,
		   				ifelse(grepl("HD753_0_25", patient_id), 0.05,
		   					ifelse(grepl("HD753_0_5", patient_id), 0.1,
		   						ifelse(grepl("HD753_1", patient_id), 0.2, 0))))) %>%
		   mutate(group = ifelse(grepl("LOD", patient_id), "synthetic", "real"))
		   
hd753_truth_join = hd753_truth_join[rep(seq_len(nrow(hd753_truth_join)), each = length(unique(dilution$dilution_factor))), ]
hd753_truth_join$dilution_factor = rep(unique(dilution$dilution_factor), len = nrow(hd753_truth_join))
hd753_truth_join = full_join(dilution, hd753_truth_join) %>%
				   dplyr::select(patient_id, chrom, hgvsp_chk, snv, indel, stock_af, dilution_factor, group) %>%
				   mutate(nominal_af = ifelse(stock_af < 6, mean(stock_af[stock_af <6])*dilution_factor, mean(stock_af[stock_af >=6])*dilution_factor))
hd753_joint.filtered = hd753_joint %>%
					   mutate(hgvsp_chk = vapply(strsplit(hgvsp, "\\|"), `[`, 1, FUN.VALUE=character(1))) %>%
					   filter(hgvsp_chk %in% hd753_truth_join$hgvsp_chk) %>%
					   full_join(hd753_truth_join) %>%
					   mutate(expected_af = stock_af * dilution_factor) %>%
					   mutate(freq = ifelse(is.na(freq), 0, freq),
					   					filter = ifelse(qualnobaq < 60 | pgtkxgdna < 0.8 | isedge | is.na(qualnobaq), "FAIL", "PASS"), call = ifelse(filter == "PASS", 1L, 0L))
sensitivity_real = hd753_joint.filtered %>%
				   filter(expected_af < 3 & group == "real")
sensitivity_syn = hd753_joint.filtered %>%
				  filter(expected_af < 3 & group == "synthetic")

regmethod = "probit"
show.data.real = fun_lodmdl(sensitivity_real, regmethod, "real")
lod95_real = show.data.real$expected_af[max(which(show.data.real$yfit <= 0.95))]
show.data.syn = fun_lodmdl(sensitivity_syn, regmethod, "synthetic")
lod95_syn = show.data.syn$expected_af[max(which(show.data.syn$yfit <= 0.95))]
show.data = full_join(show.data.real, show.data.syn)
sensitivity = full_join(sensitivity_real, sensitivity_syn)
show.data = show.data %>% filter(expected_af<1.2)

cols = c("#D7191C", "#2B83BA")
pdf(file="../res/figureS2/analytical_performance.pdf", width=7, height=7)
par(mar = c(6.1, 6, 4.1, 1))
plot(1,1, type="n", xlim=c(0,1.2), ylim=c(0,1), xlab="", ylab="", axes=FALSE, frame.plot=FALSE)
polygon(x=c(show.data[show.data$group=="real","expected_af"], rev(show.data[show.data$group=="real","expected_af"])),
		y=c(show.data[show.data$group=="real","ymax"], rev(show.data[show.data$group=="real","ymin"])),
		border = cols[1], col = transparentRgb(cols[1], 55), lwd=.5)
points(show.data[show.data$group=="real","expected_af"], show.data[show.data$group=="real","yfit"], type="l", col=cols[1], lwd=1.5)
points(show.data[show.data$group=="synthetic","expected_af"], show.data[show.data$group=="synthetic","yfit"], type="l", col=cols[2], lwd=1.5)
polygon(x=c(show.data[show.data$group=="synthetic","expected_af"], rev(show.data[show.data$group=="synthetic","expected_af"])),
		y=c(show.data[show.data$group=="synthetic","ymax"], rev(show.data[show.data$group=="synthetic","ymin"])),
		border = cols[2], col = transparentRgb(cols[2], 55), lwd=.5)
x = jitter(sensitivity$expected_af[sensitivity$group=="real" & !sensitivity$isedge], factor=.15)
y = jitter(sensitivity$call[sensitivity$group=="real" & !sensitivity$isedge], factor=.1)
y = ifelse(y>1, 1, ifelse(y<0, 0, y))
points(x, y, type="p", pch=1, cex=.85, col=transparentRgb(cols[1], 185), lwd=.65)
x = jitter(sensitivity$expected_af[sensitivity$group=="synthetic" & !sensitivity$isedge], factor=.15)
y = jitter(sensitivity$call[sensitivity$group=="synthetic" & !sensitivity$isedge], factor=.1)
y = ifelse(y>1, 1, ifelse(y<0, 0, y))
points(x, y, type="p", pch=1, cex=.85, col=transparentRgb(cols[2], 185), lwd=.65)
points(x=rep(lod95_real,2), y=c(-1,1.01), col=cols[1], type="l", lwd=1.5, lty=3)
points(x=rep(lod95_syn,2), y=c(-1,1.01), col=cols[2], type="l", lwd=1.5, lty=3)
legend(x=.97, y=0.195, col=cols, legend=c("2430X", "4577X\n(synthetic)"), title="Mean collapsed\ntarget coverage", box.lwd=-1, lty=1, pch=1, pt.cex=.75, cex=.9)
axis(1, at = NULL, cex.axis = 1.5, padj = 0.25, lwd=1.5, lwd.ticks=1.35)
axis(2, at = NULL, cex.axis = 1.5, las = 1, lwd=1.5, lwd.ticks=1.35)
mtext(side = 1, text = "Variant allele frequency (%)", line = 4, cex = 1.5)
mtext(side = 2, text = "Detection probability", line = 4, cex = 1.85)
dev.off()
