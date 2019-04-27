#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/rebuttal")) {
	dir.create("../res/rebuttal")
}

#==================================================
# Scatter plot of replicates and ddPCR
#==================================================
ddpcr_manifest = data.frame(
  patient_id = c("MSK-VB-0050",
                 "MSK-VL-0028",
                 "MSK-VL-0042", 
                 "MSK-VB-0023", 
                 "MSK-VL-0038"),
  gene = c("PIK3CA", 
           "EGFR", 
           "KRAS", 
           "PIK3CA", 
           "KRAS"),
  hgvsp  = c("p.His1047Arg", 
             "p.Leu861Gln", 
             "p.Gly12Cys", 
             "p.Glu542Lys", 
             "p.Gly12Ala"))

ddpcr = gs_url(url_ddpcr) %>%
		gs_read(ws="data", col_names = T, range = cell_rows(1:11)) %>%
		filter(!is.na(DNA))
retest = read_tsv(url_retest)
original = read_tsv(url_original)
ddpcr_filtered = ddpcr %>%
				 select(patient_id = X3, probe, af_ddPCR = `% MT`) %>%
				 full_join(ddpcr_manifest)

retest_filtered = retest %>%
				  filter(patient_id %in% ddpcr_manifest$patient_id & genename %in% ddpcr_manifest$gene & hgvsp %in% ddpcr_manifest$hgvsp) %>%
				  select(sample_id, patient_id, gene = genename, hgvsp, adnobaq, dpnobaq, qualnobaq, isedge, pgtkxgdna) %>%
				  mutate(af_grail = adnobaq/dpnobaq*100) %>%
				  mutate(Protocol = "V2")

original_filtered = original %>%
					filter(patient_id %in% ddpcr_manifest$patient_id & gene_id %in% ddpcr_manifest$gene & hgvs_p %in% ddpcr_manifest$hgvsp) %>%
					select(sample_id, patient_id, gene = gene_id, hgvsp = hgvs_p, adnobaq, dpnobaq, qualnobaq, isedge, pgtkxgdna) %>%
					mutate(af_grail = adnobaq/dpnobaq*100) %>%
					mutate(Protocol = "V1")

grail = full_join(original_filtered, retest_filtered)  
combined = left_join(grail, ddpcr_filtered) %>%
		   filter(qualnobaq >= 60) %>%
		   mutate(facet = "Tumor matched") #%>%
		   #mutate(UID = )
		   
plot.0 = ggplot(combined, aes(y = af_grail, x = af_ddPCR, shape = Protocol)) +
		 geom_abline(slope = 1, color = "goldenrod3", linetype = 1) +
		 geom_point(alpha = .8, size = 2.5, fill = "#2B83BA") +
		 scale_shape_manual(values = c(24, 21)) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="\nTargeted DNA assay (%)\n", x="\nddPCR (%)\n") +
		 scale_x_log10(
		 	breaks = function(x) { c(.05, .1, .2, .5, 1, 5, 10, 50) },
 			labels = function(x) { c(".05",".1",".2",".5","1","5","10","50") }
		 ) +
		 scale_y_log10(
		 	breaks = function(x) { c(.05, .1, .2, .5, 1, 5, 10, 50) },
 			labels = function(x) { c(".05",".1",".2",".5","1","5","10","50") }
		 ) +
		 annotation_logticks() +
		 coord_cartesian(xlim=c(0.05,50), ylim = c(0.05,50)) +
		 facet_wrap(~facet) +
		 guides(shape=guide_legend(title=c("Assay protocol"), override.aes=list(fill="black")))
		 
pdf(file="../res/rebuttal/ddPCR_Biopsy_Matched.pdf", width=6, height=6)
print(plot.0)
dev.off()


pdf(file="../res/figure1/assay_reproduce_combined_ddprc.pdf", width=7, height=7)
par(mar = c(6.1, 6, 4.1, 1))
shapes = c("V1"=21, "V2"=24)
plot(1, 1, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, xlim=c(0.05, 50), ylim=c(0.05, 50), log="xy")
points(c(0.05,45), c(0.05,45), type="l", lty=1, lwd=2, col="goldenrod3")
points(combined$af_ddPCR, combined$af_grail, pch=shapes[combined$Protocol], cex=1.5)
axis(1, at = c(.05,.1,.2,.5,1,5,10,50), labels=c(".05",".1",".2",".5","1","5","10","50"), cex.axis = 1.5, padj = 0.25, lwd=1.25, lwd.ticks=1.35)
axis(2, at = c(.05,.1,.2,.5,1,5,10,50), labels=c(".05",".1",".2",".5","1","5","10","50"), cex.axis = 1.5, las = 1, lwd=1.5, lwd.ticks=1.35)
mtext(side = 1, text = "ddPCR (%)", line = 4, cex = 1.35)
mtext(side = 2, text = "Targeted DNA assay (%)", line = 4, cex = 1.35)
mtext(side=3, text="Assay protocol", line=-1.15, at=.12, font=2)
legend(x=.05, y=45, pch=24, col="black", pt.bg="white", pt.cex=1.5, legend="V1", box.lwd=-1, cex=1)
legend(x=.05, y=30, pch=21, col="black", pt.bg="white", pt.cex=1.5, legend="V2", box.lwd=-1, cex=1)
text(x=.42, y=.1, labels="KRAS G12A")
text(x=.44, y=.18, labels="KRAS G12A")
text(x=.32, y=.9, labels="EGFR L861Q")
text(x=8.4, y=3.9, labels="KRAS G12C")
text(x=7, y=23, labels="PIK3CA H1047R")
text(x=12, y=35, labels="PIK3CA E542K")
dev.off()
