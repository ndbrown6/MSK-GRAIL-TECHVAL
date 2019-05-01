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
		   mutate(facet = "Tumor matched") %>%
		   mutate(hgvsp_short = ifelse(gene == "PIK3CA" & hgvsp == "p.Glu542Lys", "E542K", "")) %>%
		   mutate(hgvsp_short = ifelse(gene == "PIK3CA" & hgvsp == "p.His1047Arg", "H1047R", hgvsp_short)) %>%
		   mutate(hgvsp_short = ifelse(gene == "EGFR" & hgvsp == "p.Leu861Gln", "L861Q", hgvsp_short)) %>%
		   mutate(hgvsp_short = ifelse(gene == "KRAS" & hgvsp == "p.Gly12Ala", "G12A", hgvsp_short)) %>%
		   mutate(hgvsp_short = ifelse(gene == "KRAS" & hgvsp == "p.Gly12Cys", "G12C", hgvsp_short)) %>%
		   mutate(UID = paste0(gene, " ", hgvsp_short))
		   
plot.0 = ggplot(combined, aes(y = af_grail, x = af_ddPCR, shape = Protocol, label = UID)) +
		 geom_abline(slope = 1, color = "goldenrod3", linetype = 1) +
		 geom_point(alpha = .8, size = 2.5, fill = "#2B83BA") +
		 scale_shape_manual(values = c(24, 21)) +
		 geom_text_repel(size = 2.25, color = "black", segment.color = "transparent") +
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
