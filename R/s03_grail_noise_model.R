#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/figureS3")) {
	dir.create("../res/figureS3")
}

approx_posterior <- read_tsv(file="../modified_v11/Resources/Grail_noise_model/inputs/approx_posterior_by_alt.dispersion.tsv",
					col_types = cols(.default = col_character())) %>%
					type_convert()

manifest <- read_csv(file="../modified_v11/Tracker_files/s3_sample_tracker_TechVal_Merlin.csv",
			col_types = cols(.default = col_character())) %>%
			type_convert()
					
stacked_vcfs <- read_tsv(file="../modified_v11/Variants_Calls/Joined_cfDNA_IMPACT_variants/scored_merged_snvs_20171115.tsv",
				col_types = cols(.default = col_character())) %>%
				type_convert()
						
kMerlinContaminated <- data_frame(cfdna_sample_id = c(
  "MRL0200GH",
  "MRL0121CH",
  "MRL0122CH",
  "MRL0053CH",
  "MRL0081CH",
  "MRL0197GH_MRL0206GH",
  "MRL0055CH",
  "MRL0083CH",
  "MRL0089CH",
  "MRL0115CH",
  "MRL0235GH",
  "MRL0125CH",
  "SDBB-W044216563342-CH",
  "SDBB-W044216563343-CH",
  "MRL0154GH",
  "SDBB-W044216563331-CH",
  "MRL0067CH",
  "MRL0060CH",
  "MRL0062CH",
  "MRL0123CH",
  "MRL0057CH",
  "MRL0068CH",
  "MRL0245GH_MRL0144GH_MRL0215GH_MRL0244GH",
  "MRL0163GH",
  "MRL0198GH_MRL0207GH",
  "MRL0236GH",
  "MRL0234GH",
  "MRL0058CH",
  "MRL0233GH",
  "MRL0232GH",
  "MRL0222GH"
))

excl_patient_id <- manifest %>% 
                   filter(tissue %in% c("Healthy"), cfdna_sample_id %in% kMerlinContaminated$cfdna_sample_id) %>% 
                   .$patient_id
                  
unique_train_ids <- read_tsv(file="../modified_v11/Resources/Grail_noise_model/inputs/unique_train_ids.tsv",
					col_types = cols(.default = col_character()), col_names=FALSE) %>%
					mutate(cfdna_id = unlist(lapply(X1, function(i) strsplit(i, "[.]collapsed")[[1]][1]))) %>%
					mutate(patient_id = stringr::str_match(cfdna_id, "(BRH_[0-9]+){1}[-_].*")[, 2])

train_samples <- manifest %>% 
                 filter(cfdna_sample_id %in% unique_train_ids$cfdna_id)
                 
train_pileup_locs <- train_samples %>%
                     dplyr::select(patient_id, gdna_s3_pileup, cfdna_s3_pileup)
                     
cfdna_pileups <- read_tsv(file="../modified_v11/Resources/Grail_noise_model/intermediate/cfdna_pileups.tsv",
				 col_types = cols(.default = col_character())) %>%
				 type_convert()

cfdna_pileups_long <- cfdna_pileups %>%
                      dplyr::select(patient_id, CHROM, POS, DEPTH, REF, A, C, G, T, N) %>%
                      gather("ALT", "AD", -patient_id, -CHROM, -POS, -REF, -DEPTH) %>%
                      filter(!(ALT == "N")) %>%
                      filter(!(REF == ALT)) %>%
                      filter(!(AD/DEPTH > 0.3))
                        
sites_to_sample <- sample(unique(cfdna_pileups_long %>% 
                                 mutate(pos_id = str_c(CHROM, POS, REF, ALT, "_")) %>% 
                                 .$pos_id), 20000)
                                    
cfdna_pileups_noise_params <- cfdna_pileups_long %>%
                              mutate(pos_id = str_c(CHROM, POS, REF, ALT, "_")) %>%
                              filter(pos_id %in% sites_to_sample) %>%
                              dplyr::left_join(approx_posterior, by = c("CHROM", "POS", "REF", "ALT")) %>%
                              filter(!(is.na(mu)))

obs_error_rates <- cfdna_pileups_noise_params %>%
                   group_by(CHROM, POS, REF, ALT, mu, size) %>%
                   dplyr::summarise(tot_ad = sum(AD),
                                    tot_depth = sum(DEPTH),
                                    obs_mu = sum(AD)/sum(DEPTH),
                                    max_obs_mu = max(AD/DEPTH)) %>%
                   ungroup() %>%
                   mutate(tot_ad_range = case_when(tot_ad < 5 ~ as.character(tot_ad), tot_ad >= 5 ~ ">= 5")) %>%
                   mutate(tot_ad_range = factor(tot_ad_range, levels = c("0", "1", "2", "3", "4", ">= 5"))) %>%
                   mutate(facet = "Allele Frequency") %>%
                   mutate(obs_mu = ifelse(obs_mu==0, 1e-6, obs_mu))

pdf(file="../res/figureS3/Observed_vs_Posterior_Lambda_Mu.pdf", width=7.5, height=8)
par(mar = c(6.1, 7, 4.1, 1))
shapes = 21
cols = c("0"="salmon", "1"="#FDAE61", "2"="#ABDDA4", "3"="cadetblue", "4"="grey50", ">= 5"="steelblue")
plot(1, 1, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, xlim=c(1e-7,1), ylim=c(1e-6,1), log="xy")
axis(1, at=c(1e-7, 1e-5, 1e-3, 1e-1, 1), labels=c(expression(10^-7), expression(10^-5), expression(10^-3), 1e-1, 1), cex.axis = 1.75, padj = 0.25, lwd=1.85, lwd.ticks=1.75)
axis(2, at=c(1e-6, 1e-5, 1e-3, 1e-1, 1), labels=c(0, expression(10^-5), expression(10^-3), 1e-1, 1), cex.axis = 1.75, las = 1, lwd=1.85, lwd.ticks=1.75)
mtext(side = 1, text = expression("Mean posterior"~lambda[p]~"("~mu[p]~")"), line = 4, cex = 1.75)
mtext(side = 2, text = expression("Observed "~lambda[p]), line = 4.5, cex = 1.75)
points(c(1e-7,1), c(1e-7,1), type="l", lty=1, lwd=2, col="goldenrod3")
x = as.numeric(obs_error_rates$mu)
y = as.numeric(obs_error_rates$obs_mu)
z = as.character(obs_error_rates$tot_ad_range)
points(x, y, pch=shapes, col="black", bg=cols[z], cex=1.25)
log10_axis(side=1, at=c(1e-7, 1e-5, 1e-3, 1e-1, 1), lwd=0, lwd.ticks=1)
log10_axis(side=2, at=c(1e-6, 1e-5, 1e-3, 1e-1, 1), lwd=0, lwd.ticks=1)
legend(x=1e-7, y=.5, pch=shapes, col="black", pt.bg=cols, legend=names(cols), box.lwd=-1, pt.cex=1.35)
text(x=1e-7, y=.8, labels="Alternate\nallele count", pos=4, font=2)
text(x=1e-1, y=.8, labels="y = x", pos=4, col="goldenrod3", cex=1.25)
dev.off()

sites <- distinct(approx_posterior, CHROM, POS)

kRefGenomePath <- "../modified_v11/Resources/Grail_noise_model/inputs/genome.fa"

ref_genome <- Biostrings::readDNAStringSet(filepath=kRefGenomePath)

sites <- sites %>%
		 mutate(trinucleotide_context = as.character(Biostrings::subseq(ref_genome[CHROM], POS - 1, POS + 1)))

annotated <- approx_posterior %>%
			 left_join(sites) %>%
			 mutate(substitution_type = paste0("\n", REF, ">", ALT)) %>%
			 mutate(substitution_type = factor(substitution_type)) %>%
			 mutate(REF = factor(REF)) %>%
			 mutate(ALT = factor(ALT)) %>%
			 mutate(facet = "Allele Frequency")
			 
pdf(file="../res/figureS3/Posterior_Mu_by_Substitution_Type.pdf", width=6.5, height=7)
par(mar = c(6.1, 6, 4.1, 1))
shapes = 21
cols = c("\nA>C"="salmon", "\nA>G"="salmon", "\nA>T"="salmon",
		 "\nT>A"="#FDAE61", "\nT>C"="#FDAE61", "\nT>G"="#FDAE61",
		 "\nG>A"="#ABDDA4", "\nG>C"="#ABDDA4", "\nG>T"="#ABDDA4",
		 "\nC>A"="steelblue", "\nC>G"="steelblue", "\nC>T"="steelblue")
plot(1, 1, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, xlim=c(1,12), ylim=c(0,8))
for (i in 1:length(cols)) {
	z = mean_cl_normal((annotated$mu*1e5)[annotated$substitution_type==names(cols)[i]])
	points(c(i,i), z[1,2:3], type="l", col="black", lwd=1.5)
	points(i, z[1,1], pch=shapes, bg=cols[i], col="black", cex=1.25)
}
axis(1, at=1:length(cols), labels=names(cols), cex.axis = 1, las = 2, lwd=1.85, lwd.ticks=1.75)
axis(2, at=c(0, 2, 4, 6, 8), labels=c(0, 2, 4, 6, 8), cex.axis = 1.5, las = 1, lwd=1.85, lwd.ticks=1.75)
mtext(side = 2, text = expression("Mean posterior"~lambda[p]~"("~mu[p] %.% 10^5~")"), line = 3, cex = 1.5)
dev.off()

pdf(file="../res/figureS3/Posterior_Mu_by_Trinucleotide_Context.pdf", width=9, height=7)
par(mar = c(5.1, 7, 3.1, 1.1))
zz = split.screen(figs=matrix(c(0,.575, 0.425,1, 0.425,1,0.425,1, 0,.575,0,.575, .425,1,0,.575, 0,1,0,1) , nrow=5, ncol=4, byrow=TRUE))
shapes = 21
cols = c("A"="salmon",
		 "T"="#FDAE61",
		 "G"="#ABDDA4",
		 "C"="steelblue")
nuc = c("A", "T", "G", "C")
for (i in 1:length(nuc)) {
	screen(zz[i])
	plot(1, 1, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, xlim=c(1,16), ylim=c(0,40))
	for (j in 1:length(nuc)) {
		if (i!=j) {
			trinuc = sort(unique(annotated$trinucleotide_context[annotated$REF==nuc[i]]))
			for (ii in 1:length(trinuc)) {
				for (jj in 1:length(nuc)) {
					if (i!=jj) {
						z = mean_cl_normal((annotated$mu*1e5)[annotated$trinucleotide_context==trinuc[ii] & annotated$ALT==nuc[jj]])
						points(c(ii,ii), c(max(0,z[1,2]), z[1,3]), type="l", col="black", lwd=1.5)
					}
				}
				for (jj in 1:length(nuc)) {
					if (i!=jj) {
						z = mean_cl_normal((annotated$mu*1e5)[annotated$trinucleotide_context==trinuc[ii] & annotated$ALT==nuc[jj]])
						points(ii, z[1,1], pch=shapes, bg=cols[nuc[jj]], col="black", cex=1.25)
					}
				}
			}
		}
	}
	axis(1, at=1:16, labels=trinuc, cex.axis = 0.75, las = 2, lwd=1.5, lwd.ticks=1.25, line=0)
	if (i==1 | i==3) {
		axis(2, at=c(0, 10, 20, 30, 40), labels=c(0, 10, 20, 30, 40), cex.axis = 1.5, las = 1, lwd=1.5, lwd.ticks=1.25)
	} else {
		axis(2, at=c(0, 10, 20, 30, 40), labels=rep("", 5), cex.axis = 1.5, las = 1, lwd=1.5, lwd.ticks=1.25)
	}
	if (i==2) {
		text(x=13, y=38, labels="Alternate\nallele", pos=4, font=2, cex=.85)
		legend(x=13, y=37, pch=21, col="black", pt.bg=cols, legend=c("A", "T", "G", "C"), box.lwd=-1, cex=.85, pt.cex=1.15)
	}
}
screen(zz[5])
plot(1, 1, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, xlim=c(0,1), ylim=c(0,1))
mtext(side = 2, text = expression("Mean posterior"~lambda[p]~"("~mu[p] %.% 10^5~")"), line = 4, cex = 1.5)
close.screen(all.screens=TRUE)
dev.off()


scored_test_data <- read_tsv(file="../modified_v11/Resources/Grail_noise_model/inputs/scored_merged_snvs.tsv",
					col_types = cols(.default = col_character())) %>%
					type_convert()

test_samples <- scored_test_data %>% 
                filter(subj_type %in% c("Healthy")) %>%
                distinct(patient_id)
                
cfdna_pileups <- read_tsv(file="../modified_v11/Resources/Grail_noise_model/intermediate/cfdna_pileups.tsv",
				 col_types = cols(.default = col_character())) %>%
				 type_convert()
				 
cfdna_pileups_long <- cfdna_pileups %>%
					  dplyr::select(patient_id, CHROM, POS, DEPTH, REF, A, C, G, T, N) %>%
					  gather("ALT", "AD", -patient_id, -CHROM, -POS, -REF, -DEPTH) %>%
					  filter(!(ALT == "N")) %>%
					  filter(!(REF == ALT)) %>%
					  filter(AD > 0)

scored_pileups <- cfdna_pileups_long %>%
				  left_join(approx_posterior, by = c("CHROM", "POS", "REF", "ALT")) %>%
				  filter(!(is.na(mu))) %>%
				  filter(!(AD/DEPTH > 0.3)) %>%
				  filter(DEPTH > 200) %>%
				  filter(!(patient_id %in% excl_patient_id)) %>%
				  mutate(log_p = pnbinom(AD - 1,
				  						 size = size,
				  						 mu = mu * DEPTH,
				  						 lower.tail = FALSE,
				  						 log.p = TRUE)) %>%
				  mutate(score = -10 * log_p/log(10))

gdna_pileups <- read_tsv(file="../modified_v11/Resources/Grail_noise_model/intermediate/gdna_pileups.tsv",
				 col_types = cols(.default = col_character())) %>%
				 type_convert()
				 
gdna_params <- gdna_pileups %>%
			   dplyr::select(patient_id, CHROM, POS, DEPTH, REF, A, C, G, T, N) %>%
			   gather("ALT", "AD", -patient_id, -CHROM, -POS, -REF, -DEPTH) %>%
			   filter(!(ALT == "N")) %>%
			   filter(!(REF == ALT)) %>%
			   filter(AD > 0) %>%
			   inner_join(approx_posterior) %>%
			   mutate(gdna_size = AD + 1/2, gdna_mu = gdna_size / DEPTH) %>%
			   dplyr::select(-AD, -DEPTH)

test_scores <- scored_pileups %>%
               left_join(gdna_params) %>%
               mutate(gdna_score = -10 * pnbinom(AD - 1, size = gdna_size, mu = gdna_mu * DEPTH, lower.tail = FALSE, log.p = TRUE) / log(10)) %>%
               mutate(gdna_score = if_else(is.na(gdna_score), Inf, gdna_score)) %>%
               filter(gdna_score >= 60, AD / DEPTH > 3 * gdna_mu | is.na(gdna_mu))

kNumericalTol <- 1e-5
kPileupChrom <- "chr21"

panel_size <- approx_posterior %>%
			  filter(CHROM == kPileupChrom) %>%
			  distinct(CHROM, POS) %>%
			  n_distinct()

num_zeros <- test_scores %>%
			 group_by(patient_id) %>%
			 dplyr::summarize(num_zeros = panel_size - n()) %>%
			 ungroup() %>%
			 {sum(.$num_zeros)}

cdf <- ecdf(c(rep(.Machine$double.eps, num_zeros), test_scores$score))

calibration <- data_frame(Nominal = seq(0, 70, 1)) %>%
			   mutate(p_obs = 1 - cdf(Nominal)) %>%
			   mutate(Observed = -10 * log10(p_obs)) %>%
			   mutate(facet = "Calibration of quality scores")

plot.0 = ggplot(calibration,  aes(x = Nominal, y = Observed)) +
		 geom_line() +
		 geom_abline(color="goldenrod3", size=.5) +
		 facet_wrap(~facet) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=12)) +
		 labs(x="Nominal score (Phred-scaled)", y="Observed probability (Phred-scaled)\n\n") +
		 coord_cartesian(xlim=c(0,70), ylim = c(0,70)) +
		 theme(legend.justification = c(1, 0),
		 	   legend.position = c(0.2, 0.6),
		 	   legend.title = element_blank(),
		 	   legend.background = element_blank(),
		 	   legend.text=element_text(size=8))

pdf(file="../res/rebuttal/Nominal_vs_Observed_Pr.pdf", width=5.15, height=6)
print(plot.0)
dev.off()

all_snvs <- stacked_vcfs %>%
			dplyr::select(study, patient_id, chrom, position, ref, alt, gene, qualnobaq, isedge,
			is_snv, is_nonsyn, pgtkxgdna, MSK, grail, subj_type,
			clinicalaction, c_clinicalaction, ccd, c_panel, panel,
			adgdna, dpgdna) %>%
			filter(study == "TechVal" | subj_type == "Healthy") %>%
			mutate(is_driver = clinicalaction != "" | c_clinicalaction != "",
				   panel_overlap = (c_panel == 1 | panel == 1),
				   use_for_concordance = ccd == 1 & panel_overlap,
				   is_germline = adgdna / dpgdna > 0.2)

cancer_snvs <- all_snvs %>%
			   filter(subj_type != "Healthy")

qualnobaq_cutoffs <- seq(0, 3000, 5)
joint_cutoffs <- seq(0, 1, 0.01)
fixed_p <- 0.8
fixed_qualnobaq <- 60 

cancer_types = unique(cancer_snvs$subj_type)
roc_qual <- foreach(i = 1:length(qualnobaq_cutoffs),
                    .combine = rbind,
                    .packages = c('tidyverse', 'data.table')) %dopar% {
                    	qualnobaq_min <- qualnobaq_cutoffs[i]
                    	tissue_compare <- cancer_snvs %>%
                        				  filter(use_for_concordance, !patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id)) %>% 
                        				  group_by(subj_type) %>%
                        				  mutate(qual_edge_pass = (isedge != TRUE) & qualnobaq >= qualnobaq_min,
                        				  		 called = grail & !is_germline & qual_edge_pass & (pgtkxgdna >= fixed_p),
                        				  		 is_tp = called & MSK) %>%
                        				  summarize(total_pos = sum(MSK, na.rm = TRUE),
                        				  			total_called = sum(called, na.rm = TRUE),
                        				  			tp = sum(is_tp, na.rm = TRUE),
                        				  			recall = tp / total_pos,
                        				  			precision = tp / total_called) %>%
                        				  mutate(qualnobaq_min = qualnobaq_min)}

highlight_q60 <- roc_qual %>%
				 filter(abs(qualnobaq_min - 60) < kNumericalTol)

roc_joint <- foreach(i = 1:length(joint_cutoffs),
                     .combine = rbind,
                     .packages = c('tidyverse', 'data.table')) %dopar% {
                     	joint_min <- joint_cutoffs[i]
                     	tissue_compare <- cancer_snvs %>%
                     					  filter(use_for_concordance, !patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id)) %>%
                     					  group_by(subj_type) %>%
                     					  mutate(qual_edge_pass = (isedge != TRUE) & qualnobaq >= fixed_qualnobaq,
												 called = grail & !is_germline & qual_edge_pass & (pgtkxgdna >= joint_min),
												 is_tp = called & MSK) %>%
										  summarize(total_pos = sum(MSK, na.rm = TRUE),
										  			total_called = sum(called, na.rm = TRUE),
										  			tp = sum(is_tp, na.rm = TRUE),
										  			recall = tp / total_pos,
										  			precision = tp / total_called) %>%
										  mutate(joint_min = joint_min)}
										  
highlight_p <- roc_joint %>%
			   filter((subj_type == "Breast" & abs(joint_min - 0.79) < kNumericalTol) |
					  (subj_type == "Lung" & abs(joint_min - 0.82) < kNumericalTol) |
					  (subj_type == "Prostate" & abs(joint_min - 0.79) < kNumericalTol))

noncancer_snvs <- all_snvs %>%
				  filter(subj_type == "Healthy")

n_noncancer <- length(unique(noncancer_snvs$patient_id))
roc_noncancer_qual <- foreach(i = 1:length(qualnobaq_cutoffs),
							  .combine = rbind,
                              .packages = c('tidyverse', 'data.table')) %dopar% {
                              	qualnobaq_min <- qualnobaq_cutoffs[i]
                                fp <- noncancer_snvs %>%
                                	  filter(panel_overlap, !patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id)) %>%
                                	  mutate(qual_edge_pass = (isedge != TRUE) & qualnobaq >= qualnobaq_min,
                                	  		 called = grail & !is_germline & qual_edge_pass & (pgtkxgdna >= fixed_p)) %>%
                                	  summarize(sample_mean = sum(called, na.rm = TRUE) / n_noncancer) %>%
                                	  mutate(qualnobaq_min = qualnobaq_min)}
                               
noncancer_snvs <- all_snvs %>%
				  filter(subj_type == "Healthy")

roc_noncancer_joint <- foreach(i = 1:length(joint_cutoffs),
                               .combine = rbind,
                               .packages = c('tidyverse', 'data.table')) %dopar% {
                               		joint_min <- joint_cutoffs[i]
                                 	fp <- noncancer_snvs %>%
                                 		  filter(panel_overlap, !patient_id %in% c(hypermutators$patient_id, msi_hypermutators$patient_id)) %>%
                                 		  mutate(qual_edge_pass = (isedge != TRUE) & qualnobaq >= fixed_qualnobaq,
                                          		 called = grail & !is_germline & qual_edge_pass & (pgtkxgdna >= joint_min)) %>%
                                          		 summarize(sample_mean = sum(called, na.rm = TRUE) / n_noncancer) %>%
                                          		 mutate(joint_min = joint_min)}

recall_FP_qual <- roc_qual %>%
				  inner_join(roc_noncancer_qual, by = "qualnobaq_min")

highlight_recall_FP_q60 <- recall_FP_qual %>%
						   filter(abs(qualnobaq_min - 60) < kNumericalTol)

plot.0 = ggplot(recall_FP_qual %>% mutate(facet = "Recall rate"), aes(x = sample_mean, y = recall, color=subj_type)) +
		 geom_line(size=1) +
		 scale_color_manual(values=c("salmon", "#FDAE61", "#ABDDA4")) +
		 facet_wrap(~facet) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=12)) +
		 labs(x="Mean number of SNVs / sample", y="Recall rate\n\n") +
         coord_cartesian(xlim=c(0, 7), ylim = c(0.2, 1)) +
         geom_point(data = highlight_recall_FP_q60,
             aes(x = sample_mean, y = recall, colour = factor(subj_type)),
             show.legend = FALSE) +
         geom_label(data = highlight_recall_FP_q60,
             aes(label = "QUALNOBAQ >= 60"),
             nudge_x = 0.1, nudge_y = -0.1, size = 3) +
         guides(color=guide_legend(title=c("Cancer\ntype"))) +
         theme(legend.justification = c(1, 0),
		 	   legend.position = c(1, 0),
		 	   legend.title = element_blank(),
		 	   legend.background = element_blank(),
		 	   legend.text=element_text(size=8))
             
pdf(file="../res/rebuttal/Number_SNVs_Recall.pdf", width=5.15, height=6)
print(plot.0)
dev.off()

recall_FP_joint <- roc_joint %>%
				   inner_join(roc_noncancer_joint, by = "joint_min")

highlight_recall_FP_joint <- recall_FP_joint %>%
							 filter((subj_type == "Breast" & abs(joint_min - 0.79) < kNumericalTol) |
							 		(subj_type == "Lung" & abs(joint_min - 0.82) < kNumericalTol) |
							 		(subj_type == "Prostate" & abs(joint_min - 0.79) < kNumericalTol))
							 		
							 		
plot.0 = ggplot(recall_FP_joint %>% filter(subj_type=="Breast"), aes(x = sample_mean, y = recall)) +
		 geom_line(size=1, color="salmon") +
		 facet_wrap(~subj_type) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=12)) +
		 labs(x="Mean number of SNVs / sample", y="Recall rate\n\n") +
         coord_cartesian(xlim=c(0, 7), ylim = c(0.2, 1)) +
         geom_point(data = highlight_recall_FP_joint %>% filter(subj_type=="Breast"),
         		aes(x = sample_mean, y = recall),
         		show.legend = FALSE, color="salmon") +
         geom_label(data = highlight_recall_FP_joint %>% filter(subj_type=="Breast"),
         		aes(label = paste("pgtkxgdna >=", joint_min, sep = " ")),
         		nudge_x = 0.1, nudge_y = -0.1, size = 3, color="salmon")
         

pdf("../res/rebuttal/Number_SNVs_Recall_Pr_WBC_Breast.pdf", width=5, height=6)
print(plot.0)
dev.off()

plot.0 = ggplot(recall_FP_joint %>% filter(subj_type=="Lung"), aes(x = sample_mean, y = recall)) +
		 geom_line(size=1, color="#FDAE61") +
		 facet_wrap(~subj_type) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=12)) +
		 labs(x="Mean number of SNVs / sample", y="Recall rate\n\n") +
         coord_cartesian(xlim=c(0, 7), ylim = c(0.2, 1)) +
         geom_point(data = highlight_recall_FP_joint %>% filter(subj_type=="Lung"),
         		aes(x = sample_mean, y = recall),
         		show.legend = FALSE, color="#FDAE61") +
         geom_label(data = highlight_recall_FP_joint %>% filter(subj_type=="Lung"),
         		aes(label = paste("pgtkxgdna >=", joint_min, sep = " ")),
         		nudge_x = 0.1, nudge_y = -0.1, size = 3, color="#FDAE61")
         

pdf("../res/rebuttal/Number_SNVs_Recall_Pr_WBC_Lung.pdf", width=5, height=6)
print(plot.0)
dev.off()


plot.0 = ggplot(recall_FP_joint %>% filter(subj_type=="Prostate"), aes(x = sample_mean, y = recall)) +
		 geom_line(size=1, color="#ABDDA4") +
		 facet_wrap(~subj_type) +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=12)) +
		 labs(x="Mean number of SNVs / sample", y="Recall rate\n\n") +
         coord_cartesian(xlim=c(0, 7), ylim = c(0.2, 1)) +
         geom_point(data = highlight_recall_FP_joint %>% filter(subj_type=="Prostate"),
         		aes(x = sample_mean, y = recall),
         		show.legend = FALSE, color="#ABDDA4") +
         geom_label(data = highlight_recall_FP_joint %>% filter(subj_type=="Prostate"),
         		aes(label = paste("pgtkxgdna >=", joint_min, sep = " ")),
         		nudge_x = 0.1, nudge_y = -0.1, size = 3, color="#ABDDA4")
         

pdf("../res/rebuttal/Number_SNVs_Recall_Pr_WBC_Prostate.pdf", width=5, height=6)
print(plot.0)
dev.off()