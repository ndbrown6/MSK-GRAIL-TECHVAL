#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
source('config.R')

if (!dir.exists("../res/figureS7")) {
	dir.create("../res/figureS7")
}

#==================================================
# bar plot of msi scores
#==================================================
all_msi_df = read_tsv(file=url_msi_processed_data)

alias_df = read_tsv(file=url_subject_alias)

type_frame_df = data_frame(type=c("VB","VL","VP"), expand=c("Metastatic Breast","Metastatic Lung","Metastatic Prostate"), subject_type=c("Breast","Lung","Prostate"))

useful_msi_df = all_msi_df %>%
				left_join(alias_df %>% mutate(id=patient_id)) %>%
				mutate(id=alias) %>%
				dplyr::select(patient_id,original_msi,fixed_msi,subject_type,sample)
				
cols = c("Tumor"="#005A91",
		 "cfDNA"="#65A9D3")
		 
msi_tumor = useful_msi_df %>%
		    filter(sample=="tumor") %>%
		    arrange(patient_id)
msi_cfdna = useful_msi_df %>%
		    filter(sample=="cfdna") %>%
		    arrange(patient_id)

index = order(msi_tumor$original_msi)
msi_tumor = msi_tumor[index,,drop=FALSE]
msi_cfdna = msi_cfdna[index,,drop=FALSE]
msi_tumor = bind_cols(msi_tumor, "level" = factor(msi_tumor$patient_id, levels=msi_tumor$patient_id, ordered=TRUE))
msi_cfdna = bind_cols(msi_cfdna, "level" = factor(msi_cfdna$patient_id, levels=msi_cfdna$patient_id, ordered=TRUE))
msi_df = bind_rows(msi_tumor, msi_cfdna)
		 
plot.0 = msi_df %>%
		 mutate(sample = ifelse(sample=="cfdna", "cfDNA", "Tumor")) %>%
		 ggplot(aes(x=level, y=original_msi+.001, group=sample, fill=sample)) +
		 geom_bar(stat="identity", position=position_dodge(), color="black", width=.8) +
		 facet_wrap(~subject_type, scales = "free_y") +
		 theme_classic(base_size=15) +
		 coord_flip() +
		 labs(y="\nMSI score\n", x="\n") +
		 scale_fill_manual(values = cols) +
		 guides(fill=guide_legend(title=c("Sample type"))) +
		 theme(legend.title=element_text(size=14)) +
		 ylim(0, 0.3) +
		 theme(
		 	axis.title.x = element_text(size = 16),
		 	axis.text.x = element_text(size = 16),
		 	strip.background = element_rect(colour = "white", fill = "white"),
		 	strip.text.x = element_text(size = 16))

pdf(file="../res/figureS7/msi_score_original.pdf", height=12, width=16)
h = matrix(NA, nrow=length(unique(msi_df %>% .[["patient_id"]])), ncol=3, dimnames=list(unique(msi_df %>% .[["patient_id"]]), c("Tumor", "cfDNA", "subj_type")))
h[,"Tumor"] = msi_df %>%
			  filter(sample=="tumor") %>%
			  .[["original_msi"]]
h[,"cfDNA"] = msi_df %>%
			  filter(sample=="cfdna") %>%
			  .[["original_msi"]]
subj_type = msi_df %>%
			filter(sample=="tumor") %>%
			.[["subject_type"]]
z = split.screen(matrix(c(0+.05+.025,1/3+.025,0,1, 1/3+.025+.025,2/3-.025+.025,0,1, 2/3+.025,1-.05+.025,0,1), nrow=3, ncol=4, byrow=TRUE))
screen(z[1])
x = barplot(t(subset(h[,c(2:1),drop=FALSE], subj_type=="Breast")), horiz=TRUE, las=1, beside=TRUE,
		col=rev(cols), xlim=c(0,.3), cex.axis = 1.25, lwd=1.85, space=c(0,.3))
axis(side=2, at=apply(x, 2, mean), labels=rep("", ncol(x)), line=0.35, tcl=-.25, lwd=1.35)
title(main="\n\nBreast", cex.main=1.15)
screen(z[2])
x = barplot(t(subset(h[,c(2:1),drop=FALSE], subj_type=="Lung")), horiz=TRUE, las=1, beside=TRUE,
		col=rev(cols), xlim=c(0,.3), cex.axis = 1.25, lwd=1.85, space=c(0,.3))
axis(side=2, at=apply(x, 2, mean), labels=rep("", ncol(x)), line=0.35, tcl=-.25, lwd=1.35)
title(main="\n\nLung", xlab="\n\n\n\nMSI score", cex.main=1.15, cex.lab=1.5)
screen(z[3])
x = barplot(t(subset(h[,c(2:1),drop=FALSE], subj_type=="Prostate")), horiz=TRUE, las=1, beside=TRUE,
		col=rev(cols), xlim=c(0,.3), cex.axis = 1.25, lwd=1.85, space=c(0,.3))
axis(side=2, at=apply(x, 2, mean), labels=rep("", ncol(x)), line=0.35, tcl=-.25, lwd=1.35)
title(main="\n\nProstate", cex.main=1.15)
legend(x=0.2, y=70, pch=22, col="black", legend=c("Tumor", "cfDNA"), box.lwd=-1, cex=1.15, pt.cex=2.15, pt.bg=rev(cols))
close.screen(all.screens=TRUE)
dev.off()


plot.0 = msi_df %>%
		 mutate(sample = ifelse(sample=="cfdna", "cfDNA", "Tumor")) %>%
		 ggplot(aes(x=level, y=fixed_msi+.001, group=sample, fill=sample)) +
		 geom_bar(stat="identity", position=position_dodge(), , color="black", width=.8) +
		 facet_wrap(~subject_type, scales = "free_y") +
		 theme_classic(base_size=15) +
		 coord_flip() +
		 labs(y="\nMSI score\n", x="\n") +
		 scale_fill_manual(values = cols) +
		 guides(fill=guide_legend(title=c("Sample type"))) +
		 theme(legend.title=element_text(size=14)) +
		 ylim(0, 0.3) +
		 theme(
		 	axis.title.x = element_text(size = 16),
		 	axis.text.x = element_text(size = 16),
		 	strip.background = element_rect(colour = "white", fill = "white"),
		 	strip.text.x = element_text(size = 16))

pdf(file="../res/figureS7/msi_score_fixed.pdf", height=12, width=16)
h = matrix(NA, nrow=length(unique(msi_df %>% .[["patient_id"]])), ncol=3, dimnames=list(unique(msi_df %>% .[["patient_id"]]), c("Tumor", "cfDNA", "subj_type")))
h[,"Tumor"] = msi_df %>%
			  filter(sample=="tumor") %>%
			  .[["fixed_msi"]]
h[,"cfDNA"] = msi_df %>%
			  filter(sample=="cfdna") %>%
			  .[["fixed_msi"]]
subj_type = msi_df %>%
			filter(sample=="tumor") %>%
			.[["subject_type"]]
z = split.screen(matrix(c(0+.05+.025,1/3+.025,0,1, 1/3+.025+.025,2/3-.025+.025,0,1, 2/3+.025,1-.05+.025,0,1), nrow=3, ncol=4, byrow=TRUE))
screen(z[1])
x = barplot(t(subset(h[,c(2:1),drop=FALSE], subj_type=="Breast")), horiz=TRUE, las=1, beside=TRUE,
		col=rev(cols), xlim=c(0,.3), cex.axis = 1.25, lwd=1.85, space=c(0,.3))
axis(side=2, at=apply(x, 2, mean), labels=rep("", ncol(x)), line=0.35, tcl=-.25, lwd=1.35)
title(main="\n\nBreast", cex.main=1.15)
screen(z[2])
x = barplot(t(subset(h[,c(2:1),drop=FALSE], subj_type=="Lung")), horiz=TRUE, las=1, beside=TRUE,
		col=rev(cols), xlim=c(0,.3), cex.axis = 1.25, lwd=1.85, space=c(0,.3))
axis(side=2, at=apply(x, 2, mean), labels=rep("", ncol(x)), line=0.35, tcl=-.25, lwd=1.35)
title(main="\n\nLung", xlab="\n\n\n\nMSI score", cex.main=1.15, cex.lab=1.5)
screen(z[3])
x = barplot(t(subset(h[,c(2:1),drop=FALSE], subj_type=="Prostate")), horiz=TRUE, las=1, beside=TRUE,
		col=rev(cols), xlim=c(0,.3), cex.axis = 1.25, lwd=1.85, space=c(0,.3))
axis(side=2, at=apply(x, 2, mean), labels=rep("", ncol(x)), line=0.35, tcl=-.25, lwd=1.35)
title(main="\n\nProstate", cex.main=1.15)
close.screen(all.screens=TRUE)
dev.off()
