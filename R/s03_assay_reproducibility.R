#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
source('config.R')

if (!dir.exists("../res/figureS3")) {
	dir.create("../res/figureS3")
}

#==================================================
# Scatter plot of replicates for MSK-VB-0023
#==================================================
sample_tracker = read.csv(file=url_sample.tracker, header=TRUE, sep=",", stringsAsFactors=FALSE)
msk_snvs = read.csv(file=url_msk.snv, header=TRUE, sep="\t", stringsAsFactors=FALSE)
msk_indels = read.csv(file=url_msk.indel, header=TRUE, sep="\t", stringsAsFactors=FALSE)
techval_repeats = read.csv(file=url_techval.repeats, header=TRUE, sep="\t", stringsAsFactors=FALSE)

clean_target_region = read_tsv(url_target.bed, col_names = c("chrom", "start", "end"))
clean_target_region$in_target = TRUE

gdna_params = data_frame(
				subj_type = c("Healthy", "Breast", "Lung", "Prostate"),
				min_p = c(0.8, 0.79, 0.82, 0.79))

replicate2_id = techval_repeats %>%
				dplyr::select(patient_id) %>%
				distinct() %>%
				mutate(subject_id = substr(patient_id, 1, 11), replicate = "merlin")
replicate2 = left_join(replicate2_id, techval_repeats) %>%
			 filter(snv)
			 
replicate2_fixed = replicate2 %>%
				   mutate_at(vars(matches("qual")), funs(if_else(is.na(.), 0, .))) %>%
				   mutate(pedge = ifelse(is.na(pedge), 0, pedge), pgtkxgdna = ifelse(is.na(pgtkxgdna), 0, pgtkxgdna))

replicate2_fixed2 = replicate2_fixed %>%
					separate(hgvsp, into = c("p1", "p2", "p3", "p4"), sep = "\\|") %>%
					separate(is_nonsyn, into = c("t1","t2", "t3", "t4"), sep = "\\|") %>%
					separate(symbol, into = c("g1", "g2", "g3", "g4"), sep = "\\|") %>%
					mutate(gene = case_when(
									(.$t1 == TRUE | is.na(.$t2)) ~g1,
									(.$t1 != TRUE & .$t2 == TRUE) ~g2,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 == TRUE) ~g3,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 != TRUE & .$t4 == TRUE) ~g4)) %>%
					mutate(hgvs_p = case_when(
									(.$t1 == TRUE) ~p1,
									(.$t1 != TRUE & .$t2 == TRUE) ~p2,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 == TRUE) ~p3,
									(.$t1 != TRUE & .$t2 != TRUE & .$t3 != TRUE & .$t4 == TRUE) ~p4)) %>%
					mutate(hgvs_p = sub(".*:", "", hgvs_p),
						   is_nonsyn = ifelse(is.na(hgvs_p), FALSE, TRUE)) %>%
					dplyr::select(-t1, -t2, -t3, -t4, -g1, -g2, -g3, -g4, -p1, -p2, -p3, -p4)

replicate2_filtered = replicate2_fixed2 %>%
					  mutate(filter = "",
					  		 subj_type = case_when(
					  		 				grepl("VB", .$patient_id) ~ "Breast",
					  		 				grepl("VL", .$patient_id) ~ "Lung",
					  		 				grepl("VP", .$patient_id) ~ "Prostate")) %>%
					  left_join(gdna_params) %>%
					  mutate_(filter = ~update_filter(
					  		  filter = filter,
					  		  qualnobaq = qualnobaq,
					  		  pgtkxgdna = pgtkxgdna,
					  		  is_edge = isedge,
					  		  min_p = min_p)) %>%
					  dplyr::select(-min_p)

depth_lowest <- 200
replicate2_annotated = replicate2_filtered %>%
					   mutate(loc = str_c(chrom, ":", pos, "_", ref, ">", alt), af_point_cfdna = adnobaq/(dpnobaq+0.5), af_point_gdna = adgdna/(dpgdna+0.5)) %>%
					   mutate(gratio = (adnobaq+2)*(dpgdna+4)/((adgdna+2)*(dpnobaq+4)), gzero = adgdna/sqrt(adgdna+2)) %>%
					   mutate(category = case_when(
					   			.$isedge == TRUE ~ "edge",
					   			.$dpnobaq < depth_lowest ~ "low_depth",
					   			grepl("QUAL_LT", .$filter) ~ "low_qual",
					   			(.$af_point_gdna) > 0.3 ~ "germlineish",
					   			grepl("GDNA", .$filter) & .$adgdna > 0 ~ "blood",
					   			grepl("GDNA", .$filter) ~ "bloodish",
					   			.$gzero >= 1 & .$gratio < 6 ~ "bloodier",
					   			grepl("PASS", .$filter) ~ "somatic",
					   			TRUE ~ "other"))

replicate1_id = replicate2_id %>%
				mutate(patient_id = subject_id, replicate = "techval") %>%
				distinct()
				
replicate1 = left_join(replicate1_id, msk_snvs) %>%
			 filter(is_snv)

replicate1_annotated = replicate1 %>%
					   mutate(loc = str_c(chrom, ":", position, "_", ref, ">", alt), af_point_cfdna = adnobaq/(dpnobaq+0.5), af_point_gdna = adgdna/(dpgdna+0.5)) %>%
					   mutate(gratio = (adnobaq+2)*(dpgdna+4)/((adgdna+2)*(dpnobaq+4)),
					   		  gzero = adgdna/sqrt(adgdna+2)) %>%
					   mutate(category = case_when(
					   		  	.$isedge == TRUE ~ "edge",
					   		  	.$dpnobaq < depth_lowest ~ "low_depth",
					   		  	grepl("QUAL_LT", .$filter) ~ "low_qual",
					   		  	(.$af_point_gdna) > 0.3 ~ "germlineish",
					   		  	grepl("GDNA", .$filter) & .$adgdna > 0 ~ "blood",
					   		  	grepl("GDNA", .$filter) ~ "bloodish",
					   		  	.$gzero >= 1 & .$gratio < 6 ~ "bloodier",
					   		  	grepl("PASS", .$filter) ~ "somatic",
					   		  	TRUE ~ "other"))

replicate_all = full_join(replicate1_annotated, replicate2_annotated,
                          by = c("subject_id", "chrom", "position" = "pos", "ref", "alt"),
                          suffix = c(".1", ".2")) %>%
                          mutate(start = position,
                          		 end = position +1) %>%
                          genome_left_join(clean_target_region, by = c("chrom", "start", "end")) %>%
                          mutate(chrom = chrom.x) %>%
                          dplyr::select(-c(chrom.x, chrom.y, start.x, end.x, start.y, end.y))

replicate_all_filtered = replicate_all %>%
						 filter(in_target) %>%
						 mutate_at(vars(matches("^(afmean|af_point)")), funs(if_else(is.na(.), 0, .))) %>%
						 mutate_at(vars(matches("^(MSK)")), funs(if_else(is.na(.), 0L, .)))  %>%
						 mutate_at(vars(matches("^category")), funs(if_else(is.na(.), "x", .))) %>%
						 mutate(jointcat = (category.1 == "somatic" | category.2 == "somatic"),
						 jointad = pmax(adnobaq.1, adnobaq.2, na.rm=T)) %>%
						 mutate(cat_plus = case_when(
						 					.$category.1 == "x" ~ .$category.2,
						 					TRUE ~ .$category.1)) %>%
						 mutate(cat_minus = case_when(
						 					(.$category.1 == "x" | .$category.2 == "x") ~ "Not_detected",
						 					(.$category.1 == "somatic") ~ .$category.2,
						 					(.$category.2 == "somatic") ~ .$category.1,
						 					TRUE ~ .$category.1))

replicate_all_clean = replicate_all_filtered %>%
					  filter(!(category.1 == "edge" | category.2 == "edge"),
					  		 !(category.1 == "germlineish" | category.2 == "germlineish"),
					  		 !(category.1 == "low_depth" | category.2 == "low_depth"),
					  		 !(category.1 == "low_qual" | category.2 == "low_qual"),
					  		 !(cat_plus=="x"))

df_repro = replicate_all_filtered %>%
		   filter(in_target) %>%
		   filter(filter.1 == "PASS" | filter.2 == "PASS") %>%
		   filter(is_nonsyn.1 | is_nonsyn.2)
df_repro_hyper = df_repro %>%
				 filter(subject_id == "MSK-VB-0023")

df = df_repro_hyper %>%
     mutate(cat_minus = case_when(
       (.$cat_minus == "somatic") ~ "Called in both replicates",
       (.$cat_minus == "bloodier" | .$cat_minus == "bloodish" | .$cat_minus == "blood") ~ 
         "Incorrect assignment between replicates",
       (.$cat_minus == "edge" | .$cat_minus == "low_qual") ~ 
         "Not called in one replicate due to low quality",
       (.$cat_minus == "Not_detected") ~ "Not detected in one replicate")) %>%
     mutate(cat_minus = factor(cat_minus,
                              levels = c("Not detected in one replicate",
                                         "Not called in one replicate due to low quality",
                                         "Incorrect assignment between replicates",
                                         "Called in both replicates"),
                              ordered = TRUE)) %>%
     mutate(MSK = case_when(
       (.$MSK == "0") ~ "Biopsy-unmatched",
       (.$MSK == "1") ~ "Biopsy-matched")) %>%
     mutate(MSK = factor(MSK,
                         levels = c("Biopsy-unmatched", "Biopsy-matched"),
                         ordered = TRUE)) %>%
     mutate(af_point_cfdna.1 = af_point_cfdna.1*100,
            af_point_cfdna.2 = af_point_cfdna.2*100)
 
df1 = df %>% filter(MSK == "Biopsy-matched")
df2 = df %>% filter(MSK == "Biopsy-unmatched")
df3 = df %>% filter(!is.na(filter.1) & !is.na(filter.2))

repeats.fit = lm(af_point_cfdna.2~af_point_cfdna.1+0, data=df3)
fit.r2 = round(summary(repeats.fit)$r.squared, 4)
label.r2 = paste0("R^{2} == ", fit.r2)

pdf(file="../res/figureS3/assay_reproduce_combined_hyper.pdf", width=7, height=7)
par(mar = c(6.1, 6, 4.1, 1))
zz = split.screen(figs=matrix(c(0,1,0,1, 0.45,.95,0.075,.625) , nrow=2, ncol=4, byrow=TRUE))
screen(zz[1])
epsilon = 0
shapes = c("Biopsy-matched"=24,
		   "Biopsy-unmatched"=21)
cols = c("Not detected in one replicate"="#D7191C",
		 "Not called in one replicate due to low quality"="#FDAE61",
		 "Incorrect assignment between replicates"="#ABDDA4",
		 "Called in both replicates"="#2B83BA")
plot(1, 1, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, xlim=c(0, 1), ylim=c(0, 1))
axis(1, at = NULL, cex.axis = 1.5, padj = 0.25, lwd=1.25, lwd.ticks=1.35)
axis(2, at = NULL, cex.axis = 1.5, las = 1, lwd=1.5, lwd.ticks=1.35)
mtext(side = 1, text = "Replicate 1 (%)", line = 4, cex = 1.35)
mtext(side = 2, text = "Replicate 2 (%)", line = 4, cex = 1.35)
points(c(0,.99), c(0,.99)*repeats.fit$coefficients, type="l", lty=1, lwd=2, col="goldenrod3")
x = df2$af_point_cfdna.1
y = df2$af_point_cfdna.2
z1 = as.character(df2$MSK)
z2 = as.character(df2$cat_minus)
points(x, y, pch=shapes[z1], col="black", bg=cols[z2], cex=1.25)
x = df1$af_point_cfdna.1+epsilon
y = df1$af_point_cfdna.2+epsilon
z1 = as.character(df1$MSK)
z2 = as.character(df1$cat_minus)
z2[which(z1=="Biopsy-matched" & z2=="Incorrect assignment between replicates")] = "Called in both replicates"
points(x, y, pch=shapes[z1], col="black", bg=cols[z2], cex=1.25)
rect(xleft=0, ybottom=0.55, xright=0.4, ytop=1.0, col=transparent_rgb("white", 205), border=transparent_rgb("white", 205))
rect(xleft=0.475, ybottom=0, xright=1.0, ytop=0.475, col=transparent_rgb("white", 205), border=transparent_rgb("white", 205))
mtext(side=3, text="Variant category", line=-.85, at=.13, font=2)
legend(x=-.01, y=1.01, pch=21, col="black", pt.bg=cols[1], pt.cex=1.35, legend="Not detected in one replicate", box.lwd=-1, cex=.85)
legend(x=-.01, y=.97, pch=21, col="black", pt.bg=cols[2], pt.cex=1.35, legend="Not called in one replicate\ndue to low quality", box.lwd=-1, cex=.85)
legend(x=-.01, y=.89, pch=21, col="black", pt.bg=cols[3], pt.cex=1.35, legend="Incorrect assignment\nbetween replicates", box.lwd=-1, cex=.85)
legend(x=-.01, y=.80, pch=21, col="black", pt.bg=cols[4], pt.cex=1.35, legend="Called in both replicates", box.lwd=-1, cex=.85)
mtext(side=3, text="Biopsy concordance", line=-8.5, at=.165, font=2)
legend(x=-.01, y=.67, pch=21, col="black", pt.bg="black", pt.cex=1.25, legend="Biopsy unmatched", box.lwd=-1, cex=.85)
legend(x=-.01, y=.62, pch=24, col="black", pt.bg="black", pt.cex=1.25, legend="Biopsy matched", box.lwd=-1, cex=.85)
screen(zz[2])
epsilon = .05
plot(1, 1, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, xlim=c(0.05, 100), ylim=c(0.05, 100), log="xy")
rect(.045, .045, 1, 1, col="grey90", border="grey90")
axis(1, at = c(.05,1,10,50,100), labels=rep("",5), cex.axis = 1, padj = 0.15, lwd=1, lwd.ticks=1)
axis(1, at = c(.05,1,10,50,100), labels=c(0,1,10,50,100), cex.axis = 1, padj = 0.15, lwd=-1, line=-.5)
axis(1, at = 100, labels="    100", cex.axis = 1, padj = 0.15, lwd=-1, line=-.5)
axis(2, at = c(.05,1,10,50,100), labels=rep("",5), cex.axis = 1, las = 1, lwd=1, lwd.ticks=1, las=1)
axis(2, at = c(.05,1,10,50,100), labels=c(0,1,10,50,100), cex.axis = 1, lwd=-1, line=-.25, las=1)
points(c(0.05,100), c(0.05,100)*repeats.fit$coefficients, type="l", lty=1, lwd=1.25, col="goldenrod3")
x = df$af_point_cfdna.1+epsilon
y = df$af_point_cfdna.2+epsilon
z1 = as.character(df$MSK)
points(x, y,  pch=shapes[z1], col="salmon", cex=.55)
close.screen(all.screens=TRUE)
dev.off()