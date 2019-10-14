#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/figureS15")) {
	dir.create("../res/figureS15")
}

if (!dir.exists("../res/etc/Source_Data_Extended_Data_Fig_7")) {
	dir.create("../res/etc/Source_Data_Extended_Data_Fig_7")
}

'plot_log3_' <- function(x, y, axis = TRUE, ylim=c(-2.5, 2.5))
{
   	par(mar=c(6.1, 9.5, 4.1, 1.1))
   	data(CytoBand)
   	end = NULL
	for (j in 1:23) {
		end = c(end, max(CytoBand$End[CytoBand$Chromosome==j]))
	}
	end = cumsum(end)
	start = rep(0, 23)
	start[2:23] = end[1:22]+1
	for (j in 1:23) {
		y[y[,"Chromosome"]==j,"Start"] = y[y[,"Chromosome"]==j,"Start"] + start[j]
		y[y[,"Chromosome"]==j,"End"] = y[y[,"Chromosome"]==j,"End"] + start[j]
		x[x[,"Chromosome"]==j,"Position"] = x[x[,"Chromosome"]==j,"Position"] + start[j]
	}
	plot(x[,"Position"], x[,"Log2Ratio"], type="p", pch=".", cex=1, col="grey75", axes=FALSE, frame=FALSE, xlab="", ylab="", main="", ylim=ylim)
	for (j in 1:nrow(y)) {
 		lines(x=c(y[j,"Start"], y[j,"End"]), y=rep(y[j,"Log2Ratio"],2), lty=1, lwd=2.5, col="red")
 	}
  	axis(2, at = c(-3, -2, -1, 0, 1, 2, 3), labels = c(-3, -2, -1, 0, 1, 2, 3), cex.axis = 1, las = 1)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3, cex = 1.15, las=1)
	points(c(-10000000, max(x[,"Position"])), c(0,0), type="l", col="black", lty=1, lwd=1)
	if (axis) {
		axis(side=1, at=c(start, end[length(end)]), labels=rep("", length(start)+1), tcl=.5)
		axis(1, at = .5*(start+end), labels=c(1:22, "X"), tcl=-.5, lwd=0, lwd.ticks=1, tcl=-.25)
	}

}

'undo_' <- function(x, n=10)
{
	cnm = matrix(NA, nrow=nrow(x), ncol=nrow(x))
	for (j in 1:nrow(x)) {
		cnm[,j] = abs(2^x[j,"Log2Ratio"] - 2^x[,"Log2Ratio"])
	}
	cnt = hclust(as.dist(cnm), "average")
	cnc = cutree(tree=cnt, k=n)
	for (j in unique(cnc)) {
		indx = which(cnc==j)
		if (length(indx)>2) {
 			mcl = mean(x[indx,"Log2Ratio"])
			scl = sd(x[indx,"Log2Ratio"])
			ind = which(x[indx,"Log2Ratio"]<(mcl+1.96*scl) & x[indx,"Log2Ratio"]>(mcl-1.96*scl))
			x[indx[ind],"Log2Ratio"] = mean(x[indx[ind],"Log2Ratio"])
		} else {
			x[indx,"Log2Ratio"] = mean(x[indx,"Log2Ratio"])
		}
	}
	return(x)
}

#==================================================
# 2-by-2 scatterplots of technical replicates
#==================================================
clinical = read_tsv(file=clinical_file, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))
		   
snv_vars = read_tsv(snv_file$scored, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   mutate(level_2a = as.character(level_2a)) %>%
		   mutate(level_r1 = as.character(level_r1))
		   
indel_vars = read_tsv(indel_file$scored, col_types = cols(.default = col_character())) %>%
			 type_convert() %>%
			 mutate(level_2a = as.character(level_2a)) %>%
		     mutate(level_r1 = as.character(level_r1))

wbc_stack = read_tsv(wbc_variants$scored, col_types = cols(.default = col_character())) %>%
			type_convert()
			
msk_anno = read_tsv(msk_anno_joined, col_types = cols(.default = col_character())) %>%
  		   type_convert()
			
tracker_grail = read_csv(file=patient_tracker)

tracker_impact = read_csv(file=impact_tracker)

valid_patient_ids = tracker_grail %>%
					filter(patient_id %in% tracker_impact$patient_id) %>%
					.[["patient_id"]]
  
indel_vars = indel_vars %>%
			 mutate(filter = replace(filter,
             		patient_id == "MSK-VB-0001" &
             		gene == "GATA3" &
             		filter == "PASS",
             		"CSR_MATCH_ELIMINATED"),
         	 		ccd = replace(ccd,
             			   		  patient_id == "MSK-VB-0001" &
                           		  gene == "GATA3" &
                           		  filter == "CSR_MATCH_ELIMINATED",
                           		  0))

snv_plasma = snv_vars %>%
  			 filter(ccd == 1,
         			(c_panel == 1 | panel == 1),
         			study == "TechVal",
         			grail == 1 | MSK == 1,
         			patient_id %in% valid_patient_ids) %>%
			 mutate(vtype = "SNV")

indel_plasma = indel_vars %>%
			   filter(ccd == 1,
         			  (c_panel == 1 | panel == 1),
         			  study == "TechVal",
         			  grail == 1 | MSK == 1,
         			  patient_id %in% valid_patient_ids) %>%
  			   mutate(vtype = "INDEL",
         			  altenddistmedian = as.integer(altenddistmedian))
         			  
healthy_snv = snv_vars %>%
  			  filter((c_panel == 1 | panel == 1),
         			  subj_type == "Healthy",
         			  grail == 1) %>%
  			  mutate(vtype = "SNV")

healthy_indel = indel_vars %>%
  				filter((c_panel == 1 | panel == 1),
         			    subj_type == "Healthy",
         				grail == 1) %>%
  				mutate(vtype = "INDEL",
         			   altenddistmedian = as.integer(altenddistmedian))

small_vars_plasma = full_join(snv_plasma, indel_plasma) %>%
					full_join(healthy_snv) %>%
					full_join(healthy_indel)
small_vars_plasma = small_vars_plasma %>%
  					mutate(subj_type = ifelse(subj_type == "Healthy", "Control", subj_type))

small_vars_plasma = small_vars_plasma %>%
					mutate(loc = str_c(chrom, ":", position_orig, "_", ref_orig, ">", alt_orig))

all_patient_table = small_vars_plasma %>%
					distinct(subj_type, patient_id)
					
all_patient_table = cbind.data.frame(subj_type = rep(all_patient_table$subj_type, 4),
                                     patient_id = rep(all_patient_table$patient_id, 4),
                                     bio_source = rep(c("WBC_matched",
                                                        "VUSo",
                                                        "biopsy_matched",
                                                        "biopsy_only"),
                                                       each = nrow(all_patient_table)))

variants = label_bio_source(small_vars_plasma)

variants = left_join(variants, msk_anno %>% dplyr::select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
variants = variants %>%
		   mutate(bio_source = case_when(
		   					   MSK == 1 & grail == 1 ~ "biopsy_matched",
		   					   MSK == 1 & grail == 0 ~ "biopsy_only",
		   					   category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
		   					   category %in% c("germline", "germlineish") ~ "germline",
		   					   category %in% c("blood", "bloodier") ~ "WBC_matched",
		   					   category == "somatic" & `IM-T.alt_count` > bam_read_count_cutoff ~ "IMPACT-BAM_matched",
		   					   category == "somatic" ~ "VUSo",
		   					   TRUE ~ "other"),
		   		  af = ifelse(is.na(af), 0, af),
		   		  af_nobaq = round(adnobaq / dpnobaq * 100, 2),
		   		  af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq))

patient_ids = c("MSK-VB-0023")
vars_rep0 = variants %>%
			filter(patient_id %in% patient_ids)

gdna_params = data_frame(
				subj_type = c("Healthy", "Breast", "Lung", "Prostate"),
				min_p = c(0.8, 0.79, 0.82, 0.79))
clean_target_region = read_tsv(url_target.bed, col_names = c("chrom", "start", "end"))
clean_target_region$in_target = TRUE
techval_repeats = read.csv(url_techval.repeats, header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
				  filter(!patient_id %in% c("MSK-VB-0023_B", "MSK-VL-0028_B", "MSK-VL-0042_B")) %>%
				  mutate(patient_id = gsub(pattern="_A", replacement="", patient_id)) %>%
				  mutate_at(vars(matches("qual")), funs(if_else(is.na(.), 0, .))) %>%
				  mutate(pedge = ifelse(is.na(pedge), 0, pedge), pgtkxgdna = ifelse(is.na(pgtkxgdna), 0, pgtkxgdna)) %>%
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
				  dplyr::select(-t1, -t2, -t3, -t4, -g1, -g2, -g3, -g4, -p1, -p2, -p3, -p4) %>%
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
				  dplyr::select(-min_p) %>%
				  mutate(loc = str_c(chrom, ":", pos, "_", ref, ">", alt)) %>%
				  mutate(position_orig = pos, ref_orig = ref, alt_orig = alt, position = pos) %>%
				  mutate(gratio = (adnobaq+2)*(dpgdna+4)/((adgdna+2)*(dpnobaq+4)),
				  		 gzero = adgdna/sqrt(adgdna+2)) %>%
				  mutate(start = position_orig,
                         end = position_orig + 1) %>%
                  genome_left_join(clean_target_region, by = c("chrom", "start", "end")) %>%
                  mutate(chrom = chrom.x) %>%
                  dplyr::select(-c(chrom.x, chrom.y, start.x, end.x, start.y, end.y)) %>%
				  left_join(variants %>% dplyr::select(patient_id, chrom, position_orig, ref_orig, alt_orig, MSK, grail), by=c("patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
				  mutate(MSK = ifelse(is.na(MSK), 0, MSK)) %>%
				  mutate(grail = 1) %>%
				  mutate(study = "TechVal") %>%
				  filter(in_target)
				  
feature_names = intersect(colnames(small_vars_plasma), colnames(techval_repeats))
repeat_variants = bind_rows(small_vars_plasma[,feature_names,drop=FALSE] %>%
							filter(!(patient_id %in% patient_ids)),
							techval_repeats[,feature_names,drop=FALSE])

repeat_variants = label_bio_source(repeat_variants)
repeat_variants = left_join(repeat_variants, msk_anno %>% dplyr::select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
repeat_variants = repeat_variants %>%
		   				mutate(bio_source = case_when(
		   					   MSK == 1 & grail == 1 ~ "biopsy_matched",
		   					   MSK == 1 & grail == 0 ~ "biopsy_only",
		   					   category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
		   					   category %in% c("germline", "germlineish") ~ "germline",
		   					   category %in% c("blood", "bloodier") ~ "WBC_matched",
		   					   category == "somatic" & `IM-T.alt_count` > bam_read_count_cutoff ~ "IMPACT-BAM_matched",
		   					   category == "somatic" ~ "VUSo",
		   					   TRUE ~ "other"),
			   		  	af_nobaq = round(adnobaq / dpnobaq * 100, 2),
			   		  	af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq))

vars_rep1 = repeat_variants %>%
			filter(patient_id %in% patient_ids)
			
techval_repeats = read.csv(url_techval.repeats, header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
				  filter(patient_id %in% c("MSK-VB-0023_B", "MSK-VL-0028_B", "MSK-VL-0042_B")) %>%
				  mutate(patient_id = gsub(pattern="_B", replacement="", patient_id)) %>%
				  mutate_at(vars(matches("qual")), funs(if_else(is.na(.), 0, .))) %>%
				  mutate(pedge = ifelse(is.na(pedge), 0, pedge), pgtkxgdna = ifelse(is.na(pgtkxgdna), 0, pgtkxgdna)) %>%
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
				  dplyr::select(-t1, -t2, -t3, -t4, -g1, -g2, -g3, -g4, -p1, -p2, -p3, -p4) %>%
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
				  dplyr::select(-min_p) %>%
				  mutate(loc = str_c(chrom, ":", pos, "_", ref, ">", alt)) %>%
				  mutate(position_orig = pos, ref_orig = ref, alt_orig = alt, position = pos) %>%
				  mutate(gratio = (adnobaq+2)*(dpgdna+4)/((adgdna+2)*(dpnobaq+4)),
				  		 gzero = adgdna/sqrt(adgdna+2)) %>%
				  mutate(start = position_orig,
                         end = position_orig + 1) %>%
                  genome_left_join(clean_target_region, by = c("chrom", "start", "end")) %>%
                  mutate(chrom = chrom.x) %>%
                  dplyr::select(-c(chrom.x, chrom.y, start.x, end.x, start.y, end.y)) %>%
				  left_join(variants %>% dplyr::select(patient_id, chrom, position_orig, ref_orig, alt_orig, MSK, grail), by=c("patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
				  mutate(MSK = ifelse(is.na(MSK), 0, MSK)) %>%
				  mutate(grail = 1) %>%
				  mutate(study = "TechVal") %>%
				  filter(in_target)
				  
feature_names = intersect(colnames(small_vars_plasma), colnames(techval_repeats))
repeat_variants = bind_rows(small_vars_plasma[,feature_names,drop=FALSE] %>%
							filter(!(patient_id %in% patient_ids)),
							techval_repeats[,feature_names,drop=FALSE])

repeat_variants = label_bio_source(repeat_variants)
repeat_variants = left_join(repeat_variants, msk_anno %>% dplyr::select(patient_id, chrom, position, ref, alt, CASE:complex_indel_duplicate))
repeat_variants = repeat_variants %>%
		   				mutate(bio_source = case_when(
		   					   MSK == 1 & grail == 1 ~ "biopsy_matched",
		   					   MSK == 1 & grail == 0 ~ "biopsy_only",
		   					   category %in% c("artifact", "edge", "low_depth", "low_qual") ~ "noise",
		   					   category %in% c("germline", "germlineish") ~ "germline",
		   					   category %in% c("blood", "bloodier") ~ "WBC_matched",
		   					   category == "somatic" & `IM-T.alt_count` > bam_read_count_cutoff ~ "IMPACT-BAM_matched",
		   					   category == "somatic" ~ "VUSo",
		   					   TRUE ~ "other"),
			   		  	af_nobaq = round(adnobaq / dpnobaq * 100, 2),
			   		  	af_nobaq = ifelse(is.na(af_nobaq), 0, af_nobaq))

vars_rep2 = repeat_variants %>%
			filter(patient_id %in% patient_ids)
			

all_vars = full_join(vars_rep0 %>% mutate(replicate = 1),
					 vars_rep1 %>% mutate(replicate = 2), by=c("study", "subj_type", "patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
		   dplyr::select(patient_id, adnobaq.x, adnobaq.y, dpnobaq.x, dpnobaq.y, bio_source.x, bio_source.y) %>%
		   mutate(af_rep1 = 100*adnobaq.x/dpnobaq.x) %>%
		   mutate(af_rep2 = 100*adnobaq.y/dpnobaq.y) %>%
		   mutate(bio_source.x = ifelse(is.na(bio_source.x), "unmatched", bio_source.x)) %>%
		   mutate(bio_source.y = ifelse(is.na(bio_source.y), "unmatched", bio_source.y)) %>%
		   filter(bio_source.x %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo") | bio_source.y %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
		   mutate(adnobaq.x = ifelse(is.na(adnobaq.x), 0, adnobaq.x)) %>%
		   mutate(adnobaq.y = ifelse(is.na(adnobaq.y), 0, adnobaq.y)) %>%
  		   mutate(dpnobaq.x = ifelse(is.na(dpnobaq.x), 0, dpnobaq.x)) %>%
   		   mutate(dpnobaq.y = ifelse(is.na(dpnobaq.y), 0, dpnobaq.y)) %>%
   		   mutate(afnobaq.x = 100*adnobaq.x/dpnobaq.x) %>%
   		   mutate(afnobaq.x = ifelse(afnobaq.x==0 | is.na(afnobaq.x), 0.01, afnobaq.x)) %>%
   		   mutate(afnobaq.y = 100*adnobaq.y/dpnobaq.y) %>%
   		   mutate(afnobaq.y = ifelse(afnobaq.y==0 | is.na(afnobaq.y), 0.01, afnobaq.y))


cols = c(
		 "Not detected in one replicate"="#D7191C",
		 "Not called in one replicate due\nto low quality"="#FDAE61",
		 "Incorrect assignment between replicates"="#ABDDA4",
		 "Called in both replicates"="#2B83BA"
		 )
		 
tmp_vars = all_vars %>% filter(patient_id == patient_ids)

# fix biopsy_matched in y
index = tmp_vars$bio_source.x == "biopsy_matched" & (tmp_vars$bio_source.y !="biopsy_matched" & tmp_vars$bio_source.y != "unmatched")
if (sum(index)!=0) {
	tmp_vars$bio_source.y[index] = "biopsy_matched"
}

# fix biopsy_matched in y
index = tmp_vars$bio_source.x == "IMPACT-BAM_matched" & (tmp_vars$bio_source.y !="IMPACT-BAM_matched" & tmp_vars$bio_source.y != "unmatched")
if (sum(index)!=0) {
	tmp_vars$bio_source.y[index] = "IMPACT-BAM_matched"
}

tmp_vars = tmp_vars %>%
		   mutate(shape = ifelse(bio_source.x=="biopsy_matched" | bio_source.x=="IMPACT-BAM_matched" | bio_source.y=="biopsy_matched" | bio_source.y=="IMPACT-BAM_matched", "Biopsy matched", "Biopsy unmatched")) %>%
		   mutate(fill = ifelse(bio_source.x!=bio_source.y , "Incorrect assignment between replicates", "Called in both replicates")) %>%
		   mutate(fill = ifelse((bio_source.x=="unmatched" & bio_source.y!="unmatched") | (bio_source.x!="unmatched" & bio_source.y=="unmatched"), "Not detected in one replicate", "Called in both replicates")) %>%
		   mutate(fill = ifelse((bio_source.x=="noise" & (bio_source.y!="noise" & bio_source.y!="other")) | (bio_source.y=="noise" & (bio_source.x!="noise" & bio_source.x!="other")), "Not called in one replicate due\nto low quality", fill)) %>%
		   filter(dpnobaq.x>10000) %>%
		   filter(bio_source.x!="noise")
		   
pdf(file=paste0("../res/figureS15/", patient_ids, "_R1_R2_High_Depth.pdf"), width=6.5, height=6.5)
par(mar = c(6.1, 6, 4.1, 1))
epsilon = 0
shapes = c("Biopsy matched"=24,
		   "Biopsy unmatched"=21)
cols = c("Not detected in one replicate"="#D7191C",
		 "Not called in one replicate due\nto low quality"="#FDAE61",
		 "Incorrect assignment between replicates"="#ABDDA4",
		 "Called in both replicates"="#2B83BA")
plot(1, 1, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, xlim=c(0.01, 100), ylim=c(0.01, 100), log="xy")
axis(1, at=c(0.01, 0.1, 1, 10, 100), labels=c("0", "0.1", "1", "10", "100"), cex.axis = 1.75, padj = 0.25, lwd=1.85, lwd.ticks=1.75)
axis(2, at=c(0.01, 0.1, 1, 10, 100), labels=c("0", "0.1", "1", "10", "100"), cex.axis = 1.75, las = 1, lwd=1.85, lwd.ticks=1.75)
mtext(side = 1, text = "Replicate 1 (%)", line = 4, cex = 1.75)
mtext(side = 2, text = "Replicate 2 (%)", line = 4, cex = 1.75)
points(c(.01,100), c(.01,100), type="l", lty=1, lwd=2, col="goldenrod3")
x = tmp_vars$afnobaq.x
y = tmp_vars$afnobaq.y
z1 = as.character(tmp_vars$shape)
z2 = as.character(tmp_vars$fill)
points(x, y, pch=shapes[z1], col="black", bg=cols[z2], cex=1.65)
log10_axis(side=1, at=c(0.01, 0.1, 1, 10, 100), lwd=0, lwd.ticks=1)
log10_axis(side=2, at=c(0.01, 0.1, 1, 10, 100), lwd=0, lwd.ticks=1)
dev.off()

export_x = tmp_vars %>%
		   mutate(tissue = case_when(grepl("VB", patient_id) ~ "Breast",
									 grepl("VL", patient_id) ~ "Lung",
			   				 		 grepl("VP", patient_id) ~ "Prostate")) %>%
		   dplyr::select(`patient_id` = `patient_id`,
			   			 `tissue` = `tissue`,
			   			 `af_rep_0` = `afnobaq.x`,
			   			 `af_rep_1` = `afnobaq.y`,
			   			 `variant_category` = `fill`,
			   			 `biopsy_concordance` = `shape`) %>%
			   			 mutate(variant_category = ifelse(variant_category == "Not called in one replicate due\nto low quality", "Not called in one replicate due to low quality", variant_category))
write_tsv(export_x, path="../res/etc/Source_Data_Extended_Data_Fig_7/Extended_Data_Fig_7g.tsv", append=FALSE, col_names=TRUE)

all_vars = full_join(vars_rep0 %>% mutate(replicate = 1),
					 vars_rep2 %>% mutate(replicate = 3), by=c("study", "subj_type", "patient_id", "chrom", "position_orig", "ref_orig", "alt_orig")) %>%
		   dplyr::select(patient_id, adnobaq.x, adnobaq.y, dpnobaq.x, dpnobaq.y, bio_source.x, bio_source.y) %>%
		   mutate(af_rep1 = 100*adnobaq.x/dpnobaq.x) %>%
		   mutate(af_rep2 = 100*adnobaq.y/dpnobaq.y) %>%
		   mutate(bio_source.x = ifelse(is.na(bio_source.x), "unmatched", bio_source.x)) %>%
		   mutate(bio_source.y = ifelse(is.na(bio_source.y), "unmatched", bio_source.y)) %>%
		   filter(bio_source.x %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo") | bio_source.y %in% c("biopsy_matched", "WBC_matched", "IMPACT-BAM_matched", "VUSo")) %>%
		   mutate(adnobaq.x = ifelse(is.na(adnobaq.x), 0, adnobaq.x)) %>%
		   mutate(adnobaq.y = ifelse(is.na(adnobaq.y), 0, adnobaq.y)) %>%
  		   mutate(dpnobaq.x = ifelse(is.na(dpnobaq.x), 0, dpnobaq.x)) %>%
   		   mutate(dpnobaq.y = ifelse(is.na(dpnobaq.y), 0, dpnobaq.y)) %>%
   		   mutate(afnobaq.x = 100*adnobaq.x/dpnobaq.x) %>%
   		   mutate(afnobaq.x = ifelse(afnobaq.x==0 | is.na(afnobaq.x), 0.01, afnobaq.x)) %>%
   		   mutate(afnobaq.y = 100*adnobaq.y/dpnobaq.y) %>%
   		   mutate(afnobaq.y = ifelse(afnobaq.y==0 | is.na(afnobaq.y), 0.01, afnobaq.y))


tmp_vars = all_vars %>% filter(patient_id == patient_ids)
	
# fix biopsy_matched in y
index = tmp_vars$bio_source.x == "biopsy_matched" & (tmp_vars$bio_source.y !="biopsy_matched" & tmp_vars$bio_source.y != "unmatched")
if (sum(index)!=0) {
	tmp_vars$bio_source.y[index] = "biopsy_matched"
}
	
# fix biopsy_matched in y
index = tmp_vars$bio_source.x == "IMPACT-BAM_matched" & (tmp_vars$bio_source.y !="IMPACT-BAM_matched" & tmp_vars$bio_source.y != "unmatched")
if (sum(index)!=0) {
	tmp_vars$bio_source.y[index] = "IMPACT-BAM_matched"
}
	
tmp_vars = tmp_vars %>%
		   mutate(shape = ifelse(bio_source.x=="biopsy_matched" | bio_source.x=="IMPACT-BAM_matched" | bio_source.y=="biopsy_matched" | bio_source.y=="IMPACT-BAM_matched", "Biopsy matched", "Biopsy unmatched")) %>%
		   mutate(fill = ifelse(bio_source.x!=bio_source.y , "Incorrect assignment between replicates", "Called in both replicates")) %>%
		   mutate(fill = ifelse((bio_source.x=="unmatched" & bio_source.y!="unmatched") | (bio_source.x!="unmatched" & bio_source.y=="unmatched"), "Not detected in one replicate", "Called in both replicates")) %>%
		   mutate(fill = ifelse((bio_source.x=="noise" & (bio_source.y!="noise" & bio_source.y!="other")) | (bio_source.y=="noise" & (bio_source.x!="noise" & bio_source.x!="other")), "Not called in one replicate due\nto low quality", fill)) %>%
		   filter(dpnobaq.x>10000) %>%
		   filter(bio_source.x!="noise")
			   	
pdf(file=paste0("../res/figureS15/", patient_ids, "_R1_R3_High_Depth.pdf"), width=6.5, height=6.5)
par(mar = c(6.1, 6, 4.1, 1))
epsilon = 0
shapes = c("Biopsy matched"=24,
		   "Biopsy unmatched"=21)
cols = c("Not detected in one replicate"="#D7191C",
		 "Not called in one replicate due\nto low quality"="#FDAE61",
		 "Incorrect assignment between replicates"="#ABDDA4",
		 "Called in both replicates"="#2B83BA")
plot(1, 1, type="n", xlab="", ylab="", main="", axes=FALSE, frame.plot=FALSE, xlim=c(0.01, 100), ylim=c(0.01, 100), log="xy")
axis(1, at=c(0.01, 0.1, 1, 10, 100), labels=c("0", "0.1", "1", "10", "100"), cex.axis = 1.75, padj = 0.25, lwd=1.85, lwd.ticks=1.75)
axis(2, at=c(0.01, 0.1, 1, 10, 100), labels=c("0", "0.1", "1", "10", "100"), cex.axis = 1.75, las = 1, lwd=1.85, lwd.ticks=1.75)
mtext(side = 1, text = "Replicate 1 (%)", line = 4, cex = 1.75)
mtext(side = 2, text = "Replicate 3 (%)", line = 4, cex = 1.75)
points(c(.01,100), c(.01,100), type="l", lty=1, lwd=2, col="goldenrod3")
x = tmp_vars$afnobaq.x
y = tmp_vars$afnobaq.y
z1 = as.character(tmp_vars$shape)
z2 = as.character(tmp_vars$fill)
points(x, y, pch=shapes[z1], col="black", bg=cols[z2], cex=1.65)
log10_axis(side=1, at=c(0.01, 0.1, 1, 10, 100), lwd=0, lwd.ticks=1)
log10_axis(side=2, at=c(0.01, 0.1, 1, 10, 100), lwd=0, lwd.ticks=1)
dev.off()

export_x = tmp_vars %>%
		   mutate(tissue = case_when(grepl("VB", patient_id) ~ "Breast",
									 grepl("VL", patient_id) ~ "Lung",
			   				 		 grepl("VP", patient_id) ~ "Prostate")) %>%
		   dplyr::select(`patient_id` = `patient_id`,
			   			 `tissue` = `tissue`,
			   			 `af_rep_0` = `afnobaq.x`,
			   			 `af_rep_1` = `afnobaq.y`,
			   			 `variant_category` = `fill`,
			   			 `biopsy_concordance` = `shape`) %>%
			   			 mutate(variant_category = ifelse(variant_category == "Not called in one replicate due\nto low quality", "Not called in one replicate due to low quality", variant_category))
write_tsv(export_x, path="../res/etc/Source_Data_Extended_Data_Fig_7/Extended_Data_Fig_7h.tsv", append=FALSE, col_names=TRUE)

#==================================================
# log2 ratio plots grail cfdna tumor samples
#==================================================
load("../res/rebuttal/uncollapsed_bam/cnvkit/totalcopy/MSK-VB-0023-T.RData")
tmp2 = winsorize(CN, method="mad", tau=3.5, verbose=FALSE)
colnames(tmp2) = c("Chromosome","Position","Log2Ratio")
tmp = undo_(tmp, n=2)
pdf(file="../res/figureS15/MSK-VB-0023_cfDNA.pdf", width=8, height=4)
plot_log3_(x=tmp2, y=tmp, axis=TRUE, ylim=c(-3.25,3.25))
dev.off()

export_x = tmp2 %>%
		   dplyr::select(`chromosome` = `Chromosome`,
		   				 `position` = `Position`,
		   				 `log2ratio` = `Log2Ratio`)
export_y = tmp %>%
		   dplyr::select(`chromosome` = `Chromosome`,
		   				 `arm` = `Arm`,
		   				 `start` = `Start`,
		   				 `end` = `End`,
		   				 `n` = `N`,
		   				 `log2ratio` = `Log2Ratio`)
write_tsv(export_x, path="../res/etc/Source_Data_Extended_Data_Fig_7/Extended_Data_Fig_7e_1.tsv", append=FALSE, col_names=TRUE)
write_tsv(export_y, path="../res/etc/Source_Data_Extended_Data_Fig_7/Extended_Data_Fig_7e_2.tsv", append=FALSE, col_names=TRUE)


#==================================================
# log2 ratio plots msk-impact tumor samples
#==================================================
key_file = read_tsv(file=url_master_key, col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   dplyr::select(PATIENT_ID, GRAIL_ID, DMP_ID, TUMOR_ID, NORMAL_ID, GRAIL_alpha, GRAIL_psi, IMPACT_alpha, IMPACT_psi) %>%
		   filter(GRAIL_ID == "MSK-VB-0023")

load(paste0("../res/rebuttal/msk_impact/cnvkit/totalcopy/", key_file$TUMOR_ID, ".RData"))
tmp2 = winsorize(CN, method="mad", tau=4.5, verbose=FALSE)
colnames(tmp2) = c("Chromosome","Position","Log2Ratio")
tmp = pcf(data=tmp2, kmin=10, gamma=40, verbose=FALSE)[,2:7,drop=FALSE]
colnames(tmp) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
tmp = undo_(tmp, n=3)
pdf(file=paste0("../res/figureS15/", key_file$GRAIL_ID, "_Tumor.pdf"), width=8, height=4)
plot_log3_(x=tmp2, y=tmp, axis=FALSE, ylim=c(-3.25,3.25))
dev.off()

export_x = tmp2 %>%
		   dplyr::select(`chromosome` = `Chromosome`,
		   				 `position` = `Position`,
		   				 `log2ratio` = `Log2Ratio`)
export_y = tmp %>%
		   dplyr::select(`chromosome` = `Chromosome`,
		   				 `arm` = `Arm`,
		   				 `start` = `Start`,
		   				 `end` = `End`,
		   				 `n` = `N`,
		   				 `log2ratio` = `Log2Ratio`)
write_tsv(export_x, path="../res/etc/Source_Data_Extended_Data_Fig_7/Extended_Data_Fig_7f_1.tsv", append=FALSE, col_names=TRUE)
write_tsv(export_y, path="../res/etc/Source_Data_Extended_Data_Fig_7/Extended_Data_Fig_7f_2.tsv", append=FALSE, col_names=TRUE)
