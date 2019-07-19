#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source('config.R')

if (!dir.exists("../res/rebuttal")) {
	dir.create("../res/rebuttal")
}

tracker = read_csv(file=patient_tracker, col_types = cols(.default = col_character()))  %>%
 				   type_convert()
 				   
qc = read_tsv(file=url_qc.metrics, col_types = cols(.default = col_character()))  %>%
 			  type_convert()
 			  
snv = read_tsv(snv_file$scored, col_types = cols(.default = col_character())) %>%
		       type_convert() %>%
		       mutate(level_2a = as.character(level_2a)) %>%
		       mutate(level_r1 = as.character(level_r1))

ver_tidy = "../modified_v11/"
s3_target = paste0(ver_tidy, "Resources/Bed_files/pan_v2_wo_decoy_wo_iSNP_wo_CNV_IMPACT_common_regions.merged.bed")		       
bed = read_tsv(file=s3_target, col_names = c("chrom", "start", "end"), col_types = cols(.default = col_character())) %>%
		       type_convert() %>%
			   mutate(is_clean = TRUE)

vuso = read_tsv(file=paste0(ver_tidy, "Resources/MSK_additional_data/Table_S7_vuso_biorad_ddpcr_annot.tsv"), col_types = cols(.default = col_character())) %>%
		        type_convert()

leftover = read_tsv(file=url_leftover_samples, col_types = cols(.default = col_character())) %>%
		            type_convert()

results = read_csv(file=url_ddpcr_results, col_types = cols(.default = col_character())) %>%
		           type_convert()
		           
vuso_target = vuso %>%
  			  filter(!is.na(FAM_Assay_ID)) %>%
  			  mutate(chrom = paste0("chr", Chromosome)) %>%
  			  dplyr::select(patient_id = Tumor_Sample_Barcode,
  			  		 gene = Hugo_Symbol,
  			  		 chrom,
  			  		 position = Start_Position,
  			  		 ref = Reference_Allele,
  			  		 alt = Tumor_Seq_Allele2,
  			  		 HGVSc,
  			  		 HGVSp_Short,
  			  		 UUID,
  			  		 FAM_Assay_ID,
  			  		 HEX_Assay_ID,
  			  		 FAM_HEX_Assay_ID)
  			  		 
id_candidate = as.character(unique(vuso_target$patient_id))

ddpcr_target = vuso_target %>%
  			   dplyr::select(-patient_id) %>%
  			   unique()
  			   
lib_op = data.frame(Library_preparation_performed_by = c("OP1", "OP2", "OP3", "OP4", "OP5"),
                    Library_prep_operator = c("CT", "JN", "YP", "TL", "MG"))
                    
qc_short = qc %>%
  		   filter(sample_type == "cfDNA") %>%
  		   filter(Study == "TechVal" & tissue != "Healthy") %>%
  		   left_join(lib_op) %>%
  		   mutate(Library_yield_ng = ifelse(is.na(Library_yield_ng), Adapter_dimer_ng*(100-Percent_dimer)/Percent_dimer, Library_yield_ng)) %>%
  		   dplyr::select(patient_id,
  		   		  sample_id,
  		   		  tube_id,
  		   		  tissue,
                  DNA_extraction_date,
                  DNA_extraction_yield_ng,
                  Library_preparation_start_date,
                  Library_prep_operator,
                  Library_yield_ng)
                  
leftover_samples = leftover %>%
				   filter(cfDNA_evaluable) %>%
				   dplyr::select(patient_id,
				   		  tube_id,
				   		  ddPCR_candidate_sample,
				   		  striptube_label,
				   		  striptube_position,
				   		  leftover_cfDNA_volume_ul,
				   		  leftover_undiluted_library_volume_ul,
				   		  leftover_10xdiluted_library_volume_ul) %>%
				   left_join(qc_short) %>%
				   mutate(leftover_cfDNA_ng = DNA_extraction_yield_ng/60*leftover_cfDNA_volume_ul,
         				  leftover_undiluted_library_ng = Library_yield_ng/20*leftover_undiluted_library_volume_ul,
         				  leftover_10xdiluted_library_ng = Library_yield_ng/20/10*leftover_10xdiluted_library_volume_ul,
         				  leftover_library_ng = pmax(leftover_undiluted_library_ng, leftover_10xdiluted_library_ng)) %>%
  				   mutate(cfDNA_available = ifelse(leftover_cfDNA_ng > 100, TRUE, FALSE), library_available = ifelse(leftover_library_ng > 300, TRUE, FALSE)) %>%
  				   filter(ddPCR_candidate_sample | library_available)

id_leftover = as.character(unique(leftover_samples$patient_id))

snv_short = snv %>%
			mutate(start = position, end = position +1) %>%
			genome_left_join(bed, by = c("chrom", "start", "end")) %>%
			mutate(chrom = chrom.x) %>%
			dplyr::select(-c(chrom.x, chrom.y, start.x, end.x, start.y, end.y)) %>%
			filter(is_clean) %>%
			filter(patient_id %in% id_leftover) %>%
			mutate(ref = ifelse(is.na(ref), c_ref, ref),
				   alt = ifelse(is.na(alt), c_alt, alt),
				   gene = ifelse(is.na(gene), c_gene, gene),
				   hgvsc = ifelse(is.na(hgvs_c), cdna_change, hgvs_c),
				   hgvsp = ifelse(is.na(hgvs_p), cdna_change, hgvs_p),
         		   hgvsp_short = ifelse(is.na(hgvsp_short), amino_acid_change, hgvsp_short)) %>%
        	mutate(adnobaq = ifelse(is.na(adnobaq), 0, adnobaq), af_cfdna = ifelse(is.na(dpnobaq), 0, adnobaq/dpnobaq)) %>%
        	dplyr::select(patient_id,
        		   chrom,
        		   position,
        		   ref,
        		   alt,
                   gene,
                   hgvsc,
                   hgvsp,
                   hgvsp_short, 
                   MSK,
                   grail,
                   t_depth,
                   t_ref_count,
                   t_alt_count,
                   adnobaq,
                   dpnobaq,
                   af_cfdna,
                   filter) %>%
  			full_join(ddpcr_target)

df_merged = full_join(leftover_samples, snv_short)

df_merged %>% filter(ddPCR_candidate_sample) %>%
  			filter(MSK == 0) %>%
  			filter(!is.na(FAM_HEX_Assay_ID)) %>%
  			arrange(patient_id, -af_cfdna) %>%
  			mutate(`cfDNA MAF (%)` = round(af_cfdna*100, 2)) %>%
  			dplyr::select(patient_id, gene, HGVSp_Short, `cfDNA MAF (%)`, `ddPCR assay` = FAM_HEX_Assay_ID) %>%
  			pandoc.table(caption = "Candidate VUSo for ddPCR Retest", split.table = Inf)
  			
df_leftover = df_merged %>%
			  filter(ddPCR_candidate_sample) %>%
			  filter(MSK == 0) %>%
			  filter(!is.na(FAM_HEX_Assay_ID))

agg_probes = aggregate(data = df_leftover, UUID ~ patient_id, toString) %>%
			 mutate(ddPCR_assay = UUID) %>%
			 dplyr::select(-UUID)

df_leftover %>% arrange(-library_available, -cfDNA_available, patient_id) %>%
				dplyr::select(patient_id, `cfDNA (ng)` = leftover_cfDNA_ng, `library (ng)` = leftover_library_ng) %>%
				unique() %>%
				full_join(agg_probes) %>%
				pandoc.table(caption = "Estimated Leftover cfDNA and Library Amount", split.table = Inf, round = 0)

retest = df_merged %>%
		 filter(ddPCR_candidate_sample) %>%
		 filter(filter == "PASS" | MSK == 1) %>%
		 filter(library_available | cfDNA_available) %>%
		 filter(af_cfdna < 0.1) %>%
		 mutate(variant_name = paste0(gene, "_", hgvsp)) %>%
		 group_by(patient_id) %>%
		 mutate(select_ddPCR_assay = ifelse(af_cfdna == max(af_cfdna[!is.na(FAM_Assay_ID)]), TRUE, FALSE),
         		max_af_cfdna = max(af_cfdna),
         		max_af_matched = ifelse(max(af_cfdna[MSK == 1]) < 0, NA, max(af_cfdna[MSK == 1])),
         		variant_max_af_cfdna = variant_name[af_cfdna == max_af_cfdna],
         		variant_max_af_matched = ifelse(is.na(max_af_matched), NA, variant_name[af_cfdna == max_af_matched])) %>%
         ungroup()

retest %>% filter(select_ddPCR_assay) %>%
		   arrange(patient_id) %>%
		   mutate(af_cfdna = round(af_cfdna*100,2),
         		  max_af_cfdna = round(max_af_cfdna*100,2),
         		  max_af_matched = round(max_af_matched*100,2)) %>%
           dplyr::select(patient_id, `ddPCR assay` = UUID, `MAF (%)` = af_cfdna, `Highest MAF SNV` = variant_max_af_cfdna, `Highest MAF (%)` = max_af_cfdna) %>%
           pandoc.table(split.table = Inf, caption = "Selected Paient Samples for ddPCR Retest")

n_cfDNA <- length(unique(retest$patient_id[retest$cfDNA_available]))
n_library <- length(unique(retest$patient_id[retest$library_available]))
n_ddPCR_assay <- length(unique(retest$UUID[retest$select_ddPCR_assay]))
n_positive_runs <- 3*n_cfDNA + 3*n_library
n_negative_runs <- 3*n_ddPCR_assay

select_ddPCR = retest %>%
			   filter(select_ddPCR_assay) %>%
			   dplyr::select(UUID, chrom, position, ref, alt) %>%
			   unique()

id_blacklisted = df_merged %>%
				 right_join(select_ddPCR) %>%
				 dplyr::select(patient_id) %>%
				 unique()

id_neg = setdiff(id_leftover, id_blacklisted$patient_id)

neg = df_merged %>%
	  filter(library_available) %>%
	  filter(patient_id %in% id_neg) %>%
	  filter(filter == "PASS") %>%
	  mutate(variant_name = paste0(gene, "_", hgvsp)) %>%
	  group_by(patient_id) %>%
	  mutate(max_af_cfdna = max(af_cfdna),
	  		 variant_max_af_cfdna = variant_name[af_cfdna == max_af_cfdna],
	  		 n_variants = n()) %>%
	  mutate(max_af_cfdna = round(max_af_cfdna*100,2)) %>%
	  ungroup()

neg %>% arrange(tissue, max_af_cfdna) %>%
		dplyr::select(patient_id, `library (ng)` = leftover_library_ng, `Highest cfDNA MAF (%)` = max_af_cfdna, `N cfDNA variants` = n_variants,  `ddPCR target` = UUID) %>%
  		unique() %>%
  		pandoc.table(split.table = Inf, caption = "Candidate Negative Patient Samples")

ddpcr_results = results %>%
				filter(!is.na(patient_id)) %>%
				filter(Target != "WT") %>%
				mutate(patient_id = gsub("_", "-", patient_id)) %>%
				mutate(af_ddPCR = ifelse(is.na(`Fractional Abundance`), 0, `Fractional Abundance`/100)) %>%
				dplyr::select(patient_id, input_type, rep, Sample, UUID = Target, af_ddPCR)

ddpcr_results2 = df_merged %>%
				 filter(filter == "PASS" | MSK == 1) %>%
				 full_join(ddpcr_results, by = c("patient_id", "UUID")) %>%
				 filter(!is.na(af_ddPCR)) %>%
				 mutate(af_cfdna = ifelse(is.na(af_cfdna), 0, af_cfdna)) %>%
				 group_by(patient_id) %>%
				 mutate(ddPCR_retest = ifelse(n() == 1, "Negative", "Positive")) %>%
				 ungroup()

ddpcr_results_pos = ddpcr_results2 %>%
					filter(ddPCR_retest == "Positive") %>%
					group_by(patient_id, UUID, ddPCR_retest, input_type) %>%
					summarise(af_ddPCR_mean = mean(af_ddPCR), af_ddPCR_sd = sd(af_ddPCR), af_cfdna = unique(af_cfdna))
					
ddpcr_results_neg = ddpcr_results2 %>%
					filter(ddPCR_retest == "Negative") %>%
					group_by(UUID, ddPCR_retest, input_type) %>%
					summarise(patient_id = NA, af_ddPCR_mean = mean(af_ddPCR), af_ddPCR_sd = sd(af_ddPCR), af_cfdna = unique(af_cfdna))

ddpcr_results_summary = rbind(ddpcr_results_pos, ddpcr_results_neg)

ddpcr_results_summary %>% arrange(UUID, ddPCR_retest, patient_id) %>%
						  dplyr::select(`ddPCR assay` = UUID,
						  		 ddPCR_retest,
						  		 patient_id,
						  		 input_type, 
						  		 `AF ddPCR, mean` = af_ddPCR_mean,
						  		 `AF ddPCR, std` = af_ddPCR_sd,
                				 `AF GRAIL` = af_cfdna) %>%
  						  pandoc.table(split.table = Inf, caption = "ddPCR Tested Samples")
  						  
tmp = ddpcr_results_summary %>%
	  dplyr::select(patient_id, UUID, ddPCR_retest, input_type, af_ddPCR_mean, af_ddPCR_sd, af_cfdna) %>%
	  mutate(af_ddPCR_mean = ifelse(af_ddPCR_mean < 5e-5, 2e-4, af_ddPCR_mean)) %>%
	  mutate(af_cfdna = ifelse(af_cfdna == 0, 2e-4, af_cfdna)) %>%
	  mutate(af_ddPCR_mean = af_ddPCR_mean*100) %>%
	  mutate(af_cfdna = af_cfdna*100) %>%
	  mutate(facet = "VUSo") %>%
	  mutate(input_type = factor(ifelse(input_type=="library", "Library", input_type))) %>%
	  mutate(UID = gsub("_p.", " ", as.character(UUID), fixed=TRUE))

cols = c("cfDNA"="#D7191C",
		 "Library"="#2B83BA")
		 
plot.0 = ggplot(tmp, aes(y = af_cfdna, x = af_ddPCR_mean, shape = ddPCR_retest, fill = input_type, label = UID)) +
		 geom_abline(slope = 1, color = "goldenrod3", linetype = 1) +
		 geom_point(alpha = .8, size = 2.5) +
		 scale_fill_manual(values = cols) +
		 scale_shape_manual(values = c(24, 21)) +
		 geom_text_repel(size = 2.25, color = "black", segment.color = "transparent") +
		 theme_bw(base_size=15) +
		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.text=element_text(size=9), legend.title=element_text(size=10), legend.position = c(0.2, 0.75), legend.background = element_blank(), legend.key.size = unit(1, 'lines')) +
		 labs(y="\nTargeted DNA assay (%)\n", x="\nddPCR (%)\n") +
		 scale_x_log10(
		 	breaks = function(x) { c(0.02, 0.05, 0.1, 1, 10) },
 			labels = function(x) { c("ND", ".05", ".1", "1", "10") }
		 ) +
		 scale_y_log10(
		 	breaks = function(x) { c(0.02, 0.05, 0.1, 1, 10) },
 			labels = function(x) { c("ND", ".05", ".1", "1", "10") }
		 ) +
		 annotation_logticks() +
		 coord_cartesian(xlim=c(2e-2,10), ylim = c(2e-2, 10)) +
		 facet_wrap(~facet) +
		 guides(shape=guide_legend(title=c("Called in targeted\ncfDNA assay"), override.aes=list(fill="black"))) +
		 guides(fill=guide_legend(title=c("Input type")))
		 
pdf(file="../res/rebuttal/ddPCR_VUSo.pdf", width=6, height=6)
print(plot.0)
dev.off()


fill = c("cfDNA"="#D7191C", "Library"="#2B83BA")
shape = c("Positive"=24, "Negative"=21)

pdf(file="../res/rebuttal/ddPCR_VUSo.pdf", width=8, height=8)
par(mar = c(6.1, 6, 4.1, 1))
plot(1,1, type="n", xlim=c(2e-2,10), ylim=c(2e-2,10), xlab="", ylab="", axes=FALSE, frame.plot=FALSE, log="xy")
points(x=c(0.02, 9), y=c(0.02,9), type="l", lty=1, lwd=2, col="goldenrod3")
x = tmp$af_ddPCR_mean
y = tmp$af_cfdna
points(x, y, type="p", col="black", bg=fill[tmp$input_type], pch=shape[tmp$ddPCR_retest], lwd=.5, cex=1.55)
axis(1, at = c(.02, .05, .10, .25, 1.0, 2.5, 5.0, 10.0), labels = c("ND", ".05", ".1", ".25", "1", "2.5", "5", "10"), cex.axis = 1.5, padj = 0.25, lwd=1.5, lwd.ticks=1.35)
axis(2, at = c(.02, .05, .10, .25, 1.0, 2.5, 5.0, 10.0), labels = c("ND", ".05", ".1", ".25", "1", "2.5", "5", "10"), cex.axis = 1.5, las = 1, lwd=1.5, lwd.ticks=1.35)
mtext(side = 1, text = "ddPCR (%)", line = 4, cex = 1.85)
mtext(side = 2, text = "Targeted DNA assay (%)", line = 4, cex = 1.85)
text(x=0.03, y=10.5, labels="Called in targeted\ncfDNA assay")
legend(x=0.02, y=10, pch=shape, col="black", pt.bg="black", pt.cex=1.55, pt.lwd=.5, legend=paste0(" ", names(shape)), box.lwd=-1, cex=1.15)
text(x=0.03, y=4.5, labels="Input type")
legend(x=0.02, y=3.5, pch=21, col="black", pt.bg=fill, pt.cex=1.55, pt.lwd=.5, legend=paste0(" ", names(fill)), box.lwd=-1, cex=1.15)
dev.off()


