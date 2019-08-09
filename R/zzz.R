snv_file <- list(
	scored = str_c(
		"../",
		"modified_v11/",
		"Variants_Calls/",
		"Joined_cfDNA_IMPACT_variants/",
		"scored_merged_snvs_20171115.tsv"
	)
)

indel_file <- list(
	scored = str_c(
		"../",
		"modified_v11/",
		"Variants_Calls/",
		"Joined_cfDNA_IMPACT_variants/",
		"scored_merged_indels_20171115.tsv"
	)
)

hypermutators <- list(
	patient_id = c(
		"MSK-VB-0023",
		"MSK-VB-0044",
		"MSK-VB-0046",
		"MSK-VB-0050",
		"MSK-VB-0057",
		"MSK-VL-0035",
		"MSK-VL-0054",
		"MSK-VP-0031",
		"MSK-VP-0054"
	)
)

msi_hypermutators <- list(
	patient_id = c(
		"MSK-VP-0041"
	)
)

common_bed <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"Bed_files/",
	"pan_v2_wo_decoy_wo_iSNP_wo_CNV_IMPACT_common_regions.merged.bed"
)

common_bed_annotated <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"Bed_files/",
	"pan_v2_wo_decoy_wo_iSNP_wo_CNV_IMPACT_common_regions.merged_annotated.bed"
)

common_bed_annotated_w_introns <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"Bed_files/",
	"pan_v2_wo_decoy_wo_iSNP_wo_CNV_IMPACT_common_regions.merged_annotated_w_introns.bed"
)

patient_tracker <- str_c(
	"../",
	"modified_v11/",
	"Tracker_files/",
	"s3_sample_tracker_TechVal_Merlin.csv"
)

impact_tracker <- str_c(
	"../",
	"modified_v11/",
    "Tracker_files/",
    "sample_tracker_IMPACT_BAM_20170227.csv"
)

wbc_variants <- list(
	scored = str_c(
		"../",
		"modified_v11/",
		"Variants_Calls/",
		"Stacked_Scored_WBC/",
		"wbc_scored_annotated.tsv"
	),
	annotations = str_c(
		"../",
		"modified_v11/",
		"Resources/",
		"MSK_additional_data/",
		"wbc_scored_annotated.maf"
	)
)

clinical_file <- str_c(
	"../",
	"modified_v11/",
	"Variants_Calls/",
	"Joined_cfDNA_IMPACT_variants/",
	"clinical_no_dates_1_201703241549.tsv"
)

clinical_file_updated <- str_c(
	"../",
	"modified_v11/",
	"Variants_Calls/",
	"Joined_cfDNA_IMPACT_variants/",
	"clinical_no_dates_1_201703241549_updated.tsv"
)

somatic_vars_breast <- str_c(
	"../",
	"modified_v11/",
	"Variants_Calls/",
	"IMPACT_annotated_variants/",
	"normal_silent_common_intron_target.tables_merged_breast.tsv"
)

somatic_vars_lung <- str_c(
	"../",
	"modified_v11/",
	"Variants_Calls/",
	"IMPACT_annotated_variants/",
	"normal_silent_common_intron_target.tables_merged_lung.tsv"
)

somatic_vars_prostate <- str_c(
	"../",
	"modified_v11/",
	"Variants_Calls/",
	"IMPACT_annotated_variants/",
	"normal_silent_common_intron_target.tables_merged_prostate.tsv"
)

cosmic_v84 <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"External_database/",
	"cosmic_v84.tsv"
)

gnomad_r2.0.1 <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"External_database/",
	"gnomad_r2.0.1.tsv"
)

hotspot_v2 <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"External_database/",
	"cancer_hotspots_v2.tsv"
)

wbc_scored_annotated_and_clinical <- list(
	scored = str_c(
				"../",
				"modified_v11/",
				"Variants_Calls/",
				"Stacked_Scored_WBC/",
				"wbc_scored_annotated_280519.tsv"),
	clinical = str_c(
				"../",
				"modified_v11/",
				"Variants_Calls/",
				"Stacked_Scored_WBC/",
				"clinical_280519.tsv")
)

chip_genes <- c(
	"DNMT3A",
	"TET2",
	"ASXL1",
	"PPM1D",
	"TP53",
	"JAK2",
	"RUNX1",
	"SF3B1",
	"SRSF2",
	"IDH1",
	"IDH2",
	"U2AF1",
	"CBL",
	"ATM",
	"CHEK2"
)

somatic_snvs_grail <- list(
	scored = str_c(
    "../",
    "modified_v11/",
    "Variants_Calls/",
    "Joined_cfDNA_IMPACT_variants/",
    "scored_merged_snvs_20171115.tsv"
  )
)

somatic_indels_grail <- list(
  scored = str_c(
    "../",
    "modified_v11/",
    "Variants_Calls/",
    "Joined_cfDNA_IMPACT_variants/",
    "scored_merged_indels_20171115.tsv"
  )
)

msk_anno_joined <- str_c(
	"../",
	"modified_v11/",
	"Variants_Calls/",
	"Joined_cfDNA_IMPACT_variants/",
	"scripts2/",
	"20180412_MSK_TechVal_Grail_small_variants_bio_source_label.filter_vus_biopsy_only.all_cases_ccf.tsv"
)
						 
url_techval.repeats <- str_c(
	"../",
	"modified_v11/",
	"Variants_Calls/",
	"Joined_cfDNA_IMPACT_variants/",
	"annotated.tsv"
)

url_ddpcr <- str_c(
	"https://",
	"docs.google.com/",
	"spreadsheets/",
	"d/",
	"1s9PKGYZ1v3vgUiC25N6D5jXX1ai6c0vbisZsKIgwGlg/",
	"edit?usp=sharing"
)

url_retest <- str_c(
	"../",
	"modified_v11/",
	"Variants_Calls/",
	"Stacked_annotated_retest/",
	"TechVal_retest_annotated_stack.tsv"
)

url_cell.line <- str_c(
	"../",
	"modified_v11/",
	"Variants_Calls/",
	"Titrations_cell_line/",
	"scored/",
	"annotated.tsv"
)

url_HD753.truth <- str_c(
	"../",
	"modified_v11/",
	"Variants_Calls/",
	"Titrations_cell_line/",
	"hd753_manifest.tsv"
)

url_qc.metrics <- str_c(
	"../",
	"modified_v11/",
	"QC_metrics/",
	"TechVal_Merlin_QC_metrics.tsv"
)

url_ctdna_frac <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"MSK_additional_data/",
	"ctDNA_fraction.csv"
)

url_volumetric_data <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"MSK_additional_data/",
	"Volumetric_data.xls"
)

url_psa <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"MSK_additional_data/",
	"PSA_MSK_VP_0041.txt"
)

url_recist <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"MSK_additional_data/",
	"RECIST_MSK_VP_0041.txt"
)

url_msi_processed_data <- str_c(
	"../",
	"modified_v11/",
	"MSIMSK/",
	"data/",
	"20180423_msi_global_df.tsv"
)

url_subject_alias <- str_c(
	"../",
	"modified_v11/",
	"MSIMSK/",
	"data/",
	"MSK_TechVal_MS_patient_ID_alias_mapping.tsv"
)

url_prior_tx <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"MSK_additional_data/",
	"prior_tx_techval_0818.txt"
)

url_master_key <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"MSK_additional_data/",
	"master_sample_key.tsv"
)

url_cytogenetic_data <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"MSK_additional_data/",
	"hg19_cytoBandIdeo.txt"
)

url_leftover_samples <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"ddPCR/",
	"msk_techval_ddpcr_samples.tsv"
)

url_ddpcr_results <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"ddPCR/",
	"msk_techval_ddprc_results.csv"
)

url_smoking_history <- list(
	breast = str_c(
		"../",
		"modified_v11/",
		"Resources/",
		"MSK_additional_data/",
		"MSK_VB_SMOKING.csv"
	),
	prostate = str_c(
		"../",
		"modified_v11/",
		"Resources/",
		"MSK_additional_data/",
		"MSK_VP_SMOKING.csv"
	)
)

url_guardant_g360 <- str_c(
	"../",
	"modified_v11/",
	"Resources/",
	"MSK_additional_data/",
	"Guardant360_gene_list.tsv"
)

url_sample.tracker <- patient_tracker
url_msk.snv <- snv_file$scored
url_msk.indel <- indel_file$scored
germline_alpha <- 0.15
bam_read_count_cutoff <- 2
ad_cutoff <- 5
ref_genome <- BSgenome.Hsapiens.UCSC.hg19
url_original <- url_msk.snv
url_target.bed <- common_bed

cohort_cols <- c(
	"Control"	=	"#0063AA",
	"Breast"	=	"#EF145D",
	"Lung"		=	"#F98450",
	"Prostate"	=	"#009700"
)

variant_cols <- c(
	"VUSo"			=	"#4E9B3C",
	"biopsy_matched"	=	"#2C80C3",
	"IMPACT-BAM_matched"	=	"#F9CA03",
	"WBC_matched"		=	"#EA7180",
	"biopsy_only"		=	"#2C80C3"
)
