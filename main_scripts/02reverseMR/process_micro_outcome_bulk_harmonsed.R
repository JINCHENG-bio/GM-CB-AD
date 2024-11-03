process_micro_outcome_bulk_harmonsed <- function(i) {
  source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/get_bulk_harmonised_data.R")
  message("input = a table containg local microbiome GWAS summary data, names including:exposures	name	id	path")
  
  if(micro_exp$exposures[i] == "FR02") {
    sp <- fread(micro_exp$path[i]) %>% 
      mutate(bac = micro_exp$name[i], id = micro_exp$name[i]) 
    sp$p_value <- as.numeric(sp$p_value)
    sp <- sp %>% as.data.frame() %>% arrange(., p_value)
    
    sp_out_dat <- TwoSampleMR::format_data(
      dat = sp,
      type = "outcome",
      snps = NULL,
      header = TRUE,
      phenotype_col = "bac",
      snp_col = "hm_rsid",
      beta_col = "beta",
      se_col = "standard_error",
      eaf_col = "effect_allele_frequency",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "p_value",
      units_col = FALSE,
      ncase_col = FALSE,
      ncontrol_col = FALSE,
      samplesize_col = 5959,
      gene_col = FALSE,
      id_col = "id",
      min_pval = 1e-200,
      z_col = FALSE,
      info_col = FALSE,
      chr_col = "hm_chrom",
      pos_col = "hm_pos",
      log_pval = FALSE
    )
  } else if (micro_exp$exposures[i] == "Mibio") {
    sp <- fread(micro_exp$path[i]) %>% 
      mutate(id = micro_exp$name[i])
    sp$P.weightedSumZ <- as.numeric(sp$P.weightedSumZ)
    sp <- sp %>% as.data.frame() %>% arrange(., P.weightedSumZ) 
    
    sp_out_dat <- TwoSampleMR::format_data(
      dat = sp,
      type = "outcome",
      snps = NULL,
      header = TRUE,
      phenotype_col = "bac",
      snp_col = "rsID",
      beta_col = "beta",
      se_col = "SE",
      eaf_col = FALSE,
      effect_allele_col = "eff.allele",
      other_allele_col = "ref.allele",
      pval_col = "P.weightedSumZ",
      units_col = FALSE,
      ncase_col = FALSE,
      ncontrol_col = FALSE,
      samplesize_col = "N",
      gene_col = FALSE,
      id_col = "id",
      min_pval = 1e-200,
      z_col = "Z.weightedSumZ",
      info_col = FALSE,
      chr_col = "chr",
      pos_col = "bp",
      log_pval = FALSE
    )
  } else if (micro_exp$exposures[i] == "pathways") {
    sp <- fread(micro_exp$path[i]) %>% 
      mutate(pathway = micro_exp$exposure[i]) 
    sp$pval <- as.numeric(sp$pval)
    sp <- sp %>% as.data.frame() %>% arrange(., pval)
    
    sp_out_dat <- TwoSampleMR::format_data(
      dat = sp,
      type = "outcome",
      snps = NULL,
      header = TRUE,
      phenotype_col = "short",
      snp_col = "id",
      beta_col = "beta",
      se_col = "SE",
      eaf_col = "AF_Allele2",
      effect_allele_col = "alt",
      other_allele_col = "ref",
      pval_col = "pval",
      units_col = FALSE,
      ncase_col = FALSE,
      ncontrol_col = FALSE,
      samplesize_col = "num",
      gene_col = FALSE,
      id_col = "pathway",
      min_pval = 1e-200,
      z_col = FALSE,
      info_col = FALSE,
      chr_col = "chr",
      pos_col = "pos",
      log_pval = FALSE
    )
  }
  
  sp_out_dat$chr.outcome <- as.character(sp_out_dat$chr.outcome)
  
  tmp_AD_sp <- get_bulk_harmonised_data(exposure_dat = AD_clumped_exposures, outcome_dat = sp_out_dat)
  
  return(tmp_AD_sp)
}
