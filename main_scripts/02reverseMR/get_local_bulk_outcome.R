get_local_outcomes <- function(outcomes_df = outcomes_df)

{
  library(data.table)
  message("outcomes_df = dataframe from FR02, mibio,pathways;containing: exposures	name	id	path")
  raw_outcomes <- lapply(seq(1, dim(outcomes_df)[1]), function(i){
  if(outcomes_df$exposures[i] == "FR02")
  {
    sp <- fread(outcomes_df$path[i]) %>% mutate(bac = outcomes_df$name[i], id = outcomes_df$name[i]) 
    sp$p_value <- as.numeric(sp$p_value)
    sp <- sp %>% as.data.frame() %>% arrange(., p_value)
        
    sp_exp_dat <- TwoSampleMR::format_data(
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
    sp_exp_dat$chr.outcome <- as.character(sp_exp_dat$chr.outcome)
    
  }
  else if (outcomes_df$exposures[i] == "Mibio")
  {
    sp <- fread(outcomes_df$path[i]) %>% mutate(id = outcomes_df$name[i]) 
    sp$P.weightedSumZ <- as.numeric(sp$P.weightedSumZ)
    sp <- sp %>% as.data.frame()  %>% arrange(., P.weightedSumZ) 
    
    sp_exp_dat <- TwoSampleMR::format_data(
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
    sp_exp_dat$chr.outcome <- as.character(sp_exp_dat$chr.outcome)
    }
    else if (outcomes_df$exposures[i] == "pathways")
    {
       sp <- fread(outcomes_df$path[i]) %>% mutate( pathway = outcomes_df$name[i]) 
       sp$pval <- as.numeric(sp$pval)
       sp <- sp %>% as.data.frame() %>% arrange(., pval)
       
       sp_exp_dat <- TwoSampleMR::format_data(
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
       sp_exp_dat$chr.outcome <- as.character(sp_exp_dat$chr.outcome)
      }
    return(sp_exp_dat)
  })
  
  # Remove NULL elements
  raw_outcomes <- raw_outcomes[!sapply(raw_outcomes, is.null)]
    
  # Combine the list of data frames into a single data frame
  outcomes <- list()
  outcomes <- bind_rows(raw_outcomes)
  
  return(outcomes)
}
