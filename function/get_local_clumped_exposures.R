get_local_clumped_exposures <- function(exposures_df = exposures_df, p_threhold = 1e-5)

{
  message("exposures_df = dataframe from FR02, mibio,pathways;containing: exposures	name	id	path")
  message("p_threhold = 1e-5, for gut microbiota")
  raw_exposures <- lapply(seq(1, dim(exposures_df)[1]), function(i){
  if(exposures_df$exposures[i] == "FR02")
  {
    sp <- fread(exposures_df$path[i]) %>% mutate(bac = exposures_df$name[i], id = exposures_df$name[i]) 
    sp$p_value <- as.numeric(sp$p_value)
    sp <- sp %>% as.data.frame()  %>% dplyr::filter(p_value < p_threhold) %>% arrange(., p_value)
        
    sp_exp_dat <- TwoSampleMR::format_data(
                                    dat = sp,
                                    type = "exposure",
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
    sp_exp_dat$chr.exposure <- as.character(sp_exp_dat$chr.exposure)
    
  }
  else if (exposures_df$exposures[i] == "Mibio")
  {
    sp <- fread(exposures_df$path[i]) %>% mutate(id = exposures_df$name[i]) 
    sp$P.weightedSumZ <- as.numeric(sp$P.weightedSumZ)
    sp <- sp %>% as.data.frame()  %>% dplyr::filter(P.weightedSumZ < p_threhold) %>% arrange(., P.weightedSumZ) 
    
    sp_exp_dat <- TwoSampleMR::format_data(
                                  dat = sp,
                                  type = "exposure",
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
    sp_exp_dat$chr.exposure <- as.character(sp_exp_dat$chr.exposure)
    }
    else if (exposures_df$exposures[i] == "pathways")
    {
       sp <- fread(exposures_df$path[i]) %>% mutate( pathway = exposures_df$name[i]) 
       sp$pval <- as.numeric(sp$pval)
       sp <- sp %>% as.data.frame() %>% dplyr::filter(pval < p_threhold) %>% arrange(., pval)
       
       sp_exp_dat <- TwoSampleMR::format_data(
                                      dat = sp,
                                      type = "exposure",
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
       sp_exp_dat$chr.exposure <- as.character(sp_exp_dat$chr.exposure)
      }
    return(sp_exp_dat)
  })
  
  # Remove NULL elements
  raw_exposures <- raw_exposures[!sapply(raw_exposures, is.null)]
    
  # Combine the list of data frames into a single data frame
  exposures <- list()
  exposures<- bind_rows(raw_exposures)
  source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/local_clumb.R")
  
  clumped_exposures <-local_clumb(dat= exposures, clump_kb = 10000, clump_r2 = 0.001, pop = "EUR") 
  
  return(clumped_exposures)
}
