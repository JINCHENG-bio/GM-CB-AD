get_local_exposures_forMRBMA<- function(inputdf=inputdf,
                                        exposure=""
                                        )
{
  message("exposure='Mibio','FR02','pathways'")
  message("outcome='LOAD','ADproxy','abeta42','ptau'  ")
    
  input <- inputdf %>% dplyr::filter(exposures == exposure) 
  if(nrow(input) > 0)
  {
    each_exposures <- lapply(seq(1,length(input$name)), function(i) 
    {
      if(exposure == "FR02")
      {
        sp <- fread(input$path[i]) %>% mutate(bac = input$name[i], id = input$name[i], N=5959) 
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
                                    samplesize_col = "N",
                                    gene_col = FALSE,
                                    id_col = "id",
                                    min_pval = 1e-200,
                                    z_col = FALSE,
                                    info_col = FALSE,
                                    chr_col = "hm_chrom",
                                    pos_col = "hm_pos",
                                    log_pval = FALSE
                                    )
      }
      else if (exposure == "Mibio")
      {
        sp <- fread(input$path[i]) %>% mutate(id = input$name[i]) 
        sp$P.weightedSumZ <- as.numeric(sp$P.weightedSumZ)
        sp <- sp %>% as.data.frame() %>% arrange(., P.weightedSumZ) 
        
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
      }
      else if (exposure == "pathways")
      {
        sp <- fread(input$path[i]) %>% mutate( pathway = input$name[i]) 
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
      }
      return(sp_exp_dat)
    })
    return(each_exposures)
  }
  else
  {
    return(NULL)
  }
}
