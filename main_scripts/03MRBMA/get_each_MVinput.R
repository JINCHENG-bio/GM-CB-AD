get_MVinput_for_MRBMA<- function(exposure="",
                                  outcome="",
                                  outpath=""
                                  )
{
  message("exposure='mibio','FR02','pathways','metabolites'  ")
  message("outcome='LOAD','abeta42','ptau'  ")
  message("Author: jincheng li")

  #----------{00.01 load dataset and functions}-------------------
  library("tidyverse")
  library("data.table")
  library("dplyr")
  library("openxlsx")
  library(MendelianRandomization)
  library(TwoSampleMR)
  source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/get_mv_exp_data_withsamepvalue.R")
  source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/mv_harmonise_data_modified.R")
  
    taxa <- bind_rows(read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_taxa_AD_sig.xlsx",sheet = 2), read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_taxa_AD_sig.xlsx",sheet = 3), read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_taxa_AD_sig.xlsx",sheet = 4))

    mibio <- taxa %>% dplyr::filter(source == "MiBioGen") %>% pull("exposure") %>% unique()
    FR02 <- taxa %>% dplyr::filter(source == "FR02") %>% pull("exposure") %>% unique()
    FR02 <- gsub(" ", "_", FR02)
  
    pathways <- bind_rows(read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_pathways_AD_sig.xlsx",sheet=2), read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_pathways_AD_sig.xlsx",sheet=3), read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_pathways_AD_sig.xlsx",sheet=4)) %>% pull("exposure") %>% unique()

    metabolites <- bind_rows( read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_metabolites_AD_sig.xlsx",sheet=2), #read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_metabolites_AD_sig.xlsx",sheet=3),
 read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_metabolites_AD_sig.xlsx",sheet=4)) %>% pull("exposure") %>% unique()
    metabolites <- gsub(" ", "_", metabolites)

    mibio_input <- paste0("/mnt/data/lijincheng/mGWAS/result/02MRBMA/rawdataforBMA/MiBioGen/rawdata/",mibio,".txt.gz")
    FR02_input <- paste0("/mnt/data/lijincheng/mGWAS/result/02MRBMA/rawdataforBMA/FR02/rawdata/",FR02,".tsv.gz")
    pathways_input <- paste0("/mnt/data/lijincheng/mGWAS/result/02MRBMA/rawdataforBMA/Pathways/rawdata/",pathways,".tsv.gz")
    metabolites_input <- paste0("/mnt/data/lijincheng/mGWAS/result/02MRBMA/rawdataforBMA/metabolites/rawdata/", metabolites,".tsv.gz")

    mibio_refdata <- data.frame(mibio=gsub("_"," ",mibio),mibio_input=mibio_input)
    FR02_refdata <- data.frame(FR02=gsub("_"," ",FR02),FR02_input=FR02_input)
    pathways_refdata <- data.frame(pathways=gsub("_"," ",pathways),pathways_input=pathways_input)
    metabolites_refdata <- data.frame(metabolites=gsub("_"," ",metabolites),metabolites_input=metabolites_input)

    mibio_LOAD_feature <- read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_taxa_AD_sig.xlsx",sheet = 2)%>% dplyr::filter(source == "MiBioGen") %>% pull("exposure") 
    mibio_abeta42_feature <- read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_taxa_AD_sig.xlsx",sheet = 3)%>% dplyr::filter(source == "MiBioGen") %>% pull("exposure") 
    mibio_ptau_feature <- read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_taxa_AD_sig.xlsx",sheet = 4)%>% dplyr::filter(source == "MiBioGen") %>% pull("exposure")

    FR02_LOAD_feature <- read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_taxa_AD_sig.xlsx",sheet = 2)%>% dplyr::filter(source == "FR02") %>% pull("exposure")
    FR02_abeta42_feature <- read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_taxa_AD_sig.xlsx",sheet = 3)%>% dplyr::filter(source == "FR02") %>% pull("exposure")
    FR02_ptau_feature <- read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_taxa_AD_sig.xlsx",sheet = 4)%>% dplyr::filter(source == "FR02") %>% pull("exposure")

    pathways_LOAD_feature <- read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_pathways_AD_sig.xlsx",sheet=2) %>% pull("exposure")
    pathways_abeta42_feature <- read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_pathways_AD_sig.xlsx",sheet=3) %>% pull("exposure")
    pathways_ptau_feature <- read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_pathways_AD_sig.xlsx",sheet=4) %>% pull("exposure")

    metabolites_LOAD_feature <- read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_metabolites_AD_sig.xlsx",sheet=2) %>% pull("exposure")
    metabolites_ptau_feature <- read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/main_metabolites_AD_sig.xlsx",sheet=4) %>% pull("exposure")

    mibio_refdata <- mibio_refdata %>% mutate(LOAD= ifelse(mibio %in% mibio_LOAD_feature,"T","F"),
                                          abeta42 =  ifelse(mibio %in% mibio_abeta42_feature,"T","F"),
                                          ptau =  ifelse(mibio %in%mibio_ptau_feature,"T","F")
                                          )

    FR02_refdata <- FR02_refdata %>% mutate(LOAD= ifelse(FR02 %in% FR02_LOAD_feature,"T","F"),
                                          abeta42 =  ifelse(FR02 %in% FR02_abeta42_feature,"T","F"),
                                          ptau =  ifelse(FR02 %in% FR02_ptau_feature,"T","F")
                                          )

    pathways_refdata <- pathways_refdata %>% mutate(LOAD= ifelse(pathways %in% pathways_LOAD_feature,"T","F"),
                                          abeta42 =  ifelse(pathways %in% pathways_abeta42_feature,"T","F"),
                                          ptau =  ifelse(pathways %in% pathways_ptau_feature,"T","F")
                                          )

    metabolites_refdata <- metabolites_refdata %>% mutate(LOAD= ifelse(metabolites %in% metabolites_LOAD_feature,"T","F"),
                                          #abeta42 =  ifelse(metabolites %in% metabolites_abeta42_feature,"T","F"),
                                          ptau =  ifelse(metabolites %in% metabolites_ptau_feature,"T","F")
                                          )

  #-------------{00.02 load functions and outcome datasates}------------------------------
    #LOAD <- fread("/mnt/data/lijincheng/mGWAS/data/outcome/kunkle_load/Kunkle_etal_Stage1_results.txt") %>%  ##build37
    #        as.data.frame() %>%
    #        mutate(Phenotype = "LOAD") %>% 
    #        mutate(p_value = as.numeric(p_value)) %>% 
    #        arrange(.,p_value,decreasing=FALSE)
    abeta42 <- readRDS("/mnt/data/lijincheng/mGWAS/data/outcome/Jansen_Aβ_tau/final_amyloid_beta42.rds") %>% ##build 38
               mutate(id = "CSF amyloid beta42") %>% 
               mutate(p_value = as.numeric(p_value)) %>% 
               arrange(.,p_value,decreasing=FALSE) %>% as.data.frame()

    ptau <- readRDS("/mnt/data/lijincheng/mGWAS/data/outcome/Jansen_Aβ_tau/final_p_tau.rds") %>% ##build 38
            mutate(id = "CSF p-tau") %>%
            mutate(p_value = as.numeric(p_value)) %>% 
            arrange(.,p_value,decreasing=FALSE) %>% as.data.frame()
  
    abeta42_out <-format_data(
      abeta42,
      type = "outcome",
      snps = NULL,
      header = TRUE,
      phenotype_col = "trait",
      snp_col = "rsid",
      beta_col = "BETA",
      se_col = "SE",
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
      z_col = "Zscore",
      info_col = FALSE,
      chr_col = "chromosome",
      pos_col = "base_pair_location",
      log_pval = FALSE
      )

    ptau_out <-format_data(
      ptau,
      type = "outcome",
      snps = NULL,
      header = TRUE,
      phenotype_col = "trait",
      snp_col = "rsid",
      beta_col = "BETA",
      se_col = "SE",
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
      z_col = "Zscore",
      info_col = FALSE,
      chr_col = "chromosome",
      pos_col = "base_pair_location",
      log_pval = FALSE
      )
  
  
  #---------------{01 get MVinput}-----------
  ##--------------{01.01 mibio}-----------------------------
  if(exposure=="mibio")
  {
    mibio_input <- mibio_refdata %>% dplyr::filter(!!sym(outcome) == 'T')
    mibio_list <-list()
    for( i in 1:dim(mibio_input)[1]){
      sp <- fread(mibio_input[,2][i]) %>% mutate(id = mibio_input[,1][i]) 
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
                                  #eaf_col = "effect_allele_frequency",
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
      mibio_list[[i]] <- sp_exp_dat
    }
    
    mibio_exp <- get_mv_exp_data(l_full=mibio_list,
                     min_pval = 1e-200,
                     log_pval = FALSE, 
                     pval_threshold = 1e-05, 
                     clump_r2 = 0.001, 
                     clump_kb = 10000, 
                     harmonise_strictness = 2 )

    saveRDS(mibio_exp,file=paste0(outpath,"mibio_exp","_",outcome,".RDS"))
    
    if(outcome == "LOAD")
    {
      LOAD <- TwoSampleMR::extract_outcome_data(
      snps = unique(mibio_exp$SNP),
      outcomes = "ieu-b-2")
      mibio_LOAD <- TwoSampleMR::mv_harmonise_data(mibio_exp, LOAD)
      saveRDS(mibio_LOAD, file=paste0(outpath,"mibio","_",outcome,"_MVinput",".RDS"))
    }else if(outcome == "abeta42")
    {
      mibio_abeta42 <- get_mv_harmonise_data_modifed(mibio_exp, abeta42_out, harmonise_strictness = 2) 
      saveRDS(mibio_abeta42, file=paste0(outpath,"mibio","_",outcome,"_MVinput",".RDS"))
    }else if(outcome == "ptau")
    {
      mibio_ptau <- get_mv_harmonise_data_modifed(mibio_exp, ptau_out, harmonise_strictness = 2) 
      saveRDS(mibio_ptau, file=paste0(outpath,"mibio","_",outcome,"_MVinput",".RDS"))
    }
    
  }else if(exposure=="FR02")
  { ##--------------{01.02 FR02}-----------------------------
    FR02_input <- FR02_refdata %>% dplyr::filter(!!sym(outcome) == 'T')
    FR_list <-list()
    for( i in 1:dim(FR02_input)[1]){
        sp <- fread(FR02_input[,2][i]) %>% mutate(bac = FR02_input[,1][i], id = FR02_input[,1][i]) 
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
        FR_list[[i]] <- sp_exp_dat
      }
    
    FR_exp <- get_mv_exp_data(l_full=FR_list,
                         min_pval = 1e-200,
                         log_pval = FALSE, 
                         pval_threshold = 1e-05, 
                         clump_r2 = 0.001, 
                         clump_kb = 10000, 
                         harmonise_strictness = 2 )
                         
    saveRDS(FR_exp,file=paste0(outpath,"FR_exp","_",outcome,".RDS"))
    
    if(outcome == "LOAD")
    {
      LOAD_FR <- TwoSampleMR::extract_outcome_data(
                  snps = unique(FR_exp$SNP),
                  outcomes = "ieu-b-2")
      FR_LOAD <- TwoSampleMR::mv_harmonise_data(FR_exp, LOAD_FR)
      saveRDS(FR_LOAD , file=paste0(outpath,"FR","_",outcome,"_MVinput",".RDS"))
      
    }else if(outcome == 'abeta42')
    {
      FR_abeta42 <- get_mv_harmonise_data_modifed(FR_exp, abeta42_out, harmonise_strictness = 2) 
      saveRDS(FR_abeta42, file=paste0(outpath,"FR","_",outcome,"_MVinput",".RDS"))
    }else if(outcome == "ptau")
    {
      FR_ptau <- get_mv_harmonise_data_modifed(FR_exp, ptau_out, harmonise_strictness = 2) 
      saveRDS(FR_ptau, file=paste0(outpath,"FR","_",outcome,"_MVinput",".RDS"))
    }
    
  }else if(exposure=="pathways")
  { ##--------------{01.03 pathways}-----------------------------
    pathways_input <- pathways_refdata %>% dplyr::filter(!!sym(outcome) == 'T')
    pathways_list <-list()
    for( i in 1:dim(pathways_input)[1]){
        pathway <- fread(pathways_input[,2][i]) %>% mutate( pathway = pathways_input[,1][i]) 
        pathway$pval <- as.numeric(pathway$pval)
        pathway <- pathway %>% as.data.frame() %>% arrange(., pval)
        
        pathway_exp_dat <- TwoSampleMR::format_data(
                                      dat = pathway,
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
        pathways_list[[i]] <- pathway_exp_dat
    }
    
    pathways_exp <- get_mv_exp_data(l_full=pathways_list,
                                     min_pval = 1e-200,
                                     log_pval = FALSE, 
                                     pval_threshold = 1e-05, 
                                     clump_r2 = 0.001, 
                                     clump_kb = 10000, 
                                     harmonise_strictness = 2 )
    saveRDS(pathways_exp,file=paste0(outpath,"pathways_exp","_",outcome,".RDS"))
    
    if(outcome == "LOAD")
    {
      LOAD_pathways <- TwoSampleMR::extract_outcome_data(
                      snps = unique(pathways_exp$SNP),
                      outcomes = "ieu-b-2")

      pathways_LOAD <- TwoSampleMR::mv_harmonise_data(pathways_exp, LOAD_pathways)
      saveRDS(pathways_LOAD , file=paste0(outpath,"pathways","_",outcome,"_MVinput",".RDS"))
      
    }else if(outcome == "abeta42")
    {
      pathways_abeta42 <- get_mv_harmonise_data_modifed(pathways_exp, abeta42_out, harmonise_strictness = 2) 
      saveRDS(pathways_abeta42, file=paste0(outpath,"FR","_",outcome,"_MVinput",".RDS"))
    }else if(outcome == "ptau")
    {
      pathways_ptau <- get_mv_harmonise_data_modifed(pathways_exp, ptau_out, harmonise_strictness = 2) 
      saveRDS(pathways_ptau, file=paste0(outpath,"pathways","_",outcome,"_MVinput",".RDS"))
    }
  }else if(exposure=="metabolites")
  { ##--------------{01.04 metabolites}-----------------------------
    metabolites_input <- metabolites_refdata %>% dplyr::filter(!!sym(outcome) == 'T')
    metabolites_list <-list()
    for( i in 1:dim(metabolites_input)[1]){
        metab <- fread(metabolites_input[,2][i]) %>% mutate( metabolites = metabolites_input[,1][i],N=115082,id = metabolites_input[,1][i]) 
        metab$p_value <- as.numeric(metab$p_value)
        metab <- metab %>% as.data.frame() %>% arrange(., p_value)
        
        metabolites_exp_dat <- TwoSampleMR::format_data(
                                      dat = metab,
                                      type = "outcome",
                                      snps = NULL,
                                      header = TRUE,
                                      phenotype_col = "metabolites",
                                      snp_col = "variant_id",
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
                                      chr_col = "chromosome",
                                      pos_col = "base_pair_location",
                                      log_pval = FALSE
                                      )
      metabolites_list[[i]] <- metabolites_exp_dat
    }
  
    metabolites_exp <- get_mv_exp_data(l_full=metabolites_list,
                         min_pval = 1e-200,
                         log_pval = FALSE, 
                         pval_threshold = 5e-08, 
                         clump_r2 = 0.001, 
                         clump_kb = 10000, 
                         harmonise_strictness = 2 )
    saveRDS(metabolites_exp,file=paste0(outpath,"metabolites_exp","_",outcome,".RDS"))
  
  if(outcome == "LOAD")
    {
      LOAD_metabolites <- TwoSampleMR::extract_outcome_data(
        snps = unique(metabolites_exp$SNP),
        outcomes = "ieu-b-2")
      metabolites_LOAD <- TwoSampleMR::mv_harmonise_data(metabolites_exp, LOAD_metabolites)
      saveRDS(metabolites_LOAD ,  file=paste0(outpath,"metabolites","_",outcome,"_MVinput",".RDS"))
    }else if(outcome == "abeta42")
    {
      metabolites_abeta42 <- get_mv_harmonise_data_modifed(metabolites_exp, abeta42_out, harmonise_strictness = 2) 
      saveRDS(metabolites_abeta42, file=paste0(outpath,"metabolites","_",outcome,"_MVinput",".RDS"))
    }else if(outcome == "ptau")
    {
      metabolites_ptau <- get_mv_harmonise_data_modifed(metabolites_exp, ptau_out, harmonise_strictness = 2) 
      saveRDS(metabolites_ptau, file=paste0(outpath,"metabolites","_",outcome,"_MVinput",".RDS"))
    }
  }
    
}


