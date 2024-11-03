library(tidyverse)
library(data.table)

#----------{01 MRBMA for exposures}--------------------
#-------------{00.01 prepare exposures}---------------
##-------------{00.01.01 LOAD}---------------------
local_exposure_df <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposures_mediators_paths.xlsx",sheet=1)
gwas_exposure_df <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposures_mediators_paths.xlsx",sheet=2)


FRLOAD <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/FR_AD_total_sig.xlsx",sheet=1) %>% dplyr::pull("exposure")
mibioLOAD <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/mibio_AD_total_sig.xlsx",sheet=1) %>% dplyr::pull("exposure")
pathwaysLOAD <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/pathwyas_AD_total_sig.xlsx",sheet=1) %>% dplyr::pull("exposure")

LOAD_id <- c(FRLOAD,mibioLOAD,pathwaysLOAD)

LOAD_exposure_df <- local_exposure_df %>% dplyr::filter( name %in% LOAD_id )  #17



##-------------{00.01.02 ADproxy}---------------------
FRADproxy <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/FR_AD_total_sig.xlsx",sheet=2) %>% dplyr::pull("exposure")
mibioADproxy <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/mibio_AD_total_sig.xlsx",sheet=2) %>% dplyr::pull("exposure")
pathwaysADproxy <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/pathwyas_AD_total_sig.xlsx",sheet=2) %>% dplyr::pull("exposure")

ADproxy_id <- c(FRADproxy,mibioADproxy,pathwaysADproxy)

ADproxy_exposure_df <- local_exposure_df %>% dplyr::filter( name %in% ADproxy_id )  #19



##-------------{00.01.03 abeta42}---------------------
FRabeta42 <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/FR_AD_total_sig.xlsx",sheet=3) %>% dplyr::pull("exposure")
mibioabeta42 <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/mibio_AD_total_sig.xlsx",sheet=3) %>% dplyr::pull("exposure")
pathwaysabeta42 <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/pathwyas_AD_total_sig.xlsx",sheet=3) %>% dplyr::pull("exposure")

abeta42_id <- c(FRabeta42,mibioabeta42,pathwaysabeta42)

abeta42_exposure_df <- local_exposure_df %>% dplyr::filter( name %in% abeta42_id )  #14



##-------------{00.01.04 ptau}---------------------
FRptau <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/FR_AD_total_sig.xlsx",sheet=4) %>% dplyr::pull("exposure")
mibioptau <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/mibio_AD_total_sig.xlsx",sheet=4) %>% dplyr::pull("exposure")
pathwaysptau <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/pathwyas_AD_total_sig.xlsx",sheet=4) %>% dplyr::pull("exposure")

ptau_id <- c(FRptau,mibioptau,pathwaysptau)

ptau_exposure_df <- local_exposure_df %>% dplyr::filter( name %in% ptau_id )  #19


save(LOAD_exposure_df,ADproxy_exposure_df,
     abeta42_exposure_df,ptau_exposure_df,
     file="/mnt/data/lijincheng/mGWAS/result/02MRBMA/per_outcome_exposures_df.Rdata")


#-------------{00.02 prepare outcomes}---------------
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


#-------------{00.03 get exposures for MVinput}---------------
source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/get_local_exposures_forMRBMA.R")
source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/get_mv_exp_data_withsamepvalue.R")
source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/mv_harmonise_data_modified.R")
local_exposure_list <- list(LOAD_exposure_df,ADproxy_exposure_df,
                            abeta42_exposure_df,ptau_exposure_df)

exposurelist <- c('Mibio','FR02','pathways')

for (i in seq_along(local_exposure_list)) 
  {
    outpath <- "/mnt/data/lijincheng/mGWAS/result/02MRBMA/final_MRBMA/exposures_MVinput/"
    outcome <- ifelse(i == 1,"LOAD",
                      ifelse(i == 2,"ADproxy",
                             ifelse(i ==3, "abeta42","ptau")))
    input <- local_exposure_list[[i]]
    for(j in exposurelist)
    {
      tmp_exposures <- get_local_exposures_forMRBMA(inputdf=input,
                                                    exposure=j)
      if(is.null(tmp_exposures))next
      else
      {
        tmp_exp <- get_mv_exp_data(l_full=tmp_exposures,
                                   min_pval = 1e-200,
                                   log_pval = FALSE, 
                                   pval_threshold = 1e-05, 
                                   clump_r2 = 0.001, 
                                   clump_kb = 10000, 
                                   harmonise_strictness = 2)
        
        saveRDS(tmp_exp,file=paste0(outpath,outcome,"_",j,"_MVclumped.RDS"))
        
        if(outcome == "LOAD")
        {
          tmp_out <- TwoSampleMR::extract_outcome_data(
            snps = unique(tmp_exp$SNP),
            outcomes = "ieu-b-2")
          out <- TwoSampleMR::mv_harmonise_data(tmp_exp, tmp_out)
        }else if(outcome == "ADproxy")
        {
          tmp_out <- TwoSampleMR::extract_outcome_data(
            snps = unique(tmp_exp$SNP),
            outcomes = 'ebi-a-GCST90012877')
          out <- TwoSampleMR::mv_harmonise_data(tmp_exp,tmp_out)
        }else if(outcome == "abeta42")
        {
          out <- get_mv_harmonise_data_modifed(tmp_exp, abeta42_out, harmonise_strictness = 2) 
        }else if(outcome == "ptau")
        {
          out <- get_mv_harmonise_data_modifed(tmp_exp, ptau_out, harmonise_strictness = 2) 
        }
        saveRDS(out, file=paste0(outpath,outcome,"_",j,"_MVharmonsied",".RDS")) 
      }
    }
  }


#---------{mediators}----------
LOAD_mediators <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/mediator_AD_total_final.xlsx",sheet=1) %>%
                  dplyr::select(c("exposure","exposure_name","exposure_id","type")) %>%
                  dplyr::mutate(outcome = "LOAD")

ADproxy_mediators <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/mediator_AD_total_final.xlsx",sheet=2) %>%
                  dplyr::select(c("exposure","exposure_name","exposure_id","type")) %>%
                  dplyr::mutate(outcome = "ADproxy")

abeta42_mediators <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/mediator_AD_total_final.xlsx",sheet=3) %>%
                  dplyr::select(c("exposure","exposure_name","exposure_id","type")) %>%
                  dplyr::mutate(outcome = "abeta42")
                  
ptau_mediators <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/mediator_AD_total_final.xlsx",sheet=4) %>%
                  dplyr::select(c("exposure","exposure_name","exposure_id","type")) %>%
                  dplyr::mutate(outcome = "ptau")

mediators_df_ref <- bind_rows(LOAD_mediators,ADproxy_mediators,abeta42_mediators,ptau_mediators)
mediators_df_ref$exposure_id <- gsub("id:","",mediators_df_ref$exposure_id)

#mediators_df_ref$type[which(mediators_df_ref$type == "immune cell")] <- "immune_cell"

mediators_type_list <- list("immune cell",
                            "metabolites",
                            #"protein_NC",
                            "protein_Nature"
                            )

outcome_list <- c("LOAD","ADproxy","abeta42","ptau")

for(i in mediators_type_list){
  for(j in outcome_list){
    outpath <- "/mnt/data/lijincheng/mGWAS/result/02MRBMA/final_MRBMA/Mdeiators_MVinput/"
    tmp_id <- mediators_df_ref %>% dplyr::filter(outcome ==j & type == i) %>%pull("exposure_id")
    
    if(length(tmp_id) >0 )
    {
      
      tmp_exp <- NULL
      while (is.null(tmp_exp)) {
        tmp_exp <- TwoSampleMR::mv_extract_exposures(
                    id_exposure = tmp_id,
                    pval_threshold = 5e-08)
        if (is.null(tmp_exp)) {
        message("ao is still empty. Re-trying...")
        Sys.sleep(1)
        }
      }
    
      if(j == "LOAD")
      {
        tmp_out <- TwoSampleMR::extract_outcome_data(
                   snps = tmp_exp$SNP, outcomes = "ieu-b-2")
        out <- TwoSampleMR::mv_harmonise_data(tmp_exp, tmp_out)
      }
      else if(j == "ADproxy")
      {
        tmp_out <- TwoSampleMR::extract_outcome_data(
                   snps = tmp_exp$SNP, outcomes = 'ebi-a-GCST90012877')
        out <- TwoSampleMR::mv_harmonise_data(tmp_exp, tmp_out)
      }
      else if(j == "abeta42")
      {
        out <- get_mv_harmonise_data_modifed(tmp_exp, abeta42_out, harmonise_strictness = 2) 
      }
      else if(j == "patu")
      {
        out <- get_mv_harmonise_data_modifed(tmp_exp, ptau_out, harmonise_strictness = 2) 
      }
      saveRDS(out, file=paste0(outpath,j,"_",i,"_MVharmonsied",".RDS"))
    }
  }
}



#-------------------{2.MRBMA}----------------------
library(mrbma)
library(foreach)
totalinput <- list.files("/mnt/data/lijincheng/mGWAS/result/02MRBMA/final_MRBMA/Mdeiators_MVinput/",pattern="MVharmonsied.RDS$",full.names = T)
fileexist <- list.files("/mnt/data/lijincheng/mGWAS/result/02MRBMA/final_MRBMA/prior0.1/",pattern="_BMA.RDS$",full.names = F)
fileexist <- gsub("_BMA.RDS$", "", basename(fileexist))
total_basename <- gsub("_MVharmonsied.RDS$", "", basename(totalinput))
totalref <- data.frame(totalinput=totalinput,total_basename=total_basename)
finalinput <- totalref %>% dplyr::filter(! total_basename %in% fileexist)

library(foreach)
library(doParallel)

# 设置并行集群
cl <- makeCluster(5)
registerDoParallel(cl)
# 修改foreach循环
foreach(i = 1:dim(finalinput), .packages=c("tidyverse", "mrbma", "openxlsx")) %dopar% {
    input <- finalinput$totalinput[i]
    name <- finalinput$total_basename[i]  
    MVinput <- readRDS(input)
    res_bma <- mr_bma(
        MVinput,
        prior_prob = 0.1,
        prior_sigma = 0.5,
        top = 10,
        remove_outliers = TRUE,
        remove_influential = TRUE,
        calculate_p = TRUE,
        nrepeat = 1000000
    )
    saveRDS(res_bma, file = paste0("/mnt/data/lijincheng/mGWAS/result/02MRBMA/final_MRBMA/prior0.1/", name, "_BMA.RDS"))
    openxlsx::write.xlsx(x = list("best_model" = res_bma$model_best, "rank_factor" = res_bma$mip_table),
                         file = paste0("/mnt/data/lijincheng/mGWAS/result/02MRBMA/final_MRBMA/prior0.1/", name, "_BMA.xlsx"))
}

# 停止并行集群
stopCluster(cl)



####
library(mrbma)
library(foreach)
totalinput <- list.files("/mnt/data/lijincheng/mGWAS/result/02MRBMA/final_MRBMA/exposures_MVinput/",pattern="MVharmonsied.RDS$",full.names = T)
fileexist <- list.files("/mnt/data/lijincheng/mGWAS/result/02MRBMA/final_MRBMA/prior0.1/",pattern="_BMA.RDS$",full.names = F)
fileexist <- gsub("_BMA.RDS$", "", basename(fileexist))
total_basename <- gsub("_MVharmonsied.RDS$", "", basename(totalinput))
totalref <- data.frame(totalinput=totalinput,total_basename=total_basename)
finalinput <- totalref %>% dplyr::filter(! total_basename %in% fileexist)

library(foreach)
library(doParallel)

# 设置并行集群
cl <- makeCluster(5)
registerDoParallel(cl)
# 修改foreach循环
foreach(i = 1:dim(finalinput), .packages=c("tidyverse", "mrbma", "openxlsx")) %dopar% {
    input <- finalinput$totalinput[i]
    name <- finalinput$total_basename[i]  
    MVinput <- readRDS(input)
    res_bma <- mr_bma(
        MVinput,
        prior_prob = 0.1,
        prior_sigma = 0.5,
        top = 10,
        remove_outliers = TRUE,
        remove_influential = TRUE,
        calculate_p = TRUE,
        nrepeat = 1000000
    )
    saveRDS(res_bma, file = paste0("/mnt/data/lijincheng/mGWAS/result/02MRBMA/final_MRBMA/prior0.1/", name, "_BMA.RDS"))
    openxlsx::write.xlsx(x = list("best_model" = res_bma$model_best, "rank_factor" = res_bma$mip_table),
                         file = paste0("/mnt/data/lijincheng/mGWAS/result/02MRBMA/final_MRBMA/prior0.1/", name, "_BMA.xlsx"))
}

# 停止并行集群
stopCluster(cl)






















