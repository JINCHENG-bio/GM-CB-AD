library(tidyverse)
library(openxlsx)
library(TwoSampleMR)

source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/local_clumb.R")

##-----------{BBB & immune}---------
immune <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/lindbohm_immune_BBB_127.xlsx",sheet=4) 
BBB <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/lindbohm_immune_BBB_127.xlsx",sheet=5) 
immune_id <- unique(immune$id.exposure) #1140
BBB_id <- unique(BBB$id.exposure) #687

BBB_immune_total_df <- bind_rows(BBB,immune) %>% dplyr::select(c("id.exposure","exposure"))
names(BBB_immune_total_df) <- c("id","name")

#------------{00.01 prepare exposures}-----------------
immune_exposures <- extract_instruments(outcomes= immune_id)
BBB_exposures <- extract_instruments(outcomes= BBB_id)

immune_exposures <- local_clumb(immune_exposures, clump_kb = 10000, clump_r2 = 0.001, pop = "EUR") 
BBB_exposures <- local_clumb(BBB_exposures, clump_kb = 10000, clump_r2 = 0.001, pop = "EUR") 

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

#-------------{00.03 harmonise data}------------
##--------------{LOAD}------------
LOAD_immune <- TwoSampleMR::extract_outcome_data(
                 snps = immune_exposures$SNP, outcomes = "ieu-b-2")

immune_LOAD <- TwoSampleMR::harmonise_data(exposure_dat=immune_exposures, outcome_dat=LOAD_immune, action = 2)

LOAD_BBB <- TwoSampleMR::extract_outcome_data(
                 snps = BBB_exposures$SNP, outcomes = "ieu-b-2")

BBB_LOAD <- TwoSampleMR::harmonise_data(exposure_dat=BBB_exposures, outcome_dat=LOAD_BBB, action = 2)

##--------------{ADproxy}------------
ADproxy_immune <- TwoSampleMR::extract_outcome_data(
                 snps = immune_exposures$SNP, outcomes = 'ebi-a-GCST90012877')

immune_ADproxy <- TwoSampleMR::harmonise_data(exposure_dat=immune_exposures, outcome_dat=ADproxy_immune, action = 2)

ADproxy_BBB <- TwoSampleMR::extract_outcome_data(
                 snps = BBB_exposures$SNP, outcomes = 'ebi-a-GCST90012877')
BBB_ADproxy <- TwoSampleMR::harmonise_data(exposure_dat=BBB_exposures, outcome_dat=ADproxy_BBB, action = 2)


##--------------{abeta42}------------
source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/harmonise_data_modified.R")
immune_abeta42 <- harmonise_data_modifed(exposure_dat=immune_exposures, outcome_dat=abeta42_out, r2_thershold = 0.8, pop = "EUR")
BBB_abeta42 <- harmonise_data_modifed(exposure_dat=BBB_exposures, outcome_dat=abeta42_out, r2_thershold = 0.8, pop = "EUR")

##--------------{ptau}------------
source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/harmonise_data_modified.R")
immune_ptau <- harmonise_data_modifed(exposure_dat=immune_exposures, outcome_dat=ptau_out, r2_thershold = 0.8,  pop = "EUR")
BBB_ptau <- harmonise_data_modifed(exposure_dat=BBB_exposures, outcome_dat=ptau_out, r2_thershold = 0.8, pop = "EUR")


#-------------{00.04 UVMR}------------
##------------{immune}-------
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/get_six_UVMR_from_bulk_harmonised.R")
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/add_samplesize_info.R")
immune_LOAD <- add_samplesize_info(dat=immune_LOAD,samplesize="samplesize.exposure")
immune_LOAD_total_overview <- get_six_UVMR_from_bulk_harmonised(dat = immune_LOAD, 
                                              exposure = "immune",
                                              outcome = "LOAD",
                                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/")

source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/get_six_UVMR_from_bulk_harmonised.R")
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/add_samplesize_info.R")
immune_ADproxy <- add_samplesize_info(dat=immune_ADproxy,samplesize="samplesize.exposure")
immune_ADproxy_total_overview <- get_six_UVMR_from_bulk_harmonised(dat = immune_ADproxy, 
                                              exposure = "immune",
                                              outcome = "ADproxy",
                                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/")


source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/get_six_UVMR_from_bulk_harmonised.R")
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/add_samplesize_info.R")
immune_abeta42 <- add_samplesize_info(dat=immune_abeta42,samplesize="samplesize.exposure")
immune_abeta42_total_overview <- get_six_UVMR_from_bulk_harmonised(dat = immune_abeta42, 
                                              exposure = "immune",
                                              outcome = "abeta42",
                                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/")

source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/get_six_UVMR_from_bulk_harmonised.R")
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/add_samplesize_info.R")
immune_ptau <- add_samplesize_info(dat=immune_ptau,samplesize="samplesize.exposure")
immune_ptau_total_overview <- get_six_UVMR_from_bulk_harmonised(dat = immune_ptau, 
                                              exposure = "immune",
                                              outcome = "ptau",
                                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/")

##------------{BBB}-------
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/get_six_UVMR_from_bulk_harmonised.R")
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/add_samplesize_info.R")
BBB_LOAD <- add_samplesize_info(dat=BBB_LOAD,samplesize="samplesize.exposure")
BBB_LOAD_total_overview <- get_six_UVMR_from_bulk_harmonised(dat = BBB_LOAD, 
                                              exposure = "BBB",
                                              outcome = "LOAD",
                                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/")

source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/get_six_UVMR_from_bulk_harmonised.R")
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/add_samplesize_info.R")
BBB_ADproxy <- add_samplesize_info(dat=BBB_ADproxy,samplesize="samplesize.exposure")
BBB_ADproxy_total_overview <- get_six_UVMR_from_bulk_harmonised(dat = BBB_ADproxy, 
                                              exposure = "BBB",
                                              outcome = "ADproxy",
                                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/")


source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/get_six_UVMR_from_bulk_harmonised.R")
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/add_samplesize_info.R")
BBB_abeta42 <- add_samplesize_info(dat=BBB_abeta42,samplesize="samplesize.exposure")
BBB_abeta42_total_overview <- get_six_UVMR_from_bulk_harmonised(dat = BBB_abeta42, 
                                              exposure = "BBB",
                                              outcome = "abeta42",
                                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/")

source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/get_six_UVMR_from_bulk_harmonised.R")
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/add_samplesize_info.R")
BBB_ptau <- add_samplesize_info(dat=BBB_ptau,samplesize="samplesize.exposure")
BBB_ptau_total_overview <- get_six_UVMR_from_bulk_harmonised(dat = BBB_ptau, 
                                              exposure = "BBB",
                                              outcome = "ptau",
                                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/")


#-------------{00.05 get sig}------------
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/collect_UVMR_res.R")
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/get_sig_UVMR_res.R")
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/get_FDR_UVMR.R")
#source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/get_bonf_UVMR.R")

##immune-LOAD
immune_LOAD_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/immune/LOAD/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/immune/",
                              outfilename = "immune_LOAD_collect")

immune_LOAD_sig <- get_UVMR_sig(input = immune_LOAD_collect)

immune_LOAD_FDR <- get_UVMR_mainfdr(input = immune_LOAD_sig)  


##immune-ADproxy
immune_ADproxy_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/immune/ADproxy/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/immune/",
                              outfilename = "immune_ADproxy_collect")

immune_ADproxy_sig <- get_UVMR_sig(input = immune_ADproxy_collect)
immune_ADproxy_FDR <- get_UVMR_mainfdr(input = immune_ADproxy_collect) 
 

##immune-abeta42
immune_abeta42_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/immune/abeta42/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/immune/",
                              outfilename = "immune_abeta42_collect")

immune_abeta42_sig <- get_UVMR_sig(input = immune_abeta42_collect)
immune_abeta42_FDR <- get_UVMR_mainfdr(input = immune_abeta42_collect) 

##immune-ptau
immune_ptau_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/immune/ptau/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/immune/",
                              outfilename = "immune_ptau_collect")

immune_ptau_sig <- get_UVMR_sig(input = immune_ptau_collect)

openxlsx::write.xlsx(x=list("immune_LOAD_sig"=immune_LOAD_sig,
                            "immune_ADproxy_sig"=immune_ADproxy_sig,
                            "immune_abeta42_sig"=immune_abeta42_sig,
                            "immune_ptau_sig"=immune_ptau_sig),
                            file="/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/immune/immune_AD_sig.xlsx")

openxlsx::write.xlsx(x=list("immune_LOAD_collect"=immune_LOAD_collect,
                            "immune_ADproxy_collect"=immune_ADproxy_collect,
                            "immune_abeta42_collect"=immune_abeta42_collect,
                            "immune_ptau_collect"=immune_ptau_collect),
                            file="/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/immune/immune_AD_collect.xlsx")


BBB_LOAD_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/BBB/LOAD/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/BBB/",
                              outfilename = "BBB_LOAD_collect")

BBB_LOAD_sig <- get_UVMR_sig(input = BBB_LOAD_collect)


BBB_ADproxy_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/BBB/ADproxy/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/BBB/",
                              outfilename = "BBB_ADproxy_collect")

BBB_ADproxy_sig <- get_UVMR_sig(input = BBB_ADproxy_collect)

BBB_abeta42_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/BBB/abeta42/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/BBB/",
                              outfilename = "BBB_abeta42_collect")

BBB_abeta42_sig <- get_UVMR_sig(input = BBB_abeta42_collect)

BBB_ptau_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/BBB/ptau/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/BBB/",
                              outfilename = "BBB_ptau_collect")

BBB_ptau_sig <- get_UVMR_sig(input = BBB_ptau_collect)

openxlsx::write.xlsx(x=list("BBB_LOAD_sig"=BBB_LOAD_sig,
                            "BBB_ADproxy_sig"=BBB_ADproxy_sig,
                            "BBB_abeta42_sig"=BBB_abeta42_sig,
                            "BBB_ptau_sig"=BBB_ptau_sig),
                            file="/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/BBB/BBB_AD_sig.xlsx")
                            
openxlsx::write.xlsx(x=list("BBB_LOAD_collect"=BBB_LOAD_collect,
                            "BBB_ADproxy_collect"=BBB_ADproxy_collect,
                            "BBB_abeta42_collect"=BBB_abeta42_collect,
                            "BBB_ptau_collect"=BBB_ptau_collect),
                            file="/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/BBB/BBB_AD_collect.xlsx")
             
             
#-------------{00.05 get FDR}------------
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/get_sig_UVMR2.R")            
immune_LOAD_collect <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/immune/immune_AD_collect.xlsx",sheet=1)             
             
immune_LOAD_FDR <- get_UVMR_mainfdr(input = immune_LOAD_collect)             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
                            
                            
                            
