library(tidyverse)
library(openxlsx)
library(TwoSampleMR)
library(data.table)

#----------------{00 load functions}-----------------------
source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/local_clumb.R")
source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/harmonise_data_modified.R")
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/get_bulk_harmonised_data.R")
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/process_micro_outcome_bulk_harmonsed.R")
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/get_six_UVMR_from_bulk_outcome.R")
#--------------{01.load AD exposure   p < 5e-8}------------------

LOAD <- extract_instruments("ieu-b-2",p1 = 5e-08,p2 = ,clump = FALSE)
ADproxy <- extract_instruments('ebi-a-GCST90012877',p1 = 5e-08,p2 = ,clump = FALSE)

abeta42 <- readRDS("/mnt/data/lijincheng/mGWAS/data/outcome/Jansen_Aβ_tau/final_amyloid_beta42.rds") %>% 
  mutate(id = "CSF amyloid beta42") %>% 
  mutate(p_value = as.numeric(p_value)) %>% 
  dplyr::filter(p_value < 5e-8) %>%
  arrange(.,p_value,decreasing=FALSE) %>% as.data.frame()

ptau <- readRDS("/mnt/data/lijincheng/mGWAS/data/outcome/Jansen_Aβ_tau/final_p_tau.rds") %>% 
  mutate(id = "CSF p-tau") %>%
  mutate(p_value = as.numeric(p_value)) %>% 
  dplyr::filter(p_value < 5e-8) %>%
  arrange(.,p_value,decreasing=FALSE) %>% as.data.frame()

abeta42_exp <-format_data(
  abeta42,
  type = "exposure",
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

abeta42_exp$chr.exposure <- as.character(abeta42_exp$chr.exposure)

ptau_exp <-format_data(
  ptau,
  type = "exposure",
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

ptau_exp$chr.exposure <- as.character(ptau_exp$chr.exposure)


AD_exposres <- bind_rows(LOAD,ADproxy,abeta42_exp,ptau_exp)

AD_clumped_exposures <-local_clumb(dat= AD_exposres, clump_kb = 10000, clump_r2 = 0.001, pop = "EUR") 

save(AD_clumped_exposures,file="/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_clumped_exposures.Rdata") 

#-----------------{02.load outcomes (microbiome, circulating biomarkers)}--------------------------------------
load("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_clumped_exposures.Rdata")
##microbiome
micro_exp <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/exposures_mediators_paths.xlsx", sheet =1)

AD_sp_harmonised <- data.frame()

for(i in 1:nrow(micro_exp)){
  if(micro_exp$exposures[i] == "FR02"){
    sp <- fread(micro_exp$path[i],fill = TRUE) %>% mutate(bac = micro_exp$name[i], id = micro_exp$name[i]) 
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
  }else if (micro_exp$exposures[i] == "Mibio"){
    sp <- fread(micro_exp$path[i]) %>% mutate(id = micro_exp$name[i])
    sp$P.weightedSumZ <- as.numeric(sp$P.weightedSumZ)
    sp <- sp %>% as.data.frame()%>% arrange(., P.weightedSumZ) 
    
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
  }else if (micro_exp$exposures[i] == "pathways"){
    sp <- fread(micro_exp$path[i]) %>% mutate( pathway = micro_exp$exposure[i]) 
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
  
  tmp_AD_sp <- get_bulk_harmonised_data(exposure_dat = AD_clumped_exposures, outcome_dat= sp_out_dat)
  
  AD_sp_harmonised <- bind_rows(AD_sp_harmonised,tmp_AD_sp )
}

#micro_exp <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/exposures_mediators_paths.xlsx", sheet =1)
#load("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_clumped_exposures.Rdata")
#source("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/process_micro_outcome_bulk_harmonsed.R")
#AD_sp_harmonised_list <- lapply(1:nrow(micro_exp), process_micro_outcome_bulk_harmonsed )
#AD_sp_harmonised <- do.call(rbind, AD_sp_harmonised_list)

##FR02
micro_exp <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/exposures_mediators_paths.xlsx", sheet =1)
load("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_clumped_exposures.Rdata")
AD_FR_harmonsed  <- data.frame() 
micro_exp <- micro_exp %>% dplyr::filter(exposures == "FR02")
for(i in 1:nrow(micro_exp)){
    sp <- fread(micro_exp$path[i],fill = TRUE) %>% mutate(bac = micro_exp$name[i], id = micro_exp$name[i]) 
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
    
    sp_out_dat$chr.outcome <- as.character(sp_out_dat$chr.outcome)
    tmp_AD_sp <- get_bulk_harmonised_data(exposure_dat = AD_clumped_exposures, outcome_dat= sp_out_dat)
    AD_FR_harmonsed <- bind_rows(AD_FR_harmonsed,tmp_AD_sp )
}

sp <- readRDS("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/GCST90032474.rds") %>% mutate(bac = micro_exp$name[i], id = micro_exp$name[i]) 
    sp$p_value <- as.numeric(sp$p_value)
    sp <- sp %>% as.data.frame() %>% arrange(., p_value)
    
    sp_out_dat <- TwoSampleMR::format_data(
      dat = sp,
      type = "outcome",
      snps = NULL,
      header = TRUE,
      phenotype_col = "bac",
      snp_col = "rsid",
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
    
    sp_out_dat$chr.outcome <- as.character(sp_out_dat$chr.outcome)
    tmp_AD_sp <- get_bulk_harmonised_data(exposure_dat = AD_clumped_exposures, outcome_dat= sp_out_dat)
    AD_FR_harmonsed <- bind_rows(AD_FR_harmonsed,tmp_AD_sp )



save(AD_FR_harmonsed,file="/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_FR_harmonsed.Rdata")












##mibio
micro_exp <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/exposures_mediators_paths.xlsx", sheet =1)
load("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_clumped_exposures.Rdata")
AD_mibio_harmonsed  <- data.frame() 
micro_exp <- micro_exp %>% dplyr::filter(exposures == "Mibio")
for(i in 1:nrow(micro_exp)){
    sp <- fread(micro_exp$path[i]) %>% mutate(id = micro_exp$name[i])
    sp$P.weightedSumZ <- as.numeric(sp$P.weightedSumZ)
    sp <- sp %>% as.data.frame()%>% arrange(., P.weightedSumZ) 
    
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
    
    sp_out_dat$chr.outcome <- as.character(sp_out_dat$chr.outcome)
    tmp_AD_sp <- get_bulk_harmonised_data(exposure_dat = AD_clumped_exposures, outcome_dat= sp_out_dat)
    AD_mibio_harmonsed <- bind_rows(AD_mibio_harmonsed,tmp_AD_sp )
}
save(AD_mibio_harmonsed,file="/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_mibio_harmonsed.Rdata")


##pathways
micro_exp <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/exposures_mediators_paths.xlsx", sheet =1)
load("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_clumped_exposures.Rdata")
AD_pathways_harmonsed  <- data.frame() 
micro_exp <- micro_exp %>% dplyr::filter(exposures == "pathways")
for(i in 1:nrow( micro_exp )){
    sp <- fread(micro_exp$path[i]) %>% mutate( pathway = micro_exp$exposure[i]) 
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
    
    sp_out_dat$chr.outcome <- as.character(sp_out_dat$chr.outcome)
    tmp_AD_sp <- get_bulk_harmonised_data(exposure_dat = AD_clumped_exposures, outcome_dat= sp_out_dat)
    AD_pathways_harmonsed <- bind_rows(AD_pathways_harmonsed,tmp_AD_sp )
}

save(AD_pathways_harmonsed,file="/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_pathways_harmonsed.Rdata")

AD_sp_harmonised <- bind_rows(AD_FR_harmonsed,AD_mibio_harmonsed,AD_pathways_harmonsed)

#micro_exp[which(micro_exp$exposures == "FR02"),"name"]
#micro_exp[which(micro_exp$exposures == "Mibio"),"name"]
#micro_exp[which(micro_exp$exposures == "pathways"),"name"]

AD_FR_harmonsed$samplesize.outcome <- ifelse(
  is.na(AD_FR_harmonsed$samplesize.outcome),
  5959,
  AD_FR_harmonsed$samplesize.outcome
)


AD_sp_harmonised$samplesize.outcome <- ifelse(
  is.na(AD_sp_harmonised$samplesize.outcome) & AD_sp_harmonised$id.outcome %in% micro_exp[micro_exp$exposures == "FR02", "name"],
  5959,
  AD_sp_harmonised$samplesize.outcome
)

AD_sp_harmonised$id.outcome <- AD_sp_harmonised$outcome

save(AD_sp_harmonised,file="/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_sp_harmonised.Rdata")




##circulating biomarkers

load(file="/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_clumped_exposures.Rdata")
immune <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/lindbohm_immune_BBB_127.xlsx",sheet=4) 
BBB <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/lindbohm_immune_BBB_127.xlsx",sheet=5) 
immune_id <- unique(immune$id.exposure) #1140
BBB_id <- unique(BBB$id.exposure) #687

circulating_exp <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/exposures_mediators_paths.xlsx", sheet =2) %>% 
  dplyr::distinct() %>%
  dplyr::mutate(exposure_type = ifelse(id %in% immune_id, "immune",
                                      ifelse(id %in% BBB_id, "BBB", "metabolites") ))

#write.csv(circulating_exp,"/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/circulating_exp.csv",row.names=F)

metabolites_id <- circulating_exp %>% dplyr::filter(exposure_type == "metabolites") %>% pull("id")

#circulating_biomarkers_gwas <- c(immune_id,BBB_id,metabolites_id) %>% as.data.frame()



metabolites_data <- extract_outcome_data(
  snps = AD_clumped_exposures$SNP, outcomes = metabolites_id)

save(metabolites_data,file="/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_metabolites_outcome.Rdata")


BBB_data <- extract_outcome_data(
      snps = AD_clumped_exposures$SNP, outcomes = BBB_id)

immune_data <- extract_outcome_data(
      snps = AD_clumped_exposures$SNP, outcomes = immune_id)

save(immune_data,file="/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_immune_outcome.Rdata")      
      

repeat{
  try({
    BBB_immune_data <- extract_outcome_data(
      snps = AD_clumped_exposures$SNP, outcomes = c(BBB_id,immune_id))

    })
    if(exists("BBB_immune_data")) break
    Sys.sleep(2)
}

save(BBB_immune_data,file="/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_BBB_immune_outcome.Rdata")

circulating_outcome <- bind_rows(BBB_immune_data,metabolites_data)

repeat{
  try({
   AD_circulating_harmonised <- TwoSampleMR::harmonise_data(
      exposure_dat = AD_clumped_exposures,
      outcome_dat = circulating_outcome,
      action = 2
      )

    })
    if(exists("AD_circulating_harmonised")) break
    Sys.sleep(2)
}
#AD_circulating_harmonised <- TwoSampleMR::harmonise_data(
#  exposure_dat = AD_clumped_exposures,
#  outcome_dat = circulating_biomarkers_data,
#  action = 2
#)

AD_sp_circulating <- bind_rows(AD_sp_harmonised,AD_circulating_harmonised)

save(AD_circulating_harmonised,file="/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_circulating_harmonised.Rdata")

AD_circulating_harmonised[which(AD_circulating_harmonised$id.outcome == "met-d-GlycA"),"samplesize.outcome"] <- 115078
AD_circulating_harmonised[which(AD_circulating_harmonised$id.outcome == "met-d-Sphingomyelins"),"samplesize.outcome"] <- 114999

#-------------------{03 UVMR }-------------------------
##------------{circulating biomarker}-------
load(file="/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_circulating_harmonised.Rdata")
###------------------{LOAD - circulating biomarker}-------------------------

# met-d-GlycA  115078
# met-d-Sphingomyelins  	114999

LOAD_circulating_harmonsed <- AD_circulating_harmonised %>% dplyr::filter(exposure == "Alzheimer's disease || id:ieu-b-2")
ADproxy_circulating_harmonsed <- AD_circulating_harmonised %>% dplyr::filter(exposure == "Alzheimer's disease or family history of Alzheimer's disease || id:ebi-a-GCST90012877")
abeta_circulating_harmonsed <- AD_circulating_harmonised %>% dplyr::filter(exposure == "Cerebrospinal fluid amyloid beta 42 levels")
ptau_circulating_harmonsed <- AD_circulating_harmonised %>% dplyr::filter(exposure == "Cerebrospinal fluid p-tau levels")

#LOAD_circulating_harmonsed <- add_samplesize_info(dat=LOAD_circulating_harmonsed,samplesize="samplesize.outcome")
LOAD_circulating_total_overview <- get_six_UVMR_from_bulk_outcome2(dat = LOAD_circulating_harmonsed, 
                                              exposure = "LOAD",
                                              outcome = "circulating",
                                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/")


#ADproxy_circulating_harmonsed <- add_samplesize_info(dat=ADproxy_circulating_harmonsed,samplesize="samplesize.outcome")
ADproxy_circulating_total_overview <- get_six_UVMR_from_bulk_outcome2(dat = ADproxy_circulating_harmonsed, 
                                              exposure = "ADproxy",
                                              outcome = "circulating",
                                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/")
                                              
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/get_six_UVMR_from_bulk_outcome2.R")                                              
ADproxy_circulating_total_overview <- get_six_UVMR_from_bulk_outcome2(dat = ADproxy_circulating_harmonsed, 
                                              exposure = "ADproxy",
                                              outcome = "circulating",
                                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/")



#abeta_circulating_harmonsed <- add_samplesize_info(dat=abeta_circulating_harmonsed,samplesize="samplesize.outcome")
abeta_circulating_total_overview <- get_six_UVMR_from_bulk_outcome2(dat = abeta_circulating_harmonsed, 
                                              exposure = "abeta42",
                                              outcome = "circulating",
                                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/")


#ptau_circulating_harmonsed <- add_samplesize_info(dat=ptau_circulating_harmonsed,samplesize="samplesize.outcome")
ptau_circulating_total_overview <- get_six_UVMR_from_bulk_outcome2(dat = ptau_circulating_harmonsed, 
                                              exposure = "ptau",
                                              outcome = "circulating",
                                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/")

##------------{gut microbiome}-------
load(file="/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_sp_harmonised.Rdata")
###------------------{LOAD - circulating biomarker}-------------------------
LOAD_sp_harmonsed <- AD_sp_harmonised %>% dplyr::filter(exposure == "Alzheimer's disease || id:ieu-b-2")
ADproxy_sp_harmonsed <- AD_sp_harmonised %>% dplyr::filter(exposure == "Alzheimer's disease or family history of Alzheimer's disease || id:ebi-a-GCST90012877")
abeta_sp_harmonsed <- AD_sp_harmonised %>% dplyr::filter(exposure == "Cerebrospinal fluid amyloid beta 42 levels")
ptau_sp_harmonsed <- AD_sp_harmonised %>% dplyr::filter(exposure == "Cerebrospinal fluid p-tau levels")

LOAD_sp_total_overview <- get_six_UVMR_from_bulk_outcome(dat = LOAD_sp_harmonsed, 
                                              exposure = "LOAD",
                                              outcome = "sp",
                                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/")


ADproxy_sp_total_overview <- get_six_UVMR_from_bulk_outcome(dat = ADproxy_sp_harmonsed, 
                                              exposure = "ADproxy",
                                              outcome = "sp",
                                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/")


abeta_sp_total_overview <- get_six_UVMR_from_bulk_outcome(dat = abeta_sp_harmonsed, 
                                              exposure = "abeta42",
                                              outcome = "sp",
                                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/")


ptau_sp_total_overview <- get_six_UVMR_from_bulk_outcome(dat = ptau_sp_harmonsed, 
                                              exposure = "ptau",
                                              outcome = "sp",
                                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/")


openxlsx::write.xlsx(x=list("AD_sp_harmonised"=AD_sp_harmonised),file="/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_sp_harmonised.xlsx")


#-------------------{04 get sig }-------------------------
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/collect_UVMR_res.R")
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/get_sig_UVMR_res.R")
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/get_FDR_UVMR.R")
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/get_bonf_UVMR.R")

micro_exp <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/exposures_mediators_paths.xlsx", sheet =1)

##-------------------{abeta 42}-------------------------
abeta42_circulating_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/abeta42/circulating/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/abeta42/",
                              outfilename = "abeta42_circulating_collect")

abeta42_sp_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/abeta42/sp/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/abeta42/",
                              outfilename = "abeta42_sp_collect")

abeta42_circulating_sig <- get_UVMR_sig(input = abeta42_circulating_collect)
abeta42_sp_sig <- get_UVMR_sig(input = abeta42_sp_collect)
abeta42_sp_fdr <- get_UVMR_mainfdr(input = abeta42_sp_collect)
abeta42_sp_bonf <- get_UVMR_mainbonf(input = abeta42_sp_collect)


##-------------------{LOAD}-------------------------
LOAD_circulating_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/LOAD/circulating/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/LOAD/",
                              outfilename = "LOAD_circulating_collect")

LOAD_sp_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/LOAD/sp/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/LOAD/",
                              outfilename = "LOAD_sp_collect")

LOAD_circulating_sig <- get_UVMR_sig(input = LOAD_circulating_collect)
LOAD_sp_sig <- get_UVMR_sig(input = LOAD_sp_collect)
LOAD_sp_fdr <- get_UVMR_mainfdr(input = LOAD_sp_collect)
LOAD_sp_bonf <- get_UVMR_mainbonf(input = LOAD_sp_collect)

##-------------------{ptau}-------------------------
ptau_circulating_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/ptau/circulating/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/ptau/",
                              outfilename = "ptau_circulating_collect")

ptau_sp_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/ptau/sp/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/ptau/",
                              outfilename = "ptau_sp_collect")

ptau_circulating_sig <- get_UVMR_sig(input = ptau_circulating_collect)
ptau_sp_sig <- get_UVMR_sig(input = ptau_sp_collect)



##-------------------{ADproxy}-------------------------
ADproxy_circulating_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/ADproxy/circulating/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/ADproxy/",
                              outfilename = "ADproxy_circulating_collect")

ADproxy_sp_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/ADproxy/sp/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/ADproxy/",
                              outfilename = "ADproxy_sp_collect")

ADproxy_circulating_sig <- get_UVMR_sig(input = ADproxy_circulating_collect)
ADproxy_sp_sig <- get_UVMR_sig(input = ADproxy_sp_collect)


openxlsx::write.xlsx(x=list("LOAD_circulating_sig"=LOAD_circulating_sig,
                            "LOAD_sp_sig"=LOAD_sp_sig,
                            "ADproxy_circulating_sig"=ADproxy_circulating_sig,
                            "ADproxy_sp_sig"=ADproxy_sp_sig,
                            "abeta42_circulating_sig"=abeta42_circulating_sig ,
                            "abeta42_sp_sig"=abeta42_sp_sig,
                            "ptau_circulating_sig"=ptau_circulating_sig,
                            "ptau_sp_sig"=ptau_sp_sig),
                    file="/mnt/data/lijincheng/mGWAS/result/01UVMR_reverse/AD_circulating_sp_total_sig.xlsx")









































