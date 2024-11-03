/lijincheng/mGWAS/result/02MRBMA/hyprcoloc/")

rm(list = ls())
library(tidyverse)
library(openxlsx)

#get the significant me

##get significant result
LOAD_med <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/mediation_MR/mirobe_circulator_AD.xlsx",sheet = 1)
LOAD_med_sig <- LOAD_med %>% dplyr::filter( ((EO_indirect_or_lci95 < 1) & (EO_indirect_or_uci95 < 1)) | ((EO_indirect_or_lci95 > 1) & (EO_indirect_or_uci95 > 1))  )

ADproxy_med <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/mediation_MR/mirobe_circulator_AD.xlsx",sheet = 2)
ADproxy_med_sig <- ADproxy_med %>% dplyr::filter( ((EO_indirect_or_lci95 < 1) & (EO_indirect_or_uci95 < 1)) | ((EO_indirect_or_lci95 > 1) & (EO_indirect_or_uci95 > 1))  )

abeta42_med <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/mediation_MR/mirobe_circulator_AD.xlsx",sheet = 3)
abeta42_med_sig <- abeta42_med %>% dplyr::filter( ((EO_indirect_or_lci95 < 1) & (EO_indirect_or_uci95 < 1)) | ((EO_indirect_or_lci95 > 1) & (EO_indirect_or_uci95 > 1))  )

ptau_med <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/mediation_MR/mirobe_circulator_AD.xlsx",sheet = 4)
ptau_med_sig <- ptau_med %>% dplyr::filter( ((EO_indirect_or_lci95 < 1) & (EO_indirect_or_uci95 < 1)) | ((EO_indirect_or_lci95 > 1) & (EO_indirect_or_uci95 > 1))  )


## collect sig
AD_med_sig_total <- bind_rows(LOAD_med_sig,ADproxy_med_sig,abeta42_med_sig,ptau_med_sig)

## same direction

AD_med_coloc <- AD_med_sig_total %>% dplyr::filter(med_prop > 0)
openxlsx::write.xlsx(x=list("AD_med_sig_total" = AD_med_sig_total,
                            "AD_med_coloc" = AD_med_coloc),
                     file = "/mnt/data/lijincheng/mGWAS/result/02MRBMA/hyprcoloc/AD_med_filtered.xlsx")

# hyprcoloc for med 

AD_med_coloc <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/hyprcoloc/AD_med_filtered.xlsx", sheet = 2)

local_exp_path <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposures_mediators_paths.xlsx",sheet=1)

exp_path_need <- local_exp_path %>% dplyr::filter(name %in% AD_med_coloc$exposure)


source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/mv_harmonise_data_modified.R")
source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/get_mv_exp_data_withsamepvalue.R")

input_df <- dplyr::left_join(AD_med_coloc[,c(1:3)],exp_path_need[,c(1,2,4)],by=c("exposure"="name"))

library(data.table)
library(openxlsx)
library(MendelianRandomization)
library(TwoSampleMR)
library(hyprcoloc)


abeta42 <- readRDS("/mnt/data/lijincheng/mGWAS/data/outcome/Jansen_Aβ_tau/final_amyloid_beta42.rds") %>% 
  mutate(id = "CSF amyloid beta42") %>% 
  mutate(p_value = as.numeric(p_value)) %>% 
  arrange(.,p_value,decreasing=FALSE) %>% as.data.frame()

ptau <- readRDS("/mnt/data/lijincheng/mGWAS/data/outcome/Jansen_Aβ_tau/final_p_tau.rds") %>% 
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

med_path <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/hyprcoloc/med_full_dat.xlsx",sheet=1)
outpath <- "/mnt/data/lijincheng/mGWAS/result/02MRBMA/hyprcoloc/MVinput_res/"
for (i in 1:nrow(input_df)){
  exp_name <- input_df$exposure[i]
  
  med_name <- strsplit(input_df$mediators[i]," \\|\\| id:")[[1]][1]
  
  out_name <- ifelse(input_df$outcome[i] == "Alzheimer's disease || id:ieu-b-2", "LOAD",
                     ifelse(input_df$outcome[i] == "Alzheimer's disease or family history of Alzheimer's disease || id:ebi-a-GCST90012877","ADproxy",
                           ifelse(input_df$outcome[i] == "Cerebrospinal fluid amyloid beta 42 levels","abeta42","ptau")))
  
  if(input_df$exposures[i] == "FR02"){
    sp <- fread(input_df$path[i]) %>% mutate(bac = input_df$exposure[i], id = input_df$exposure[i]) 
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
  }else if(input_df$exposures[i] == "Mibio"){
    sp <- fread(input_df$path[i]) %>% mutate(id = input_df$exposure[i]) 
    sp$P.weightedSumZ <- as.numeric(sp$P.weightedSumZ)
    sp <- sp %>% as.data.frame()%>% arrange(., P.weightedSumZ) 
    
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
  }else if(input_df$exposures[i] == "pathways"){
    sp <- fread(input_df$path[i]) %>% mutate( pathway = input_df$exposure[i]) 
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
  sp_exp_dat$chr.outcome <- as.character(sp_exp_dat$chr.outcome)
  
  #med <- strsplit(input_df$mediators[i],"id:")[[1]][2]
  med <- input_df$mediators[i]
  tmp_med_path <- med_path %>% dplyr::filter(mediators == med) %>% pull("med_path")
  if (med %in% c("Programmed cell death 1 ligand 2 || id:prot-a-2215","Reticulon-4 receptor || id:prot-a-2609") )
  {
    mediator <- fread(tmp_med_path)  %>% mutate(n = 3301)
  }else if( med %in% c("Average diameter for LDL particles || id:ebi-a-GCST90092889","Glutamine levels || id:ebi-a-GCST90092818","Remnant cholesterol (non-HDL, non-LDL -cholesterol) || id:ebi-a-GCST90092943"))
  {
    mediator <- fread(tmp_med_path)  %>% mutate(n = 115082)
  }else
  {
    mediator <- fread(tmp_med_path) 
  }
  
  mediator$pval <- as.numeric(mediator$pval)
  mediator <- mediator %>% as.data.frame() %>% arrange(., p_value) %>% mutate(phe = med, id = strsplit(med,"id:")[[1]][2])
  med_type <- med_path %>% dplyr::filter(mediators == med) %>% pull("med_type")
  med_exp_dat <- TwoSampleMR::format_data(
      dat = mediator,
      type = "outcome",
      snps = NULL,
      header = TRUE,
      phenotype_col = "phe",
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
      samplesize_col = "n",
      gene_col = FALSE,
      id_col = "id",
      min_pval = 1e-200,
      z_col = FALSE,
      info_col = FALSE,
      chr_col = "chromosome",
      pos_col = "base_pair_location",
      log_pval = FALSE
    )
  med_exp_dat$chr.outcome <- as.character(med_exp_dat$chr.outcome)
  
  tmp_mv_exp_list <- list(sp_exp_dat,med_exp_dat)
  
  tmp_mv_exp <- get_mv_exp_data(l_full=tmp_mv_exp_list,
                         min_pval = 1e-200,
                         log_pval = FALSE, 
                         pval_threshold = 1, 
                         clump_r2 = 0.001, 
                         clump_kb = 10000, 
                         harmonise_strictness = 2 )
  
  if(out_name == "LOAD")
    {
      LOAD_out <- TwoSampleMR::extract_outcome_data(
        snps = unique(tmp_mv_exp$SNP),
        outcomes = "ieu-b-2")
      tmp_mv_out <- TwoSampleMR::mv_harmonise_data(tmp_mv_exp, LOAD_out)
    }else if(out_name == "ADproxy")
    {
      ADproxy_out <- TwoSampleMR::extract_outcome_data(
        snps = unique(tmp_mv_exp$SNP),
        outcomes = "ebi-a-GCST90012877")
      tmp_mv_out <- TwoSampleMR::mv_harmonise_data(tmp_mv_exp, ADproxy_out) 
    }else if(out_name == "abeta42")
    {
      tmp_mv_out <- get_mv_harmonise_data_modifed(tmp_mv_exp, abeta42_out, harmonise_strictness = 2) 
    }else if(out_name == "ptau")
    {
      tmp_mv_out <- get_mv_harmonise_data_modifed(tmp_mv_exp, ptau_out, harmonise_strictness = 2) 
    }
  saveRDS(tmp_mv_out,  file=paste0(outpath,exp_name,"_",med_name,"_",out_name,"_MVinput",".RDS"))
  
}




med_path <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/hyprcoloc/med_full_dat.xlsx",sheet=1)
outpath <- "/mnt/data/lijincheng/mGWAS/result/02MRBMA/hyprcoloc/MVinput_res/"




totalinput <- input_df$exposure
fileexist <- list.files(outpath,pattern="_MVinput.RDS$",full.names = F)
fileexist <- gsub("_MVinput.RDS$", "", basename(fileexist))
fileexist <- lapply(fileexist, function(x) gsub("_.*", "", x))
fileexist <- unlist(fileexist)


final_input_df <- input_df %>% dplyr:: filter(! (exposure %in% fileexist))
##并行
library(foreach)
library(doParallel)
# 设置并行集群
cl <- makeCluster(5)
registerDoParallel(cl)
# 修改foreach循环
foreach(i = 1:nrow(final_input_df), .packages=c("tidyverse", "data.table", "openxlsx","MendelianRandomization","TwoSampleMR")) %dopar% {
  exp_name <- final_input_df$exposure[i]
  
  med_name <- strsplit(final_input_df$mediators[i]," \\|\\| id:")[[1]][1]
  
  out_name <- ifelse(final_input_df$outcome[i] == "Alzheimer's disease || id:ieu-b-2", "LOAD",
                     ifelse(final_input_df$outcome[i] == "Alzheimer's disease or family history of Alzheimer's disease || id:ebi-a-GCST90012877","ADproxy",
                           ifelse(final_input_df$outcome[i] == "Cerebrospinal fluid amyloid beta 42 levels","abeta42","ptau")))
  
  if(final_input_df$exposures[i] == "FR02"){
    sp <- fread(final_input_df$path[i]) %>% mutate(bac = final_input_df$exposure[i], id = final_input_df$exposure[i]) 
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
  }else if (final_input_df$exposures[i] == "Mibio"){
    sp <- fread(final_input_df$path[i]) %>% mutate(id = final_input_df$exposure[i]) 
    sp$P.weightedSumZ <- as.numeric(sp$P.weightedSumZ)
    sp <- sp %>% as.data.frame()%>% arrange(., P.weightedSumZ) 
    
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
  }else if (final_input_df$exposures[i] == "pathways"){
    sp <- fread(final_input_df$path[i]) %>% mutate( pathway = final_input_df$exposure[i]) 
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
  
  sp_exp_dat$chr.outcome <- as.character(sp_exp_dat$chr.outcome)
  
  #med <- strsplit(final_input_df$mediators[i],"id:")[[1]][2]
  med <- final_input_df$mediators[i]
  tmp_med_path <- med_path %>% dplyr::filter(mediators == med) %>% pull("med_path")
  if (med %in% c("Programmed cell death 1 ligand 2 || id:prot-a-2215","Reticulon-4 receptor || id:prot-a-2609") )
  {
    mediator <- fread(tmp_med_path)  %>% mutate(n = 3301)
  }else if( med %in% c("Average diameter for LDL particles || id:ebi-a-GCST90092889","Glutamine levels || id:ebi-a-GCST90092818","Remnant cholesterol (non-HDL, non-LDL -cholesterol) || id:ebi-a-GCST90092943"))
  {
    mediator <- fread(tmp_med_path)  %>% mutate(n = 115082)
  }else
  {
    mediator <- fread(tmp_med_path) 
  }
  
  mediator$pval <- as.numeric(mediator$pval)
  mediator <- mediator %>% as.data.frame() %>% arrange(., p_value) %>% mutate(phe = med, id = strsplit(med,"id:")[[1]][2])
  med_type <- med_path %>% dplyr::filter(mediators == med) %>% pull("med_type")
  med_exp_dat <- TwoSampleMR::format_data(
      dat = mediator,
      type = "outcome",
      snps = NULL,
      header = TRUE,
      phenotype_col = "phe",
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
      samplesize_col = "n",
      gene_col = FALSE,
      id_col = "id",
      min_pval = 1e-200,
      z_col = FALSE,
      info_col = FALSE,
      chr_col = "chromosome",
      pos_col = "base_pair_location",
      log_pval = FALSE
    )
  med_exp_dat$chr.outcome <- as.character(med_exp_dat$chr.outcome)
  
  tmp_mv_exp_list <- list(sp_exp_dat,med_exp_dat)
  
  
  tmp_mv_exp <- tryCatch(
    get_mv_exp_data(l_full=tmp_mv_exp_list,
                         min_pval = 1e-200,
                         log_pval = FALSE, 
                         pval_threshold = 1, 
                         clump_r2 = 0.001, 
                         clump_kb = 10000, 
                         harmonise_strictness = 2 )
  )
  
  if(out_name == "LOAD")
    {
      LOAD_out <- TwoSampleMR::extract_outcome_data(
        snps = unique(tmp_mv_exp$SNP),
        outcomes = "ieu-b-2")
      tmp_mv_out <- TwoSampleMR::mv_harmonise_data(tmp_mv_exp, LOAD_out)
    }else if(out_name == "ADproxy")
    {
      ADproxy_out <- TwoSampleMR::extract_outcome_data(
        snps = unique(tmp_mv_exp$SNP),
        outcomes = "ebi-a-GCST90012877")
      tmp_mv_out <- TwoSampleMR::mv_harmonise_data(tmp_mv_exp, ADproxy_out) 
    }else if(out_name == "abeta42")
    {
      tmp_mv_out <- get_mv_harmonise_data_modifed(tmp_mv_exp, abeta42_out, harmonise_strictness = 2) 
    }else if(out_name == "ptau")
    {
      tmp_mv_out <- get_mv_harmonise_data_modifed(tmp_mv_exp, ptau_out, harmonise_strictness = 2) 
    }
  saveRDS(tmp_mv_out,  file=paste0(outpath,exp_name,"_",med_name,"_",out_name,"_MVinput",".RDS"))
 
}

# 停止并行集群
stopCluster(cl)














