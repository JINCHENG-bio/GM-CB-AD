#------{00 load functions and data input}----------
library(data.table)
library(dplyr)
library(openxlsx)
library(MendelianRandomization)
library(TwoSampleMR)
library(stringr)
library(grDevices)
library(readxl)
library(MRcML)
library(tidyverse)
library(openxlsx)
rm(list=ls())


source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/collect_UVMR_res.R")
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/get_sig_UVMR_res.R")

#-----------------{01 FR02}------------------------
FR_LOAD_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/04replication/Kunkle/FR/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/",
                              outfilename = "FR_LOAD_collect")

FR_ADproxy_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR/01exposure_outcome/FR02/final_MRconmix/with_conmix/",
                                        outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/",
                                        outfilename = "FR_ADproxy_collect")

FR_abeta_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR/01exposure_outcome/FR02/results_ADpatho/Abeta/",
                                        outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/",
                                        outfilename = "FR_abeta_collect")

FR_ptau_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR/01exposure_outcome/FR02/results_ADpatho/tau/",
                                        outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/",
                                        outfilename = "FR_ptau_collect")

openxlsx::write.xlsx(x=list("FR_LOAD_collect" = FR_LOAD_collect,
                            "FR_ADproxy_collect" = FR_ADproxy_collect,
                            "FR_abeta_collect" = FR_abeta_collect,
                            "FR_ptau_collect" = FR_ptau_collect),
                     file = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/FR02_AD_total_collect.xlsx" )


FRsp213 <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/species213_inputdf.xlsx",sheet=1) %>% dplyr::pull("shortname")

FR_LOAD_sig <- get_UVMR_sig(input = FR_LOAD_collect) %>% dplyr::filter(exposure %in% FRsp213)

FR_ADproxy_sig <- get_UVMR_sig(input = FR_ADproxy_collect) %>% dplyr::filter(exposure %in% FRsp213)

FR_abeta_sig <- get_UVMR_sig(input = FR_abeta_collect) %>% dplyr::filter(exposure %in% FRsp213)

FR_ptau_sig <- get_UVMR_sig(input = FR_ptau_collect) %>% dplyr::filter(exposure %in% FRsp213)

openxlsx::write.xlsx(x=list("FR_LOAD_sig" = FR_LOAD_sig,
                            "FR_ADproxy_sig" = FR_ADproxy_sig,
                            "FR_abeta_sig" = FR_abeta_sig,
                            "FR_ptau_sig" = FR_ptau_sig),
                     file = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/FR_AD_total_sig.xlsx" )


#------------------{02 MibioGen}-------------------
mibio_LOAD_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/04replication/Kunkle/mibiogen/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/",
                              outfilename = "Mibio_LOAD_collect")

mibio_ADproxy_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR/01exposure_outcome/MiBioGen/final_with_conmix/with_conmix/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/",
                              outfilename = "Mibio_ADproxy_collect")
                              

mibio_abeta_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR/01exposure_outcome/MiBioGen/MiBioGen_ADpathno/Abeta/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/",
                              outfilename = "Mibio_abeta_collect")

mibio_ptau_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR/01exposure_outcome/MiBioGen/MiBioGen_ADpathno/tau/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/",
                              outfilename = "Mibio_ptau_collect")
                              

mibio_LOAD_sig <- get_UVMR_sig(input = mibio_LOAD_collect)

mibio_ADproxy_sig <- get_UVMR_sig(input = mibio_ADproxy_collect)

mibio_abeta_sig <- get_UVMR_sig(input = mibio_abeta_collect)

mibio_ptau_sig <- get_UVMR_sig(input = mibio_ptau_collect)
                              
openxlsx::write.xlsx(x=list("mibio_LOAD_sig" = mibio_LOAD_sig,
                            "mibio_ADproxy_sig" = mibio_ADproxy_sig,
                            "mibio_abeta_sig" = mibio_abeta_sig,
                            "mibio_ptau_sig" = mibio_ptau_sig),
                     file = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/mibio_AD_total_sig.xlsx" )
                              


#------------------{03 pathways}-------------------
pathwyas_LOAD_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/04replication/Kunkle/pathways/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/",
                              outfilename = "pathwyas_LOAD_collect")

pathwyas_ADproxy_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR/01exposure_outcome/pathways/result/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/",
                              outfilename = "pathwyas_ADproxy_collect")

                              
pathwyas_abeta_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR/01exposure_outcome/pathways/results_pathway_ADpathno/Abeta/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/",
                              outfilename = "pathwyas_abeta_collect")

pathwyas_ptau_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR/01exposure_outcome/pathways/results_pathway_ADpathno/tau/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/",
                              outfilename = "pathwyas_ptau_collect")

pathwyas_LOAD_sig <- get_UVMR_sig(input = pathwyas_LOAD_collect)

pathwyas_ADproxy_sig <- get_UVMR_sig(input = pathwyas_ADproxy_collect)

pathwyas_abeta_sig <- get_UVMR_sig(input = pathwyas_abeta_collect)

pathwyas_ptau_sig <- get_UVMR_sig(input = pathwyas_ptau_collect)
                              
openxlsx::write.xlsx(x=list("pathwyas_LOAD_sig" = pathwyas_LOAD_sig,
                            "pathwyas_ADproxy_sig" = pathwyas_ADproxy_sig,
                            "pathwyas_abeta_sig" = pathwyas_abeta_sig,
                            "pathwyas_ptau_sig" = pathwyas_ptau_sig),
                     file = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/pathwyas_AD_total_sig.xlsx" )



#------------------{04 metabolites}-------------------
metabolites_LOAD_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/04replication/Kunkle/metabolites/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/",
                              outfilename = "metabolites_LOAD_collect") %>% dplyr::filter(!grepl("ratio | in",exposure)) 

metabolites_ADproxy_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR/02mediator_outcome/final_with_conmix/with_conmix/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/",
                              outfilename = "metabolites_ADproxy_collect") %>% dplyr::filter(!grepl("ratio | in",exposure)) 
                              
metabolites_abeta_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR/02mediator_outcome/results_metabolits_ADpathno/Abeta/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/",
                              outfilename = "metabolites_abeta_collect") %>% dplyr::filter(!grepl("ratio | in",exposure)) 

metabolites_ptau_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR/02mediator_outcome/results_metabolits_ADpathno/tau/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/",
                              outfilename = "metabolites_ptau_collect") %>% dplyr::filter(!grepl("ratio | in",exposure)) 

openxlsx::write.xlsx(x=list("metabolites_LOAD_collect" = metabolites_LOAD_sig,
                            "metabolites_ADproxy_collect" = metabolites_ADproxy_sig,
                            "metabolites_abeta_collect" = metabolites_abeta_sig,
                            "metabolites_ptau_collect" = metabolites_ptau_sig),
                     file = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/metabolites_AD_filtered_collect.xlsx" )


metabolites_LOAD_sig <- get_UVMR_sig(input = metabolites_LOAD_collect)

metabolites_ADproxy_sig <- get_UVMR_sig(input = metabolites_ADproxy_collect)

metabolites_abeta_sig <- get_UVMR_sig(input = metabolites_abeta_collect)

metabolites_ptau_sig <- get_UVMR_sig(input = metabolites_ptau_collect)
                              
openxlsx::write.xlsx(x=list("metabolites_LOAD_sig" = metabolites_LOAD_sig,
                            "metabolites_ADproxy_sig" = metabolites_ADproxy_sig,
                            "metabolites_abeta_sig" = metabolites_abeta_sig,
                            "metabolites_ptau_sig" = metabolites_ptau_sig),
                     file = "/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/metabolites_AD_total_sig.xlsx" )

















































