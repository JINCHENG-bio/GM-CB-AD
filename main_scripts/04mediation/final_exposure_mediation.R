##file path context 

library(tidyverse)
library(openxlsx)

#-----------{00.00 prepare datainput}--------------
##-----------------{FR}-----------------
FRsp213 <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/species213_inputdf.xlsx",sheet=1)
FRlist <- list.files("/mnt/data/lijincheng/mGWAS/data/exposure/FR02_all/FR02/",pattern=".tsv.gz$",full.name=TRUE)
FRid <- gsub(".tsv.gz", "", basename(FRlist))

FRdf <- data.frame(id = FRid,path = FRlist) 
FR213_df <- dplyr::left_join(FRsp213[,c(1,3)],FRdf,by=c("Study.accession"="id"))
FRdf <- FR213_df %>% dplyr::mutate(exposure= "FR02")
names(FRdf) <- c("id","name","path","exposures")
FRdf <- FRdf %>% dplyr::select(c("exposures","name","id","path"))

##------------------{mibio}-------------------
mibiolist <- list.files("/mnt/data/lijincheng/mGWAS/data/exposure/MiBioGen/final_raw211/",pattern=".summary.txt.gz$",full.name=TRUE)
mibio_df <- data.frame(exposures= "Mibio",name = gsub(".summary.txt.gz", "", basename(mibiolist)),id=gsub(".summary.txt.gz", "", basename(mibiolist)),path = mibiolist) 

#---------------{pathways}-----------------
pathwayslist <- list.files("/mnt/data/lijincheng/Dutch7738/pathways/",pattern=".tsv.gz$",full.name=TRUE)
id <- gsub("GWAS_","",gsub(".tsv.gz", "", basename(pathwayslist)))
pathways_df <- data.frame(exposures= "pathways",name= id, id = id, path = pathwayslist)

local_exposure_df <- bind_rows(FRdf,mibio_df,pathways_df)

##---------------{metabolites}--------------
metabo78 <- read.csv("/mnt/data/lijincheng/mGWAS/result/02MRBMA/metabolite_accession.csv")
names( metabo78) <- c("name","id")
metabo78 <- metabo78 %>% dplyr::filter(! (grepl("ratio ",name) | grepl(" in ",name) ) ) 

metabo78df <- metabo78 %>%
              dplyr::mutate(id = paste0("ebi-a-",`id`)) %>% 
              dplyr::select(c("id","name"))

##-----------{BBB & immune}---------
immune <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/lindbohm_immune_BBB_127.xlsx",sheet=4) 
BBB <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/lindbohm_immune_BBB_127.xlsx",sheet=5) 
immune_id <- unique(immune$id.exposure) #1140
BBB_id <- unique(BBB$id.exposure) #687

BBB_immune_total_df <- bind_rows(BBB,immune) %>% dplyr::select(c("id.exposure","exposure"))
names(BBB_immune_total_df) <- c("id","name")

gwas_exposure_df <- bind_rows(metabo78df,BBB_immune_total_df)

openxlsx::write.xlsx(x=list("local_exposure_df"=local_exposure_df,
                            "gwas_exposure_df"=gwas_exposure_df),
                            file="/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposures_mediators_paths.xlsx")



#-------------{00.01 prepare exposures}---------------
##-------------{00.01.01 LOAD}---------------------
FRLOAD <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/FR_AD_total_sig.xlsx",sheet=1) %>% dplyr::pull("exposure")
mibioLOAD <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/mibio_AD_total_sig.xlsx",sheet=1) %>% dplyr::pull("exposure")
pathwaysLOAD <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/pathwyas_AD_total_sig.xlsx",sheet=1) %>% dplyr::pull("exposure")

LOAD_id <- c(FRLOAD,mibioLOAD,pathwaysLOAD)

LOAD_exposure_df <- local_exposure_df %>% dplyr::filter( name %in% LOAD_id )  #17

source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/get_local_clumped_exposures.R")
LOAD_exposures <- get_local_clumped_exposures(exposures_df = LOAD_exposure_df, p_threhold = 1e-5)

##-------------{00.01.02 ADproxy}---------------------
FRADproxy <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/FR_AD_total_sig.xlsx",sheet=2) %>% dplyr::pull("exposure")
mibioADproxy <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/mibio_AD_total_sig.xlsx",sheet=2) %>% dplyr::pull("exposure")
pathwaysADproxy <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/pathwyas_AD_total_sig.xlsx",sheet=2) %>% dplyr::pull("exposure")

ADproxy_id <- c(FRADproxy,mibioADproxy,pathwaysADproxy)

ADproxy_exposure_df <- local_exposure_df %>% dplyr::filter( name %in% ADproxy_id )  #19

source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/get_local_clumped_exposures.R")
ADproxy_exposures <- get_local_clumped_exposures(exposures_df = ADproxy_exposure_df , p_threhold = 1e-5)


##-------------{00.01.03 abeta42}---------------------
FRabeta42 <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/FR_AD_total_sig.xlsx",sheet=3) %>% dplyr::pull("exposure")
mibioabeta42 <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/mibio_AD_total_sig.xlsx",sheet=3) %>% dplyr::pull("exposure")
pathwaysabeta42 <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/pathwyas_AD_total_sig.xlsx",sheet=3) %>% dplyr::pull("exposure")

abeta42_id <- c(FRabeta42,mibioabeta42,pathwaysabeta42)

abeta42_exposure_df <- local_exposure_df %>% dplyr::filter( name %in% abeta42_id )  #14

source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/get_local_clumped_exposures.R")
abeta42_exposures <- get_local_clumped_exposures(exposures_df = abeta42_exposure_df, p_threhold = 1e-5)

##-------------{00.01.04 ptau}---------------------
FRptau <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/FR_AD_total_sig.xlsx",sheet=4) %>% dplyr::pull("exposure")
mibioptau <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/mibio_AD_total_sig.xlsx",sheet=4) %>% dplyr::pull("exposure")
pathwaysptau <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/pathwyas_AD_total_sig.xlsx",sheet=4) %>% dplyr::pull("exposure")

ptau_id <- c(FRptau,mibioptau,pathwaysptau)

ptau_exposure_df <- local_exposure_df %>% dplyr::filter( name %in% ptau_id )  #19

source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/get_local_clumped_exposures.R")
ptau_exposures <- get_local_clumped_exposures(exposures_df = ptau_exposure_df, p_threhold = 1e-5)

save(LOAD_exposures,ADproxy_exposures,abeta42_exposures,ptau_exposures,file="/mnt/data/lijincheng/mGWAS/result/02MRBMA/mediation_exposure_per_outcome.Rdata")

#-------------{00.02 prepare mediators}---------------
##-------------{00.02.01 LOAD}---------------------
metabolites_LOAD <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/metabolites_AD_total_sig.xlsx",sheet=1) %>% dplyr::pull("exposure")
BBB_LOAD <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/BBB/BBB_AD_total_sig.xlsx",sheet=1) %>% dplyr::pull("exposure")
immune_LOAD <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/immune/immune_LOAD_sig.xlsx",sheet=1) %>% dplyr::pull("exposure")

LOAD_mediators_name <- unique(c(metabolites_LOAD,BBB_LOAD,immune_LOAD)) #17
LOAD_mediators_id<- gwas_exposure_df %>% dplyr::filter(name %in% LOAD_mediators_name)%>% distinct() %>% pull("id") 

LOAD_mediators <- extract_outcome_data(snps = LOAD_exposures$SNP,  outcomes = LOAD_mediators_id, maf_threshold = 0.01)

LOAD_harmonsed_input <- harmonise_data(exposure_dat = LOAD_exposures,outcome_dat = LOAD_mediators)

##-------------{00.02.02 ADproxy}---------------------
#gwas_exposure_df <-  openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposures_mediators_paths.xlsx",sheet=2)
ADproxy_mediators_name <-  openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/mediator_AD_total_final.xlsx",sheet=2) %>% dplyr::pull("exposure")
ADproxy_mediators_exposureid <-  openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/mediator_AD_total_final.xlsx",sheet=2) %>% dplyr::pull("exposure_id")
ADproxy_mediators_id <- gsub("id:","",ADproxy_mediators_exposureid)

##ADproxy_mediators_id <- gwas_exposure_df %>% dplyr::filter(name %in% ADproxy_mediators_name)%>% distinct() %>% pull("id") 

ADproxy_mediators <- TwoSampleMR::extract_outcome_data(snps = ADproxy_exposures$SNP,  outcomes = ADproxy_mediators_id, maf_threshold = 0.01)
ADproxy_harmonsed_input <- harmonise_data(exposure_dat = ADproxy_exposures,outcome_dat = ADproxy_mediators)

##-------------{00.02.03 abeta42}---------------------
#metabolites = 0
BBB_abeta42 <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/BBB/BBB_AD_total_sig.xlsx",sheet=2) %>% dplyr::pull("exposure")
immune_abeta42 <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/immune/immune_abeta42_sig.xlsx",sheet=1) %>% dplyr::pull("exposure")

abeta42_mediators_name <- unique(c(BBB_abeta42,immune_abeta42)) #7
abeta42_mediators_id <- gwas_exposure_df %>% dplyr::filter(name %in% abeta42_mediators_name) %>% distinct() %>% pull("id") 

abeta42_mediators <- extract_outcome_data(snps = abeta42_exposures$SNP,  outcomes = abeta42_mediators_id, maf_threshold = 0.01)

abeta42_harmonsed_input <- harmonise_data(exposure_dat = abeta42_exposures,outcome_dat = abeta42_mediators)



##-------------{00.02.04 ptau}---------------------
metabolites_ptau <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_micro_pathways_metabolite/metabolites_AD_total_sig.xlsx",sheet=4) %>% dplyr::pull("exposure")
BBB_ptau <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/BBB/BBB_AD_total_sig.xlsx",sheet=3) %>% dplyr::pull("exposure")
immune_ptau <- openxlsx::read.xlsx("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/immune/immune_ptau_sig.xlsx",sheet=1) %>% dplyr::pull("exposure")

ptau_mediators_name <- unique(c(metabolites_ptau,BBB_ptau,immune_ptau)) #13
ptau_mediators_id<- gwas_exposure_df %>% dplyr::filter(name %in% ptau_mediators_name)%>% distinct() %>% pull("id") 

ptau_mediators <- extract_outcome_data(snps = ptau_exposures$SNP,  outcomes = ptau_mediators_id, maf_threshold = 0.01)

ptau_harmonsed_input <- harmonise_data(exposure_dat = ptau_exposures,outcome_dat = ptau_mediators)


save(LOAD_harmonsed_input,ADproxy_harmonsed_input, abeta42_harmonsed_input, ptau_harmonsed_input,file="/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposure_mediators_harmonised.Rdata")

#-------------{00.03 UVMR}---------------
##------------{00.03.01 LOAD}----------------
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/add_samplesize_info.R")
source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/get_UVMR_exp_mediator.R")
LOAD_harmonsed_input <- add_samplesize_info(dat=LOAD_harmonsed_input,samplesize="samplesize.outcome")
LOAD_harmonsed_input <- LOAD_harmonsed_input %>%
  mutate(samplesize.exposure = ifelse(is.na(samplesize.exposure), 5959, samplesize.exposure))
LOAD_harmonsed_input_overview <- get_UVMR_exp_mediator(dat = LOAD_harmonsed_input, 
                                                        outcome = "LOAD",
                                                        outpath = "/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposure_mediation/") 
openxlsx::write.xlsx(x=list(">3SNP_harmonsied"=LOAD_harmonsed_input_overview[[1]],
                            ">3SNP_for_MR"=LOAD_harmonsed_input_overview[[2]]),
                    file = "/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposure_mediation/LOAD_harmonsed_input_overview.xlsx")

##------------{00.03.02 ADproxy}----------------
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/add_samplesize_info.R")
source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/get_UVMR_exp_mediator.R")
ADproxy_harmonsed_input <- add_samplesize_info(dat=ADproxy_harmonsed_input,samplesize="samplesize.outcome")
ADproxy_harmonsed_input <- ADproxy_harmonsed_input %>%
  mutate(samplesize.exposure = ifelse(is.na(samplesize.exposure), 5959, samplesize.exposure))
ADproxy_harmonsed_input_overview <- get_UVMR_exp_mediator(dat = ADproxy_harmonsed_input, 
                                                        outcome = "ADproxy",
                                                        outpath = "/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposure_mediation/") 
openxlsx::write.xlsx(x=list(">3SNP_harmonsied"=ADproxy_harmonsed_input_overview[[1]],
                            ">3SNP_for_MR"=ADproxy_harmonsed_input_overview[[2]]),
                    file = "/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposure_mediation/ADproxy_harmonsed_input_overview.xlsx")

##------------{00.03.03 abeta42}----------------
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/add_samplesize_info.R")
source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/get_UVMR_exp_mediator.R")
abeta42_harmonsed_input <- add_samplesize_info(dat=abeta42_harmonsed_input,samplesize="samplesize.outcome")
abeta42_harmonsed_input_overview <- get_UVMR_exp_mediator(dat = abeta42_harmonsed_input, 
                                                        outcome = "abeta42",
                                                        outpath = "/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposure_mediation/") 
                                                        
openxlsx::write.xlsx(x=list(">3SNP_harmonsied"=abeta42_harmonsed_input_overview[[1]],
                            ">3SNP_for_MR"=abeta42_harmonsed_input_overview[[2]]),
                    file = "/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposure_mediation/abeta42_harmonsed_input_overview.xlsx")

##------------{00.03.04 ptau}----------------
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/add_samplesize_info.R")
source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/get_UVMR_exp_mediator.R")
ptau_harmonsed_input <- add_samplesize_info(dat=ptau_harmonsed_input,samplesize="samplesize.outcome")
ptau_harmonsed_input <- ptau_harmonsed_input %>%
  mutate(samplesize.exposure = ifelse(is.na(samplesize.exposure), 5959, samplesize.exposure))
  
ptau_harmonsed_input_overview <- get_UVMR_exp_mediator(dat = ptau_harmonsed_input, 
                                                        outcome = "ptau",
                                                        outpath = "/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposure_mediation/") 

openxlsx::write.xlsx(x=list(">3SNP_harmonsied"=ptau_harmonsed_input_overview[[1]],
                            ">3SNP_for_MR"=ptau_harmonsed_input_overview[[2]]),
                    file = "/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposure_mediation/ptau_harmonsed_input_overview.xlsx")

#------------{00.04-05 collect and get sig}----------------
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/collect_UVMR_res.R")
source("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/get_sig_UVMR_res.R")

##----{LOAD}----
LOAD_exposure_mediator_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposure_mediation/LOAD/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposure_mediation/",
                              outfilename = "LOAD_exposure_mediator_collect")
LOAD_exposure_mediator_sig <- get_UVMR_sig(input = LOAD_exposure_mediator_collect)

##----{ADproxy}----
ADproxy_exposure_mediator_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposure_mediation/ADproxy/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposure_mediation/",
                              outfilename = "ADproxy_exposure_mediator_collect")
ADproxy_exposure_mediator_sig <- get_UVMR_sig(input = ADproxy_exposure_mediator_collect)

##----{abeta42}----
abeta42_exposure_mediator_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposure_mediation/abeta42/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposure_mediation/",
                              outfilename = "abeta42_exposure_mediator_collect")
abeta42_exposure_mediator_sig <- get_UVMR_sig(input = abeta42_exposure_mediator_collect)

##----{ptau}----
ptau_exposure_mediator_collect <- get_UVMR_collected(inputpath = "/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposure_mediation/ptau/",
                              outpath = "/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposure_mediation/",
                              outfilename = "ptau_exposure_mediator_collect")
ptau_exposure_mediator_sig <- get_UVMR_sig(input = ptau_exposure_mediator_collect)

openxlsx::write.xlsx(x=list("LOAD_exposure_mediator_sig"=LOAD_exposure_mediator_sig,
                            "ADproxy_exposure_mediator_sig"=ADproxy_exposure_mediator_sig,
                            "abeta42_exposure_mediator_sig"=abeta42_exposure_mediator_sig,
                            "ptau_exposure_mediator_sig"=ptau_exposure_mediator_sig),
                    file="/mnt/data/lijincheng/mGWAS/result/02MRBMA/exposure_mediation/AD_exposures_mediators_sig.xlsx")














































