###整理结果
library(tidyr)
library(tidyverse)
library(openxlsx)
out_all_res <- data.frame()
out_null_res <- data.frame()

inputpath <- "/mnt/data/lijincheng/mGWAS/result/04replication/Kunkle/FR/"
outlist <- list.files(inputpath)
###
for(i in 1:length(outlist)){
   inputpath <- "/mnt/data/lijincheng/mGWAS/result/04replication/Kunkle/FR/"
    setwd(inputpath)
     if(file.exists(paste0(inputpath,outlist[i],"/",outlist[i],".xlsx"))){
        ### F statistics and Isq
        F_sheet <- read.xlsx(paste0(inputpath,outlist[i],"/",outlist[i],".xlsx"), sheet=3) 
        F <- mean(F_sheet$F_statistic)
        I_sheet <-   read.xlsx(paste0(inputpath,outlist[i],"/",outlist[i],".xlsx"), sheet=9)
        Isq <- I_sheet[1,1]
        F_I <- data.frame(F_statistic=F,Isq=Isq)

        ##7MR
        seven_mr <- read.xlsx(paste0(inputpath,outlist[i],"/",outlist[i],".xlsx"), sheet=2) 
        
        ###存在有的特征小于4个SNP，因而没有做所有的7种方法，因此应当手动添加每一格为"-"
        allmethod<-c("MR Egger","Weighted median","Inverse variance weighted","Simple mode","Weighted mode","Contamination mixture","cML-MA-BIC")
        missmethod <- base::setdiff(allmethod,seven_mr$method)
        if(length(base::setdiff(allmethod,allmethod)) >0){
          miss <- data.frame(
          id.exposure = rep(c(seven_mr$`id.exposure`[1]),time=length(missmethod)),
          id.outcome= rep(c(seven_mr$`id.outcome`[1]),time=length(missmethod)),
          outcome = rep(c(seven_mr$`outcome`[1]),time=length(missmethod)),
          exposure = rep(c(seven_mr$`exposure`[1]),time=length(missmethod)),
          method = missmethod,
          nsnp = NA,
          b = NA,
          se = NA,
          pval = NA,
          lo_ci = NA,
          up_ci = NA,
          or = NA,
          or_lci95 = NA,
          or_uci95 = NA,
          estimate = NA
        )
        
        seven_mr <- bind_rows(seven_mr,miss)
        }else{
            seven_mr=seven_mr
          }
        
        
        
        seven_mr <- seven_mr%>% 
           tidyr::pivot_wider(., id_cols=c("id.outcome","outcome","exposure"),names_from="method",names_vary="slowest",
                        values_from=c("nsnp","b","se","pval")) 
        #head(seven_mr)
    
        pleiotropy <- read.xlsx(paste0(inputpath,outlist[i],"/",outlist[i],".xlsx"), sheet=4) %>% dplyr::select(c("egger_intercept","se","pval"))
    
        heterogeneity1 <- read.xlsx(paste0(inputpath,outlist[i],"/",outlist[i],".xlsx"), sheet=5) %>%
           tidyr::pivot_wider(., id_cols=c("id.outcome","outcome","exposure"),names_from="method",names_vary="slowest",
                        values_from=c("Q","Q_df","Q_pval"))


        MRPRESSO_global <- read.xlsx(paste0(inputpath,outlist[i],"/",outlist[i],".xlsx"), sheet=6) %>% slice_tail(n=1) %>% dplyr::select(c("globalP"))
    
        reverseMR <-  read.xlsx(paste0(inputpath,outlist[i],"/",outlist[i],".xlsx"), sheet=8) %>% dplyr::select(c(5:8))
        tmp_all_res <- bind_cols(seven_mr,F_I,pleiotropy,heterogeneity1[,-c(1:3)],MRPRESSO_global
                             ,reverseMR)  ##%>% dplyr::select(2:3,32:33,4:31,34:47)
        
        out_all_res$globalP <- as.character(out_all_res$globalP)
        tmp_all_res$globalP <- as.character(tmp_all_res$globalP)
        out_all_res <- bind_rows(out_all_res,tmp_all_res)

        }else{
            null_all_res <- data.frame(id=outlist[i],value=NA)
            out_null_res <- bind_rows(out_null_res,null_all_res)
            }
}


out <- out_all_res %>% mutate(n_significant_simplein=rowSums(select(.,
                                                            `pval_MR Egger`, `pval_Weighted median`, `pval_Inverse variance weighted`,
                                                            `pval_Simple mode`,`pval_Weighted mode`,`pval_Contamination mixture`,
                                                            `pval_cML-MA-BIC`) < 0.05),
                                                            n_significant_simpleout=rowSums(select(.,
                                                            `pval_MR Egger`, `pval_Weighted median`, `pval_Inverse variance weighted`,
                                                            `pval_Weighted mode`,`pval_Contamination mixture`,
                                                            `pval_cML-MA-BIC`) < 0.05))
                                                            
dim(out)

na <-c("outcome","exposure","F_statistic","Isq","n_significant_simplein","n_significant_simpleout","nsnp_MR Egger","b_MR Egger","se_MR Egger","pval_MR Egger","nsnp_Weighted median","b_Weighted median","se_Weighted median",
"pval_Weighted median","nsnp_Inverse variance weighted","b_Inverse variance weighted","se_Inverse variance weighted","pval_Inverse variance weighted","nsnp_Simple mode","b_Simple mode","se_Simple mode","pval_Simple mode","nsnp_Weighted mode",
"b_Weighted mode","se_Weighted mode","pval_Weighted mode","nsnp_Contamination mixture","b_Contamination mixture","se_Contamination mixture","pval_Contamination mixture","nsnp_cML-MA-BIC","b_cML-MA-BIC","se_cML-MA-BIC","pval_cML-MA-BIC",
"egger_intercept","se","pval","Q_MR Egger","Q_df_MR Egger","Q_pval_MR Egger","Q_Inverse variance weighted","Q_df_Inverse variance weighted","Q_pval_Inverse variance weighted","globalP","snp_r2.exposure","snp_r2.outcome","correct_causal_direction","steiger_pval")

out <- out %>% dplyr::select(na)
library(writexl)
writexl::write_xlsx(x = list("all_main_res"=out),
    path =  "/mnt/data/lijincheng/mGWAS/result/04replication/Kunkle/kunkle_FR_collect473.xlsx")




