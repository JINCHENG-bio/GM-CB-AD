library("tidyverse")
library("data.table")
library("readxl")
library("dplyr")
library("TwoSampleMR")
library("stringr")
library("openxlsx")
library("grDevices")
library(readxl)
library(dplyr)
library(MRcML)
library(MendelianRandomization)
###loop
library("doParallel")      #加载doParallel包用于之后注册进程
library("foreach")         #导入foreach包

##############################
###########################3##--------------1.6/7种MR
##############################

filelist <- list.files("/mnt/data/lijincheng/Dutch7738/pathways/",pattern = ".tsv.gz")

pathways_ADproxy205 <- function(i){
    # 使用正则表达式提取字符串
    id <-  sub(".*_(.*?)\\.tsv\\.gz", "\\1", filelist[i])
    out_path <- "/mnt/data/lijincheng/mGWAS/result/04replication/Kunkle/pathways/"
    setwd(out_path)
    
    tryCatch({
    input_path <- "/mnt/data/lijincheng/Dutch7738/pathways/"
    input <- fread(paste0(input_path,filelist[i]), header = TRUE,data.table = FALSE)
    dir.create(id)
    out_each <- paste0(out_path,id,"/")
    ### select p < 1e-5
    df <- input %>% dplyr::filter(pval< 1e-5)
    df$pval <- as.numeric(df$pval)
    df <- df[order(df$pval, decreasing = FALSE),]
    
    ### transform format
    pathway_exp_dat <- TwoSampleMR::format_data(
        dat = df,
        type = "exposure",
          snps = NULL,
          header = TRUE,
          phenotype_col = "feature",
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
          id_col = FALSE,
          min_pval = 1e-200,
          z_col = FALSE,
          info_col = FALSE,
          chr_col = "chr",
          pos_col = "pos",
          log_pval = FALSE
          )
          
    ###clumb
    pathway_exp_dat <- TwoSampleMR::clump_data(dat = pathway_exp_dat,
                          clump_r2 = 0.001,
                          clump_kb = 10000, ##10Mb
                          pop = "EUR")
    
    AD_out_dat <-  extract_outcome_data(snps = pathway_exp_dat$SNP, outcomes = 'ieu-b-2', ###kunkle
                                       maf_threshold = 0.01)
     ###-----{2}harmonise-----
    dat <- TwoSampleMR::harmonise_data(exposure_dat = pathway_exp_dat,
                          outcome_dat = AD_out_dat)
    
    dat$r.exposure <- get_r_from_bsen(b = dat$beta.exposure,
                                  dat$se.exposure,
                                 dat$samplesize.exposure)
    
    ###由于综合考虑多种MR方法的结果，所以这里保留满足INSIDE假设的SNP
    dat <- subset(dat,mr_keep)    ###将满足mr_keep假设的snp保留，作为下一步的输入
    
     ###------{3}重要封装函数的读取----- 
        source("/mnt/data/lijincheng/mGWAS/result/01UVMR/03final_exposure_mediator/function/res_MRPRESSO.R")
        source("/mnt/data/lijincheng/mGWAS/result/01UVMR/03final_exposure_mediator/function/seven_mr_res.R")
        
        ###-----{4}主要分析过程------
        ###-----{4.1}MRPRESSO判断是否存在outlier，若存在且distortion失真测试显著，则剔除离群值再进行后续分析------
        ###-----为了确保不存在SNP水平的异质性，在这一步，通过遍历，直到经过MRPRESSO分析之后无离群值outlier之后
        ###再进行后续的分析

        set.seed(123)
        res_presso_step1 = res_MRPRESSO(dat = dat,NbD = 5000,SignifThreshold = 0.05)
        outlier <- res_presso_step1$Outlier[2][[1]]
        
        library(tidyverse)
        dat_input <- dat
        while(!is.na(outlier[1])){
            ###判断直到没有outlier之后跳出循环
            dat_input <- dat_input %>% filter(!SNP %in% outlier )
            set.seed(123)
            tmp <- res_MRPRESSO(dat = dat_input,NbD = 5000,SignifThreshold = 0.05)
            res_presso_step1$globalP  <- as.character(res_presso_step1$globalP)
            tmp$globalP <- as.character(tmp$globalP)
            res_presso_step1$`Causal Estimate`  <- as.character(res_presso_step1$`Causal Estimate`)
            tmp$`Causal Estimate` <- as.character(tmp$`Causal Estimate`)
      
            res_presso_step1 <- tryCatch(
                                     bind_rows(res_presso_step1,tmp),
                            error = function(e) {
                                                 res_presso_step1$`distortionP`  <- as.character(res_presso_step1$`distortionP`)
                                                 tmp$`distortionP` <- as.character(tmp$`distortionP`)
                                                 bind_rows(res_presso_step1,tmp)
                            }
         )
        outlier <- tmp$Outlier[2][[1]]
        }
        
        ###-----{4.2}在无离群值的情况下进行7种MR分析------
        seven_mr_res <- seven_mr_res(dat_input)
        res <- seven_mr_res[[1]]
        CM <- seven_mr_res[[2]]
        cML_MA_BIC <- seven_mr_res[[3]]       
        
    ###-----{5}main_result------
        results <- TwoSampleMR::generate_odds_ratios(res)
        results[which(results$method =="Contamination mixture"),c("or_lci95","or_uci95")] <- c(exp(CM[1,"low"]),exp(CM[1,"up"]))
        # OR
        results$estimate <- paste0(
                            format(round(results$or, 3), nsmall = 2), " (", 
                            format(round(results$or_lci95, 3), nsmall = 2), "-",
                            format(round(results$or_uci95, 3), nsmall = 2), ")")
        resdata <- dat_input
        R2 <- (dat_input$r.exposure)^2
        N <- dat_input$samplesize.exposure
        resdata$F_statistic <- (R2/(1- R2))*(N-2)
        # Assumption 1 and 3
        names(resdata)
        Assumption13 <- subset(resdata,mr_keep==TRUE,
                               select = c("SNP","pval.exposure",
                                          "pval.outcome", "F_statistic",
                                          "mr_keep"))
        
        # -----{6}--Sensitive_analysis------
        res_hete <- TwoSampleMR::mr_heterogeneity(dat_input)
        res_plei <- TwoSampleMR::mr_pleiotropy_test(dat_input)
        res_leaveone <- mr_leaveoneout(dat_input)  # 
        # 

        set.seed(123)
        res_presso = res_MRPRESSO(dat = dat_input,NbD = 5000,SignifThreshold = 0.05)
    
        ##Steiger fltering
        res_dir  <-  directionality_test(dat_input)  ##dat$r.outcome 默认用get_r_from_pn计算
    
        #### Isq statistics
        I2.exposure <- TwoSampleMR::Isq(dat_input$beta.exposure,dat_input$se.exposure)
        I2.outcome <- TwoSampleMR::Isq(dat_input$beta.outcome,dat_input$se.outcome)
        res_isq = data.frame(I2.exposure=I2.exposure,I2.outcome=I2.outcome)
    
        writexl::write_xlsx(x = list(
        "step1_MRPRESSO" = res_presso_step1,
        "main"=results,
        "Assumption13"=Assumption13,
        "pleiotropy"=res_plei,
        "heterogeneity"=res_hete,
        "final_MRPRESSO"=res_presso,
        "leaveone"=res_leaveone,
        "ReverseMR_Steiger_fltering"=res_dir,
        "Isq_statistic" = res_isq,
        "Conmix"=CM,
        "cML_MA_BIC"=cML_MA_BIC,
        "dat_input"=dat_input),
        path =  paste0(out_each,id,".xlsx"))
        
        # -----{7}--Sensitive_analysis_output------
        p1 <- mr_scatter_plot(res, dat_input)
        #p1[[1]]
        pdf(paste0(out_each,id,"_scatter.pdf"))
        print(p1[[1]])
        dev.off()
    
        p11 <- mr_scatter_plot(res[which(res$method != "Simple mode"),], dat_input)
        #p1[[1]]
        pdf(paste0(out_each,id,"_scatter1.pdf"))
        print(p11[[1]])
        dev.off()
        
        res_single <- mr_singlesnp(dat_input)
        p2 <- mr_forest_plot(res_single)
        pdf(paste0(out_each,id,"_forest.pdf"))
        print(p2[[1]])
        dev.off()
        # 
        p3 <- mr_funnel_plot(res_single)
        pdf(paste0(out_each,id,"_funnel.pdf"))
        print(p3[[1]])
        dev.off()
        # 
        res_loo <- mr_leaveoneout(dat_input)
        pdf(paste0(out_each,id,"_leave_one_out.pdf"))
        print(mr_leaveoneout_plot(res_loo))
        dev.off()    
  }, error = function(e) {
    # 处理报错信息
    message(paste("Error in iteration", i, ":", e$message))
  })

}  

system.time({
  cl<- makeCluster(10)      
  registerDoParallel(cl)       #进行进程注册
  mydata1 <- foreach(
              i=1:length(filelist),          #输入等待请求的参数
              .packages = c("tidyverse", "data.table","readxl","dplyr","TwoSampleMR","stringr","openxlsx","grDevices",
                            "MRcML","MendelianRandomization") 
              #多个进程共享的系统环境
  ) %dopar% pathways_ADproxy205(i)
  stopCluster(cl)
})







