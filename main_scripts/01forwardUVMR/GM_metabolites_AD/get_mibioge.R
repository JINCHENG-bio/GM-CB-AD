library("tidyverse")
library("data.table")
library("stringr")
library("openxlsx")
library("grDevices")
library(readxl)
library(dplyr)
library(MRcML)
library(MendelianRandomization)
library(TwoSampleMR)
library("doParallel")      #加载doParallel包用于之后注册进程
library("foreach")         #导入foreach包


##02 mibio kunkle

mibio <- fread("/mnt/data/lijincheng/mGWAS/data/exposure/MiBioGen/MBG.allHits.p1e4.txt") %>% 
    group_by(bac)

gut <- unique(mibio$bac)


mibio_Kunkle <- function(i){
    ###-----{1}明确exposure，outcome以及输出路径------
    out <- "/mnt/data/lijincheng/mGWAS/result/04replication/Kunkle/mibiogen/"
    setwd(out)
    id <- gut[i]
    dir.create(id)
    out_each <- paste0(out,gut[i],"/")
    tryCatch({
    microbe <- mibio %>% filter(bac == gut[i]) %>% filter(P.weightedSumZ < 1e-5)
    
    mibio_exp_dat <- TwoSampleMR::format_data(
          dat = microbe,
          type = "exposure",
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
          id_col = FALSE,
          min_pval = 1e-200,
          z_col = "Z.weightedSumZ",
          info_col = FALSE,
          chr_col = "chr",
          pos_col = "bp",
          log_pval = FALSE
    )

    mibio_exp_dat <- TwoSampleMR::clump_data(dat = mibio_exp_dat,
                          clump_r2 = 0.001,
                          clump_kb = 10000, ##10Mb
                          pop = "EUR")
   # if(AD_out[j] == "schwart"){
    AD_out_dat <- extract_outcome_data(snps = mibio_exp_dat$SNP, outcomes = 'ieu-b-2', ###kunkle
                                       maf_threshold = 0.01)
    ###-----{2}harmonise-----
    dat <- TwoSampleMR::harmonise_data(exposure_dat = mibio_exp_dat,
                          outcome_dat = AD_out_dat)
    
    
    dat$r.exposure <- get_r_from_bsen(b = dat$beta.exposure,
                                  dat$se.exposure,
                                 dat$samplesize.exposure)
    ###由于综合考虑多种MR方法的结果，所以这里保留满足INSIDE假设的SNP
    dat <- subset(dat,mr_keep)    ###将满足mr_keep假设的snp保留，作为下一步的输入
    
    ####dat$r.outcome 默认用get_r_from_pn计算
    ### LOAD缺少EAF信息所以用，get_r_from_bsen估计R2
    #dat$r.outcome <- get_r_from_bsen(b = dat$beta.outcome,
    #                                 dat$se.outcome,
    #                                dat$samplesize.outcome)
    
    ###------{3}重要封装函数的读取----- 
    source("/mnt/data/lijincheng/mGWAS/result/01UVMR/03exposure_mediator/res_MRPRESSO.R")
    source("/mnt/data/lijincheng/mGWAS/result/01UVMR/03exposure_mediator/seven_mr_res.R")
    
    ###-----{4}主要分析过程------
    ###-----{4.1}MRPRESSO判断是否存在outlier，若存在且distortion失真测试显著，则剔除离群值再进行后续分析------
    ###-----为了确保不存在SNP水平的异质性，在这一步，通过遍历，直到经过MRPRESSO分析之后无离群值outlier之后
    ###再进行后续的分析

    set.seed(i)
    res_presso_step1 = res_MRPRESSO(dat = dat,NbD = 5000,SignifThreshold = 0.05)
    outlier <- res_presso_step1$Outlier[2][[1]]
        
    library(tidyverse)
    dat_input <- dat
    while(!is.na(outlier[1])){
      ###判断直到没有outlier之后跳出循环
      dat_input <- dat_input %>% filter(!SNP %in% outlier )
      set.seed(i)
      tmp <- res_MRPRESSO(dat = dat_input,NbD = 5000,SignifThreshold = 0.05)
      res_presso_step1$globalP  <- as.character(res_presso_step1$globalP)
      tmp$globalP <- as.character(tmp$globalP)
      res_presso_step1 <- bind_rows(res_presso_step1,tmp)
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

    set.seed(i)
    res_presso = res_MRPRESSO(dat = dat_input,NbD = 5000,SignifThreshold = 0.05)
    
    ##Steiger fltering
    res_dir  <-  directionality_test(dat_input)  ##dat$r.outcome 默认用get_r_from_pn计算
    
    #### Isq statistics
    I2.exposure <- TwoSampleMR::Isq(dat_input$beta.exposure,dat_input$se.exposure)
    I2.outcome <- TwoSampleMR::Isq(dat_input$beta.outcome,dat_input$se.outcome)
    res_isq = data.frame(I2.exposure=I2.exposure,I2.outcome=I2.outcome)
    
    # [["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
    #res_presso[[1]][[2]][[1]][["Pvalue"]]
    #sink(paste0(out_each,id,"_PRESSO.txt"),append=FALSE,split = FALSE) 
    #print(res_presso)
    #sink()
    #print(res_presso)
    
    # ------{4}Merge_results--------
    # Export
    #openxlsx::write.xlsx(x = list(
    #"step1_MRPRESSO" = res_presso_step1,
    #"main"=results,
    #"Assumption13"=Assumption13,
    #"pleiotropy"=res_plei,
    #"heterogeneity"=res_hete,
    #"final_MRPRESSO"=res_presso,
    #"leaveone"=res_leaveone,
    #"ReverseMR_Steiger_fltering"=res_dir,
    #"Isq_statistic" = res_isq,
    #"Conmix"=CM,
    #"cML_MA_BIC"=cML_MA_BIC,
    #"dat_input"=dat_input),
    #overwrite = TRUE,
    #paste0(out_each,id,".xlsx"))
    
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
    
    # -----{5}--Sensitive_analysis------
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
  cl<- makeCluster(5)      
  registerDoParallel(cl)       #进行进程注册
  mydata1 <- foreach(
              i=1:length(gut),          #输入等待请求的参数
              .packages = c("tidyverse", "data.table","readxl","dplyr","TwoSampleMR","stringr","openxlsx","grDevices",
                           "MRcML","MendelianRandomization") 
              #多个进程共享的系统环境
  ) %dopar% mibio_Kunkle(i)
  stopCluster(cl)
})
