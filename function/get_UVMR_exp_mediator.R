get_UVMR_exp_mediator <- function(dat = dat, 
                                  outcome = "",
                                  outpath = "") 
{
  source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/seven_mr_res.R")
  source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/res_MRPRESSO.R")
  message("dat= bulk harmonised data")
  message("dat= bulk harmonised data(multiple exposures and multiple outcomes)")
  
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
  
  dat$r.exposure <- TwoSampleMR::get_r_from_bsen(b = dat$beta.exposure,
                                                 dat$se.exposure,
                                                 dat$samplesize.exposure)
  dat$r.outcome <- TwoSampleMR::get_r_from_bsen(b = dat$beta.outcome,
                                                 dat$se.outcome,
                                                 dat$samplesize.outcome)
  df1 <- dat %>%
    dplyr::group_by(id.exposure,id.outcome) %>%
    dplyr::filter(n() > 2) %>%
    dplyr::ungroup()
  
  ##no enough SNPs
  #df2 <- dat %>%
  #  dplyr::anti_join(df1, by = "id.exposure")
  
  df3 <- data.frame()
  
  id_exposure <- unique(df1$id.exposure)
  id_outcome <- unique(df1$id.outcome)

  each_input <- lapply(id_exposure, function(i) {
    each_outcome <- lapply(id_outcome, function(j) {
      input <- df1 %>% dplyr::filter(id.exposure == i & id.outcome == j) %>% as.data.frame()
      mediator <- j
      exposure <- i
      ###
      input <- subset(input, mr_keep)  ###将满足mr_keep假设的snp保留，作为下一步的输入
    
      set.seed(123)
      res_presso_step1 = res_MRPRESSO(dat = input, NbD = 5000, SignifThreshold = 0.05)
      outlier <- res_presso_step1$Outlier[2][[1]]
      
      dat_input <- input
      while (!is.na(outlier[1])) {
        ###判断直到没有outlier之后跳出循环
        dat_input <- dat_input %>% filter(!SNP %in% outlier)
        set.seed(123)
        tmp <- res_MRPRESSO(dat = dat_input, NbD = 5000, SignifThreshold = 0.05)
        res_presso_step1$globalP <- as.character(res_presso_step1$globalP)
        tmp$globalP <- as.character(tmp$globalP)
      
        res_presso_step1$`Causal Estimate` <- as.character(res_presso_step1$`Causal Estimate`)
        tmp$`Causal Estimate` <- as.character(tmp$`Causal Estimate`)
      
        res_presso_step1 <- bind_rows(res_presso_step1, tmp)
        outlier <- tmp$Outlier[2][[1]]
      }
      
      if (nrow(dat_input) < 3) {
        return(NULL)
      } else {
        out_each <- paste0(outpath, outcome, "/", exposure,"_", mediator, "/")
        dir.create(out_each, recursive = TRUE)
        
       
        seven_mr_res <- seven_mr_res(dat_input)
        res <- seven_mr_res[[1]]
        CM <- seven_mr_res[[2]]
        cML_MA_BIC <- seven_mr_res[[3]]
      
      
        results <- TwoSampleMR::generate_odds_ratios(res)
        results[which(results$method == "Contamination mixture"), c("or_lci95", "or_uci95")] <- c(exp(CM[1, "low"]), exp(CM[1, "up"]))
        # OR
        results$estimate <- paste0(
        format(round(results$or, 3), nsmall = 2), " (",
        format(round(results$or_lci95, 3), nsmall = 2), "-",
        format(round(results$or_uci95, 3), nsmall = 2), ")"
        )
        resdata <- dat_input
        R2 <- (dat_input$r.exposure)^2
        N <- dat_input$samplesize.exposure
        resdata$F_statistic <- (R2 / (1 - R2)) * (N - 2)
        # Assumption 1 and 3
        names(resdata)
        Assumption13 <- subset(
        resdata, mr_keep == TRUE,
          select = c("SNP", "pval.exposure",
                     "pval.outcome", "F_statistic",
                     "mr_keep")
        )
        # -------Sensitive_analysis------
        res_hete <- TwoSampleMR::mr_heterogeneity(dat_input)
        res_plei <- TwoSampleMR::mr_pleiotropy_test(dat_input)
        res_leaveone <- mr_leaveoneout(dat_input)  #
        #
      
        set.seed(123)
        res_presso = res_MRPRESSO(dat = dat_input, NbD = 5000, SignifThreshold = 0.05)
        
        ##Steiger fltering
        res_dir <- directionality_test(dat_input) 
        
        #### Isq statistics
        I2.exposure <- TwoSampleMR::Isq(dat_input$beta.exposure, dat_input$se.exposure)
        I2.outcome <- TwoSampleMR::Isq(dat_input$beta.outcome, dat_input$se.outcome)
        res_isq = data.frame(I2.exposure = I2.exposure, I2.outcome = I2.outcome)
        
        # ------{4}Merge_results--------
        # Export
        writexl::write_xlsx(
          x = list(
            "step1_MRPRESSO" = res_presso_step1,
            "main" = results,
            "Assumption13" = Assumption13,
            "pleiotropy" = res_plei,
            "heterogeneity" = res_hete,
            "final_MRPRESSO" = res_presso,
            "leaveone" = res_leaveone,
            "ReverseMR_Steiger_fltering" = res_dir,
            "Isq_statistic" = res_isq,
            "Conmix" = CM,
            "cML_MA_BIC" = cML_MA_BIC,
            "dat_input" = dat_input
          ),
          path = paste0(out_each, i, ".xlsx")
        )
        
        # -------Sensitive_analysis------
        p1 <- mr_scatter_plot(res, dat_input)
        #p1[[1]]
        pdf(paste0(out_each, exposure,"_", mediator, "_scatter.pdf"))
        print(p1[[1]])
        dev.off()
        
        p11 <- mr_scatter_plot(res[which(res$method != "Simple mode"),], dat_input)
        #p1[[1]]
        pdf(paste0(out_each, exposure,"_", mediator, "_scatter1.pdf"))
        print(p11[[1]])
        dev.off()
        
        res_single <- mr_singlesnp(dat_input)
        p2 <- mr_forest_plot(res_single)
        pdf(paste0(out_each, exposure,"_", mediator, "_forest.pdf"))
        print(p2[[1]])
        dev.off()
        # 
        p3 <- mr_funnel_plot(res_single)
        pdf(paste0(out_each, exposure,"_", mediator, "_funnel.pdf"))
        print(p3[[1]])
        dev.off()
        # 
        res_loo <- mr_leaveoneout(dat_input)
        pdf(paste0(out_each, exposure,"_", mediator, "_leave_one_out.pdf"))
        print(mr_leaveoneout_plot(res_loo))
        dev.off()
        
        return(dat_input)
        }
      })
    dfinput <- bind_rows(each_outcome)
    return(dfinput)
    })
    df3 <- bind_rows(each_input)
  out_list <- list(df1,df3)
  return(out_list)
}
