get_UVMR_mainbonf <- function(input = input )
{
  message("input = output of get_UVMR_collected function, a dataframe ")
  #significant association criteria:
  ##1. cML as primary analysis for quantify the effect due to its robustness to both correlated and uncorrelated pleiotropic effects and consideration of potential bais of population stratification and sample overlaping, so that fdr_cML-MA-BIC < 0.05;
  ##2. at least 3 methods significant and all 6 methods with same direction;
  ##3. no Horizontal pleiotropy or heterogeneity, so that Q_pval_Inverse variance weighted > 0.05 & Q_pval_MR Egger > 0.05 & globalP > 0.05(if SNP >3)
  ##4. pass direction test, correct_causal_direction == TRUE
  
  library(tidyverse)
  
  or_name <- paste0("or_",c("MR Egger","Weighted median","Inverse variance weighted","Simple mode","Weighted mode","Contamination mixture","cML-MA-BIC"))
  out<- input %>% 
        dplyr::rowwise() %>%
        dplyr::mutate(all_positive = all(c_across(or_name) > 1),
                      all_negative = all(c_across(or_name) < 1))
  out$mainbonf<- p.adjust(out$`pval_cML-MA-BIC`,method="bonferroni")                   
  
  out<- out %>% 
        dplyr::filter(`mainbonf` < 0.05 & n_significant_simpleout >2 & (all_positive == TRUE | all_negative == TRUE) ) %>%
        dplyr::filter(`Q_pval_Inverse variance weighted` > 0.05 & `Q_pval_MR Egger` > 0.05) %>% 
        dplyr::filter(!  ( (!is.na(`globalP`)) & `globalP` < 0.05) ) %>%
        dplyr::filter(`correct_causal_direction` == TRUE)
  return(out)
}