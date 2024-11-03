get_bulk_expANDout_harmonised_data <- function(exposure_data = bulkexposure, 
                                      outcome_data = bulkoutcome
                                    ) {
  message("exposure_data = bulk exposure data input, dataframe")
  message("outcome_data = bulk outcome data input, dataframe")
  
  source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/harmonise_data_modified.R")
  library(tidyverse)
  
  id_exposure <- unique(exposure_data$id.exposure)
  id_outcome <- unique(outcome_data$id.outcome)

  exp_list <- lapply(id_exposure, function(i) {
    out_list <- lapply(id_outcome,function(j){
      tmp_exp <- exposure_data %>% dplyr::filter(id.exposure == i)
      tmp_out <- outcome_data %>% dplyr::filter(id.outcome == j)
      
      tmp_out <- harmonise_data_modifed(
        exposure_dat = tmp_exp,
        outcome_dat = tmp_out
      )
    
      if (nrow(tmp_out) > 0) {
        return(tmp_out)
      } else {
        return(NULL)
      }
    })
   each_exp <- do.call(rbind, out_list[sapply(out_list, Negate(is.null))])
   return(each_exp)
  })

  out <- do.call(rbind, exp_list[sapply(exp_list, Negate(is.null))])
  
  return(out)
}
