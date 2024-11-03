get_bulk_harmonised_data <- function(exposure_data = bulkexposure, 
                                      outcome_data = outcome
                                    ) {
  message("exposure_data = bulk exposure data input, dataframe")
  message("outcome_data = one specific outcome data, dataframe")
  
  source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/harmonise_data_modified.R")
  library(tidyverse)
  
  id_exposure <- unique(exposure_data$id.exposure)

  out_list <- lapply(id_exposure, function(i) {
    tmp <- exposure_data %>% dplyr::filter(id.exposure == i)
    
    tmp_out <- harmonise_data_modifed(
      exposure_dat = tmp,
      outcome_dat = outcome_data
    )
    
    if (nrow(tmp_out) > 0) {
      return(tmp_out)
    } else {
      return(NULL)
    }
  })

  out <- do.call(rbind, out_list[sapply(out_list, Negate(is.null))])
  
  return(out)
}
