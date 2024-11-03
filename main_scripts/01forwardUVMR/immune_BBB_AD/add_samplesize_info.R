add_samplesize_info <- function(dat=dat,samplesize=samplesize)
  {
  load("/mnt/data/lijincheng/mGWAS/result/01UVMR_immune_BBB/immune_BBB_info.Rdata")
  message("dat= bulk harmonised data")
  message("samplesize= 'samplesize.exposure' or 'samplesize.outcome' ")
  
  library(dplyr)
  library(tidyverse)
  
  BBB_immune_info <- BBB_immune_info %>% dplyr::select(c("id","sample_size"))
  col_names <- names(dat)
  
  if(samplesize == "samplesize.exposure")
  {
    missing_sam_rows <- dat %>% filter(is.na(samplesize.exposure))

    # 根据 B 中的对应 id 列的值填充缺失的 sam 列
    filled_sam <- missing_sam_rows %>%
    left_join(.,BBB_immune_info, by = c("id.exposure" = "id")) %>%
    mutate(samplesize.exposure = coalesce(samplesize.exposure, sample_size)) %>%
    select(col_names)
  
  # 
    result <- bind_rows(dat %>% filter(!is.na(samplesize.exposure)), filled_sam)
  }
  else if(samplesize == "samplesize.outcome")
  {
    missing_sam_rows <- dat %>% filter(is.na(samplesize.outcome))

    # 根据 B 中的对应 id 列的值填充缺失的 sam 列
    filled_sam <- missing_sam_rows %>%
    left_join(.,BBB_immune_info, by = c("id.exposure" = "id")) %>%
    mutate(samplesize.exposure = coalesce(samplesize.outcome, sample_size)) %>%
    select(col_names)
  
  # 
    result <- bind_rows(dat %>% filter(!is.na(samplesize.outcome)), filled_sam)
    }
  
  return(result)
  }