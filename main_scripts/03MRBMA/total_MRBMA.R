library(tidyverse)
library(data.table)

source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/get_each_MVinput.R")

#--------------------{1.getMVinput}-----------------------
##-------------{01 MiBioGen}---------------------------
get_MVinput_for_MRBMA(exposure="mibio",
                      outcome="LOAD",
                      outpath="/mnt/data/lijincheng/mGWAS/result/02MRBMA/foreachoutcome/")

get_MVinput_for_MRBMA(exposure="mibio",
                      outcome="abeta42",
                      outpath="/mnt/data/lijincheng/mGWAS/result/02MRBMA/foreachoutcome/")

get_MVinput_for_MRBMA(exposure="mibio",
                      outcome="ptau",
                      outpath="/mnt/data/lijincheng/mGWAS/result/02MRBMA/foreachoutcome/")

##-------------{02 FR02}---------------------------
get_MVinput_for_MRBMA(exposure="FR02",
                      outcome="LOAD",
                      outpath="/mnt/data/lijincheng/mGWAS/result/02MRBMA/foreachoutcome/")


get_MVinput_for_MRBMA(exposure="FR02",
                      outcome="ptau",
                      outpath="/mnt/data/lijincheng/mGWAS/result/02MRBMA/foreachoutcome/")


##-------------{03 pathways}---------------------------
get_MVinput_for_MRBMA(exposure="pathways",
                      outcome="LOAD",
                      outpath="/mnt/data/lijincheng/mGWAS/result/02MRBMA/foreachoutcome/")

get_MVinput_for_MRBMA(exposure="pathways",
                      outcome="abeta42",
                      outpath="/mnt/data/lijincheng/mGWAS/result/02MRBMA/foreachoutcome/")

get_MVinput_for_MRBMA(exposure="pathways",
                      outcome="ptau",
                      outpath="/mnt/data/lijincheng/mGWAS/result/02MRBMA/foreachoutcome/")


##-------------{04 metabolites}---------------------------
get_MVinput_for_MRBMA(exposure="metabolites",
                      outcome="LOAD",
                      outpath="/mnt/data/lijincheng/mGWAS/result/02MRBMA/foreachoutcome/")


get_MVinput_for_MRBMA(exposure="metabolites",
                      outcome="ptau",
                      outpath="/mnt/data/lijincheng/mGWAS/result/02MRBMA/foreachoutcome/")

#-------------------{2.MRBMA}----------------------
library(mrbma)
library(foreach)
totalinput <- list.files("/mnt/data/lijincheng/mGWAS/result/02MRBMA/foreachoutcome/",pattern="_MVinput.RDS$",full.names = T)
fileexist <- list.files("/mnt/data/lijincheng/mGWAS/result/02MRBMA/foreachoutcome/MRBMA_rawres/",pattern="_BMA.RDS$",full.names = F)
fileexist <- gsub("_BMA.RDS$", "", basename(fileexist))
total_basename <- gsub("_MVinput.RDS$", "", basename(totalinput))
totalref <- data.frame(totalinput=totalinput,total_basename=total_basename)
finalinput <- totalref %>% dplyr::filter(! total_basename %in% fileexist)

library(foreach)
library(doParallel)

# 设置并行集群
cl <- makeCluster(3)
registerDoParallel(cl)
# 修改foreach循环
foreach(i = 1:dim(finalinput), .packages=c("tidyverse", "mrbma", "openxlsx")) %dopar% {
    input <- finalinput$totalinput[i]
    name <- finalinput$total_basename[i]
    MVinput <- readRDS(input)
    res_bma <- mr_bma(
        MVinput,
        prior_prob = 0.1,
        prior_sigma = 0.5,
        top = 10,
        remove_outliers = TRUE,
        remove_influential = TRUE,
        calculate_p = TRUE,
        nrepeat = 10000
    )
    saveRDS(res_bma, file = paste0("/mnt/data/lijincheng/mGWAS/result/02MRBMA/foreachoutcome/MRBMA_rawres/", name, "_BMA.RDS"))
    openxlsx::write.xlsx(x = list("best_model" = res_bma$model_best, "rank_factor" = res_bma$mip_table),
                         file = paste0("/mnt/data/lijincheng/mGWAS/result/02MRBMA/foreachoutcome/MRBMA_rawres/", name, "_BMA.xlsx"))
}

# 停止并行集群
stopCluster(cl)

