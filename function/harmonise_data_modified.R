#This part of the code references the harmonise_data function of the TwosampleMR package and the Linkage Disequilibrium information information of 1000G(http://grch37.rest.ensembl.org/documentation/info/ld_id_get and http://rest.ensembl.org/documentation/info/ld_id_get), as well as the functions by Jianfeng Lin(https://github.com/linjf15/MR_tricks)
source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/snp_replace_proxy.R")
source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/get_proxy.R")
source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/get_SNP_eaf.R")
harmonise_data_modifed <- function(exposure_dat, outcome_dat, r2_thershold = 0.8, build = "37", pop = "EUR")
{
  stopifnot(build %in% c("37","38"))
  
  library(TwoSampleMR)

  outcome_from <- ifelse(is.character(outcome_dat), "ieugwas", "local")
  outcome_id <- ifelse(is.character(outcome_dat), outcome_dat, NA)
  snp_in_exposure_dat <- exposure_dat$SNP
  
  if(outcome_from == "ieugwas")
  {
    outcome_dat <- extract_outcome_data(snp_in_exposure_dat, outcome_id,
                                        proxies = TRUE, rsq = r2_thershold,
                                        palindromes = 1)
  }
  
  harmonised_dat <- harmonise_data(exposure_dat, outcome_dat)
  all_snp <- snp_in_exposure_dat[grepl("rs\\d+", snp_in_exposure_dat)]
  
  if (nrow(harmonised_dat) > 0 && any(harmonised_dat$mr_keep == "TRUE")) 
  {
    harmonised_dat <- harmonised_dat[harmonised_dat$mr_keep == "TRUE",]
    harmonised_dat$data_source.exposure <- "input"
    
    snp_available <- harmonised_dat$SNP[grepl("rs\\d+", harmonised_dat$SNP)]
    snp_na <- all_snp[!all_snp %in% snp_available]
  } else {
    snp_na <- all_snp
  }
  
  if (!is.null(snp_na)) {
    missed_snps_data <- lapply(snp_na, function(missed_snp) {
      missed_snp_proxy <- get_proxy(missed_snp, r2_thershold, build, pop)
      
      if (length(missed_snp_proxy) == 0){
        return(NULL)
      }else{
      time_search <- 0
      Sys.sleep(0.1)
      print(paste0("Proxies for ", missed_snp, ": ", paste0(missed_snp_proxy, collapse = ", ")))
      
      exposure_dat_i <- exposure_dat[exposure_dat$SNP == missed_snp,]
      
      proxy_data <- NULL  # 初始化为 NULL
      
      for (proxy_snp_i in missed_snp_proxy) {
        time_search <- time_search + 1
        if (time_search > length(missed_snp_proxy)) break
        
        if(outcome_from == "ieugwas") 
        {
          outcome_dat_i <- extract_outcome_data(proxy_snp_i, outcome_id, proxies = F)
          if(is.null(outcome_dat_i)) next
        }
        
        if(outcome_from == "local") 
        {
          if (!proxy_snp_i %in% outcome_dat$SNP) next
          
          outcome_dat_i <- outcome_dat[outcome_dat$SNP == proxy_snp_i,]
          if (is.na(outcome_dat_i[["eaf.outcome"]])){
          columns <- names(outcome_dat_i)
          outcome_dat_i <- snp_add_eaf(outcome_dat_i, build, pop)
          outcome_dat_i[["eaf.outcome"]] <- outcome_dat_i[["eaf"]]
          outcome_dat_i <- outcome_dat_i[,c(columns)]
          }
          outcome_dat_i <- snp_replace_proxy(outcome_dat_i, missed_snp, type = "outcome", build, pop)
        }
        
        if (is.null(outcome_dat_i)) next
        
        #if (is.na(outcome_dat_i$effect_allele.outcome) || is.na(outcome_dat_i$other_allele.outcome)) next
        if (any(is.na(outcome_dat_i$effect_allele.outcome)) || any(is.na(outcome_dat_i$other_allele.outcome))) next

        dat_supp_i <- harmonise_data(exposure_dat_i, outcome_dat_i)
        
        if (nrow(dat_supp_i) == 0) next
        else {
          if (dat_supp_i$mr_keep == "TRUE") {
            dat_supp_i$data_source.exposure <- proxy_snp_i
            proxy_data <- dat_supp_i  
            break
          }
        }
      }
      
      return(proxy_data)
      
      }
    })
    
    # Remove NULL elements
    missed_snps_data <- missed_snps_data[!sapply(missed_snps_data, is.null)]
    
    # Combine the list of data frames into a single data frame
    if (length(missed_snps_data) > 0) {
      harmonised_dat <- do.call(rbind, c(list(harmonised_dat), missed_snps_data))
    }
  }
  
  return(harmonised_dat)
}

