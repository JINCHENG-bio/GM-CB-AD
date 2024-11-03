#  function
get_mv_exp_data <- function(l_full=l_full,
                     min_pval = 1e-200,
                     log_pval = FALSE, 
                     pval_threshold = 5e-08, 
                     clump_r2 = 0.001, 
                     clump_kb = 10000, 
                     harmonise_strictness = 2 )
{
  source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/local_clumb.R")
  l_inst <- list()
  for (i in 1:length(l_full)) {
    l_inst[[i]] <- subset(l_full[[i]], pval.outcome < pval_threshold)
    l_inst[[i]] <- convert_outcome_to_exposure(l_inst[[i]])
    l_inst[[i]] <- subset(l_inst[[i]], pval.exposure < pval_threshold)
    l_inst[[i]] <- local_clumb(l_inst[[i]], clump_p1 = pval_threshold, 
                              clump_r2 = clump_r2, clump_kb = clump_kb, pop = "EUR")
  }
  exposure_dat <- dplyr::bind_rows(l_inst)
  id_exposure <- unique(exposure_dat$id.exposure)
  temp <- exposure_dat
  temp$id.exposure <- 1
  temp <- temp[order(temp$pval.exposure, decreasing = FALSE), ]
  temp <- subset(temp, !duplicated(SNP))
  temp <- local_clumb(temp, clump_p1 = pval_threshold, 
                              clump_r2 = clump_r2, clump_kb = clump_kb, pop = "EUR")
  exposure_dat <- subset(exposure_dat, SNP %in% temp$SNP)
  d1 <- lapply(l_full, function(x) {
    subset(x, SNP %in% exposure_dat$SNP)
  }) %>% dplyr::bind_rows()
  stopifnot(length(unique(d1$id)) == length(unique(id_exposure)))
  d1 <- subset(d1, mr_keep.outcome)
  d2 <- subset(d1, id.outcome != id_exposure[1])
  d1 <- convert_outcome_to_exposure(subset(d1, id.outcome == 
                                             id_exposure[1]))
  #d <- harmonise_data(d1, d2, action = harmonise_strictness)
  source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/harmonise_data_modified.R") ##consider proxy
  d <- harmonise_data_modifed(exposure_dat=d1, outcome_dat=d2, r2_thershold = 0.8, build = "38", pop = "EUR")
  tab <- table(d$SNP)
  keepsnps <- names(tab)[tab == length(id_exposure) - 1] ### Only keep SNPs that are present in all
  d <- subset(d, SNP %in% keepsnps)
  dh1 <- subset(d, id.outcome == id.outcome[1], select = c(SNP, 
                                                           exposure, id.exposure, effect_allele.exposure, other_allele.exposure, 
                                                           eaf.exposure, beta.exposure, se.exposure, pval.exposure))
  dh2 <- subset(d, select = c(SNP, outcome, id.outcome, effect_allele.outcome, 
                              other_allele.outcome, eaf.outcome, beta.outcome, se.outcome, 
                              pval.outcome))
  names(dh2) <- gsub("outcome", "exposure", names(dh2))
  dh <- rbind(dh1, dh2)
  return(dh)
}