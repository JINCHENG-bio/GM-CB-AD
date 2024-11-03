### All the SNPs in the output from the exposure data will be queried against the requested outcomes
### Locally or in remote database using API
### several functions are required for harmonizing data: TwoSampleMR package, get_proxy, snp_replace_proxy


get_mv_harmonise_data_modifed <- function (exposure_dat, outcome_dat, harmonise_strictness = 2) 
{
    stopifnot(all(c("SNP", "id.exposure", "exposure", "effect_allele.exposure", 
        "beta.exposure", "se.exposure", "pval.exposure") %in% 
        names(exposure_dat)))
    nexp <- length(unique(exposure_dat$id.exposure))
    stopifnot(nexp > 1)
    tab <- table(exposure_dat$SNP)
    keepsnp <- names(tab)[tab == nexp]
    exposure_dat <- subset(exposure_dat, SNP %in% keepsnp)
    exposure_mat <- reshape2::dcast(exposure_dat, SNP ~ id.exposure, 
        value.var = "beta.exposure")
        
    #dat <- harmonise_data(subset(exposure_dat, id.exposure == 
    #    exposure_dat$id.exposure[1]), outcome_dat, action = harmonise_strictness)
    source("/mnt/data/lijincheng/mGWAS/result/02MRBMA/MRBMA_function/function/harmonise_data_modified.R")
    dat <-  harmonise_data_modifed(subset(exposure_dat, id.exposure == 
        exposure_dat$id.exposure[1]), outcome_dat, r2_thershold = 0.8, build = "38", pop = "EUR") ### consider SNP proxy for local data
        
    dat <- subset(dat, mr_keep)
    dat$SNP <- as.character(dat$SNP)
    exposure_beta <- reshape2::dcast(exposure_dat, SNP ~ id.exposure, 
        value.var = "beta.exposure")
    exposure_beta <- subset(exposure_beta, SNP %in% dat$SNP)
    exposure_beta$SNP <- as.character(exposure_beta$SNP)
    exposure_pval <- reshape2::dcast(exposure_dat, SNP ~ id.exposure, 
        value.var = "pval.exposure")
    exposure_pval <- subset(exposure_pval, SNP %in% dat$SNP)
    exposure_pval$SNP <- as.character(exposure_pval$SNP)
    exposure_se <- reshape2::dcast(exposure_dat, SNP ~ id.exposure, 
        value.var = "se.exposure")
    exposure_se <- subset(exposure_se, SNP %in% dat$SNP)
    exposure_se$SNP <- as.character(exposure_se$SNP)
    index <- match(exposure_beta$SNP, dat$SNP)
    dat <- dat[index, ]
    stopifnot(all(dat$SNP == exposure_beta$SNP))
    exposure_beta <- as.matrix(exposure_beta[, -1])
    exposure_pval <- as.matrix(exposure_pval[, -1])
    exposure_se <- as.matrix(exposure_se[, -1])
    rownames(exposure_beta) <- dat$SNP
    rownames(exposure_pval) <- dat$SNP
    rownames(exposure_se) <- dat$SNP
    outcome_beta <- dat$beta.outcome
    outcome_se <- dat$se.outcome
    outcome_pval <- dat$pval.outcome
    expname <- subset(exposure_dat, !duplicated(id.exposure), 
        select = c(id.exposure, exposure))
    outname <- subset(outcome_dat, !duplicated(id.outcome), select = c(id.outcome, 
        outcome))
    return(list(exposure_beta = exposure_beta, exposure_pval = exposure_pval, 
        exposure_se = exposure_se, outcome_beta = outcome_beta, 
        outcome_pval = outcome_pval, outcome_se = outcome_se, 
        expname = expname, outname = outname))
}
