library(data.table)
library(tidyverse)

#source("/mnt/data/lijincheng/mGWAS/data/outcome/Jansen_Aβ_tau/get_rsid_from_chro_pos.R")


##https://www.ebi.ac.uk/gwas/publications/36066633
###   GCST90129599_buildGRCh38.tsv     Cerebrospinal fluid amyloid beta 42 levels
###   GCST90129600_buildGRCh38.tsv     Cerebrospinal fluid p-tau levels

amyloid_beta42 <- fread("/mnt/data/lijincheng/mGWAS/data/outcome/Jansen_Aβ_tau/GCST90129599_buildGRCh38.tsv") %>% mutate(trait = "Cerebrospinal fluid amyloid beta 42 levels")

p_tau <- fread("/mnt/data/lijincheng/mGWAS/data/outcome/Jansen_Aβ_tau/GCST90129600_buildGRCh38.tsv") %>% mutate(trait = "Cerebrospinal fluid p-tau levels")

###由于需要查找的SNP较多，因此直接与hg38匹配

library(BSgenome)
BSgenome::available.SNPs()

##BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")

library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38

##amyloid_beta42
for (i in unique(amyloid_beta42$chromosome)){
  my_pos <- amyloid_beta42$base_pair_location[amyloid_beta42$chromosome == i]
  chr_snps <- snpsBySeqname(snps, as.character(i))
  idx <- match(my_pos,pos(chr_snps))
  rsids <- mcols(chr_snps)$RefSNP_id[idx]
  amyloid_beta42$rsid[amyloid_beta42$chromosome == i] <- rsids
  print(paste(as.character(i),"is ok"))
}

##p_tau
for (i in unique(p_tau$chromosome)){
  my_pos <- p_tau$base_pair_location[p_tau$chromosome == i]
  chr_snps <- snpsBySeqname(snps, as.character(i))
  idx <- match(my_pos,pos(chr_snps))
  rsids <- mcols(chr_snps)$RefSNP_id[idx]
  p_tau$rsid[p_tau$chromosome == i] <- rsids
  print(paste(as.character(i),"is ok"))
}

saveRDS(amyloid_beta42, file = "/mnt/data/lijincheng/mGWAS/data/outcome/Jansen_Aβ_tau/amyloid_beta42.rds")
saveRDS(p_tau, file = "/mnt/data/lijincheng/mGWAS/data/outcome/Jansen_Aβ_tau/p_tau.rds")


final_amyloid_beta42 <- readRDS("/mnt/data/lijincheng/mGWAS/data/outcome/Jansen_Aβ_tau/amyloid_beta42.rds")
final_p_tau <- readRDS("/mnt/data/lijincheng/mGWAS/data/outcome/Jansen_Aβ_tau/p_tau.rds")

###添加EAF，计算beta，se
source("/mnt/data/lijincheng/mGWAS/data/outcome/Jansen_Aβ_tau/get_MAF_from_EAF.R")

##  BETA: $Z/sqrt((2*$MAF*(1-$MAF))*($N+ ($Z)^2))
##  SE: 1/sqrt((2*$MAF*(1-$MAF))*($N+($Z)^2))

final_amyloid_beta42 <- final_amyloid_beta42 %>%
                        dplyr::mutate(MAF = eaf2maf(effect_allele_frequency),
                                      BETA = Zscore/sqrt((2*MAF*(1-MAF))*(N+ (Zscore)^2))   ,
                                      SE = 1/sqrt((2*MAF*(1-MAF))*(N+(Zscore)^2))  ,
                                      )

final_p_tau <- final_p_tau %>%
                dplyr::mutate(MAF = eaf2maf(effect_allele_frequency),
                              BETA = Zscore/sqrt((2*MAF*(1-MAF))*(N+ (Zscore)^2))   ,
                              SE = 1/sqrt((2*MAF*(1-MAF))*(N+(Zscore)^2))  ,
                              )

saveRDS(final_amyloid_beta42, file = "/mnt/data/lijincheng/mGWAS/data/outcome/Jansen_Aβ_tau/final_amyloid_beta42.rds")
saveRDS(final_p_tau, file = "/mnt/data/lijincheng/mGWAS/data/outcome/Jansen_Aβ_tau/final_p_tau.rds")






























