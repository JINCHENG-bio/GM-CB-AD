library(tidyverse)
library(hyprcoloc)

total_result <- data.frame()

inputpath <- "/mnt/data/lijincheng/mGWAS/result/02MRBMA/hyprcoloc/MVinput_res/"
inputfiles <- list.files(inputpath,pattern="_MVinput.RDS$",full.names = T)
full_name <- gsub("_MVinput.RDS$", "", basename(inputfiles))

for(i in 1:length(inputfiles)){
  input <-  readRDS(inputfiles[i])
  
  input_beta <- input$exposure_beta %>% as.data.frame()
  input_se <- input$exposure_se %>% as.data.frame()
  
  if(grepl("LOAD",full_name[i])){
    input_beta$LOAD <- input$outcome_beta
    input_se$LOAD <- input$outcome_se
    binary.traits = c(0,0,1)
    
  }else if(grepl("ADproxy",full_name[i])){
    input_beta$ADproxy <- input$outcome_beta
    input_se$ADproxy <- input$outcome_se
    binary.traits = c(0,0,1)
    
  }else if(grepl("abeta42",full_name[i])){
    input_beta$abeta42 <- input$outcome_beta
    input_se$abeta42 <- input$outcome_se
    
    binary.traits = c(0,0,0)
    
  }else if(grepl("ptau",full_name[i])){
    input_beta$ptau <- input$outcome_beta
    input_se$ptau <- input$outcome_se
    binary.traits = c(0,0,0)
  }
  
  input_beta <- as.matrix(input_beta)
  input_se <- as.matrix(input_se)
  traits <- colnames(input_beta)
  rsid <- rownames(input_beta);

  res <- hyprcoloc(input_beta, input_se, trait.names=traits, snp.id=rsid, binary.outcomes = binary.traits);

  result <- res$result %>% as.data.frame()  
  total_result <- bind_rows(total_result,result)
}

openxlsx::write.xlsx(x=list("total_result13"=total_result),
                    file="/mnt/data/lijincheng/mGWAS/result/02MRBMA/hyprcoloc/hyprcoloc_result13.xlsx")


