load_data = function(reconciled, remove_inflam){ 
  
  # COMPARE WX AND JC RESULTS
  if(reconciled){
    # ... after reconciling
    # signatures = readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_21042020_signature.rds") # from Wenjia
    # dat = readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_21042020.rds")
    signatures <- readRDS("/home/share/preprocessed_two_batches/recon_25.08.2020/dt_batches1_2_steve_waves_25.08.2020_signature.rds")
    dat <- readRDS("/home/share/preprocessed_two_batches/recon_25.08.2020/dt_batches1_2_steve_waves_25.08.2020.rds")
    
  } else {  
    # the analyses *before* ensuring all the differently normalized datasets have identical genes
    # dat = readRDS("/home/share/preprocessed_two_batches/wx/dt_batches1_2_steve_waves_21042020.rds")
    # signatures = readRDS("/home/share/preprocessed_two_batches/for_cecilia/dt_batches1_2_steve_waves_21042020_signature.rds") # from Wenjia
    # after adding two new signatures and several pheno data
    signatures <- readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_22.03.2021_signature.rds")
    dat <- readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_22.03.2021.rds")
  }
  
  if(remove_inflam) {
    # remove inflamation genes in each signatures
    signatures$outcome_set = map(signatures$outcome_set, ~setdiff(.x, signatures$outcome_set$inflam1k_mRNA))
  }
  
  signature_names = signatures$outcome_set %>% names(.)
  
  # POSSIBLY?
  signatures$outcome_set = signatures$outcome_set %>% map(str_replace_all, "-", "_") # otherwise problems downstream
  featureNames(dat) <-  featureNames(dat) %>% str_replace_all("-", "_")
  
  # Loose funny names, very coarse hack
  funny_names = signatures$outcome_set %>% map(str_subset, "-")
  signatures$outcome_set = list(signatures$outcome_set, funny_names) %>% pmap(setdiff)
  
  list(dat = dat,
       signatures = signatures,
       signature_names = signature_names) %>% 
    list2env(.GlobalEnv)
  
  if(0) {
    
    library(sas7bdat)
    fromSAS = read.sas7bdat("/data/addhealth/addhealthdata/wave5/biomarkers/w5biocovars.sas7bdat")
    
    library(SASxport)
    crp = read.xport("/data/addhealth/addhealthdata/waves_1_4/W4 Supplemental Files-Biomarker/Measures of EBV and hsCRP/crp_ebv.xpt")
    crp = crp %>% mutate(high_crp = as.numeric((CRP > 3))) %>% 
      select(high_crp, AID) %>%
      mutate(high_crp = na_if(high_crp, 998), 
             high_crp = na_if(high_crp, 999))
    
    # For CRP, we used the established clinical cutoff (3 mg/L [to convert to
    # nanomoles per liter, multiply by 9.524]) to identify participants with high
    # CRP levels; thus, high CRP level indicates a CRP level greater than 3 mg/L
    # (n = 287 [20.6%]).
    # https://jamanetwork.com/journals/jamapediatrics/fullarticle/2754102
    
  }
}