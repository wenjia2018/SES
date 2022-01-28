load_data = function(reconciled, remove_inflam){ 
  # raw and outcome_set_full is for DE analysis from raw and depends on design matrix
  raw = readRDS("/home/share/preprocessing/from_Brandt/all.batches.expression.set.070121.Rds")
  # as the resulting signature maybe different for different treatment control, therefore creating this all possible signatures
  # which includes all the genes from raw, only for m12
  outcome_set_full = readRDS("/home/share/preprocessing/preprocessed_two_batches/allpossiblegene.rds")
 
    signatures <- readRDS("/home/share/preprocessing/preprocessed_two_batches/all.batches.expression.set.tmm_waves_01.09.2021_signature.rds")
   # ask Ravi which file is the newly updated one with "flagged" subjects removed
     dat <- readRDS("/home/share/preprocessing/preprocessed_two_batches/all.batches.expression.set.tmm_waves_01.09.2021.rds")
  
  
  if(remove_inflam) {
    # remove inflamation genes in each signatures
    outcome_set_full$outcome_set = map(outcome_set_full$outcome_set, ~setdiff(.x, outcome_set_full$outcome_set$inflam1k_mRNA))
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
       signature_names = signature_names,
       raw = raw,
       outcome_set_full = outcome_set_full) %>% 
    list2env(.GlobalEnv)
}