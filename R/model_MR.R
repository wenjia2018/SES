model_MR = function(gene_set_name, treatment, IV, controls){
  extract_t = function(x, var) broom::tidy(x) %>% filter(str_detect(term, var))
  source("R/utils.R", local = TRUE)
  gene_set = pluck(signatures, "outcome_set", gene_set_name)
  
  # for genetic(PGSBMI) associations with the exposure(w5bmi)  
  data = pData(dat) %>% dplyr::select(all_of(controls), all_of(treatment), IV) %>% 
    dplyr::rename(treatment = treatment, IV = IV) 
  outbx = lm(treatment ~ ., data = data) %>% extract_t("IV")
    
 # for genetic(PGSBMI) associations with the outcomes(PCs) 
  datt =
    prepro(gene_set, IV, controls) %>% 
    dplyr::rename(IV = IV) %>% 
    remove_diseased_subjects_from_datt(gene_set_name, controls)
  
    outby = fit_m7(datt, gene_set, "oblimin")$fit %>% map(extract_t, "IV") 
  
 MR = function(outbx, outby){
   MRInputObject <- mr_input(bx = outbx$estimate,
                             bxse = outbx$std.error,
                             by = outby$estimate,
                             byse = outby$std.error)
   IVWObject <- mr_ivw(MRInputObject,
                       model = "default",
                       robust = FALSE,
                       penalized = FALSE,
                       correl = FALSE,
                       weights = "simple",
                       psi = 0,
                       distribution = "normal",
                       alpha = 0.05)
   
   return(IVWObject) 
 } 
out = outby %>% map(~ MR(outbx, .x))
}