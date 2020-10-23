model_MR = function(gene_set_name, treatment, IV, controls){
  extract_t = function(x, var) broom::tidy(x) %>% filter(str_detect(term, var))
  source("R/utils.R", local = TRUE)
  gene_set = pluck(signatures, "outcome_set", gene_set_name)
  
  # perform TSLS using 'ivreg()'
  library(AER)
  
  datt =
    prepro(gene_set, IV, controls) %>% 
    dplyr::rename(IV = IV) %>% 
    remove_diseased_subjects_from_datt(gene_set_name, controls)
  
  pca_rotated = fit_pca_util(datt, gene_set, "oblimin") 
  datt = dplyr::select(datt, -gene_set) 
  gene_set = colnames(pca_rotated$scores)
  workflow_reg = 
    function(outcome, datt){
      
      datt_pca = bind_cols(datt, !!outcome := pca_rotated$scores[, outcome])
      AER::ivreg(outcome ~ treatment + controls | )
      
      (rec =
          str_c(outcome) %>%
          str_c(" ~ .") %>%
          as.formula() %>%
          recipe(datt_pca))
      (mod = linear_reg() %>%
          set_engine("lm"))
      (wf =
          workflow() %>%
          add_recipe(rec) %>%
          add_model(mod) %>%
          fit(datt_pca))
    }
  
  
  
  gene_set %>% 
    set_names() %>% 
    map(workflow_reg, datt = datt) 
  
  cig_ivreg <- ivreg(log(packs) ~ log(rprice) | salestax, data = c1995)
  
  
  
  # for genetic(PGSBMI) associations with the outcomes(PCs) 
  
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

