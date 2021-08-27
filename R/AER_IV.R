AER_IV = function(gene_set_name, treatment, IV, controls){
  # https://cran.r-project.org/web/packages/ivreg/vignettes/Diagnostics-for-2SLS-Regression.html
  # https://www.r-bloggers.com/2013/09/detecting-weak-instruments-in-r/
  # https://stats.stackexchange.com/questions/134789/interpretation-of-ivreg-diagnostics-in-r
  # https://bookdown.org/ccolonescu/RPoE4/random-regressors.html
  library(AER)
  
  # source("R/utils.R", local = TRUE)
  gene_set = pluck(signatures, "outcome_set", gene_set_name)
  
  # for genetic(PGSBMI) associations with the exposure(w5bmi)  
  datt =
    prepro(gene_set, treatment, c(controls,IV)) %>% 
    dplyr::rename(treatment = treatment) %>% 
    remove_diseased_subjects_from_datt(gene_set_name, controls)
  
  pca_rotated = fit_pca_util(datt, gene_set, rotate = "oblimin") 
  datt = dplyr::select(datt, -gene_set) 
  gene_set = colnames(pca_rotated$scores)
  
      datt_pca = bind_cols(datt, pca_rotated$scores[, gene_set])
      
      keep = datt_pca %>% complete.cases()
      datt_pca = datt_pca[keep, ]
      # remove invariant columns
      datt_pca = Filter(function(x) length(unique(x))!=1, datt_pca)

  
  en2  = datt_pca %>% colnames() %>% str_subset("color_byinterviewer")

  p1 = c(controls, "treatment") %>% str_c(collapse=" + ")
  p2 = c(controls %>% setdiff(en2), IV) %>% str_c(collapse=" + ")
  form = gene_set %>% map(~ str_c(.x," ~ ", p1 ,"|", p2))
  iv_fun = function(form, datt_pca){
    fm_skincolor =  AER::ivreg(form, data = datt_pca)
    summary_out = summary(fm_skincolor, vcov = sandwich, df = Inf, diagnostics = TRUE)
    temp = NULL
    temp$diagnostics = summary_out$diagnostics
    temp$coef = summary_out$coefficients %>% as.data.frame() %>% rownames_to_column("term") %>% filter(term== "treatment")
    temp
  }
 out = form %>% map(~ iv_fun(.x, datt_pca))
  return(out)
}
