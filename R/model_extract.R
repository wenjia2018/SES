
############################################################
# EACH FUNCTION SHOULD RETURN THREE VALUES detail, p, other
############################################################

extract_m1 = function(m, out = NULL){
  
  out$detail =
    capture.output(car::Manova(m))  %>% 
    str_subset("Df|treatment")
  
  nmbrs = capture.output(car::Manova(m))  %>% 
    str_subset("treatment") %>% 
    str_split(" ", simplify = TRUE) %>%
    str_subset("[0-9]")
  
  print(str_c("ok = ", length(nmbrs)  == 6))
  
  out$p = nmbrs %>%
    pluck(6) %>% 
    as.numeric()
  
  out$other$Pillai = nmbrs %>% 
    pluck(2) %>% 
    as.numeric()
  
  return(out = out)
  
}

extract_m2 = function(m, out = NULL){
  
  extract_anova = function(x) anova(x) %>% tidy %>% filter(str_detect(term, "treatment"))
  
  # Univariate parametric
  out$detail = extract_anova(m)
  out$p = extract_anova(m)$p.value
  out$other = "empty"
  
  return(out = out)
}

extract_m3 = extract_m2

extract_m4 = function(m, out = NULL){
  
  out$detail = "permutation test"
  out$p = mean(m$permutation_statistics$teststat > m$unpermuted_teststat)
  out$out = "empty"
  
  return(out = out)
}

extract_m5 = function(m, out = NULL){
  
  out$detail = "too big" # m
  out$p = min(m$adj.P.Val)
  out$other = "empty"
  
  return(out = out)
}

extract_m5b = function(m, out = NULL){
  out$detail = "too big" # m
  # out$p = any(m$P.Value < (0.05 / length(m$P.Value))) # does any P.value exceed Bonferonni threshold?
  out$p = min(m$P.Value * length(m$P.Value)) # smallest FWER corrected pvalue
  out$other = "empty"
  
  return(out = out)
}

extract_m6 = extract_m1

extract_m7 =  function(m, out = NULL){
  
  extract_anova = function(x) anova(x) %>% tidy %>% filter(str_detect(term, "treatment"))
  extract_t = function(x) broom::tidy(x) %>% filter(str_detect(term, "treatment"))
  
  # Univariate parametric
  out$detail$anova = m$fit %>% map(extract_anova)
  out$detail$t = m$fit %>% map(extract_t)
  out$p = out$detail$anova %>% map_dbl(pluck("p.value"))
  out$other$varexplained = m$varexplained 
  out$other$loadings = m$loadings
  if(exists("ftest_v")){
  extract_hyp = function(x) car::linearHypothesis(x, str_subset(names(coef(x)), ftest_v), verbose=F)
  out$all = m$fit %>% map(extract_hyp)
  }
  return(out = out)
}


extract_m8_fwer = function(m, out = NULL){
  
  out$detail = NA # m
  out$p = min( min(m$P.Value * length(m$P.Value)), 1) # smallest FWER corrected pvalue
  out$other = "empty"
  
  return(out = out)
}

extract_m8_fdr = function(m, out = NULL){
  
  out$detail = m %>% slice(which.min(adj.P.Val)) # m
  out$p = min(m$adj.P.Val) # smallest FWER corrected pvalue
  out$other = list(m = m, avg_logFC = mean(m$logFC))
  out$sig_genes = m %>% filter(adj.P.Val < 0.05) %>% pull(gene)

  return(out = out)
}


extract_m98 = extract_m7

extract_m99 = function(m, out = NULL){
  
  extract_med = function(x, y) {
    m %>% 
      summary %>% 
      capture.output() %>% 
      str_subset(x) %>% 
      str_split(" ", simplify = TRUE) %>%
      str_subset("[0-9]") %>% 
      pluck(y)
    
  }
  
  
  out$detail = "too big" #m
  out$p = extract_med("ACME", 4)
  out$other$med_prop = extract_med("Prop. Mediated", 1)
  out$other$med_ACME = extract_med("ACME", 1)
  out$other$med_ADE = extract_med("ADE", 1)
  out$other$med_ADE_p = extract_med("ADE", 4)
  return(out = out)
}
