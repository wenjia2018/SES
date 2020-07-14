
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
  
  out$detail = m
  out$p = min(m$adj.P.Val)
  out$other = "empty"
  
  return(out = out)
}

extract_m6 = extract_m1

extract_m7 =  function(m, out = NULL){
  
  extract_anova = function(x) anova(x) %>% tidy %>% filter(str_detect(term, "treatment"))
  
  # Univariate parametric
  out$detail = m %>% map(extract_anova)
  out$p = out$detail %>% map_dbl(pluck("p.value"))
  out$other = "empty"
  return(out = out)
}

extract_m98 = extract_m7

extract_m99 = function(m, out = NULL){
  
  extract_p = function(x) {
    m %>% 
      summary %>% 
      capture.output() %>% 
      str_subset(x) %>% 
      str_split(" ", simplify = TRUE) %>%
      str_subset("[0-9]") %>% 
      pluck(4)
    
  }
  out$detail = "see code details" # m # m is quite big
  out$p = extract_p("ACME")
  out$other = "empty"
  return(out = out)
}
