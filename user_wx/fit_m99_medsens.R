
fit_m99 = function(datt, gene_set, out = NULL){
  
  datt = bind_cols(dplyr::select(datt, -gene_set),
                   gene_set = apply(datt[gene_set], 1, mean))
  gene_set = "gene_set"
  
  if (is.numeric(datt$mediator)){
    
    # gene_set outcome is dependent measure
    (rec =
       str_c(gene_set) %>%
       str_c(" ~ .") %>%
       as.formula() %>%
       recipe(datt))
    (mod = linear_reg() %>%
        set_engine("lm"))
    (wf =
        workflow() %>%
        add_recipe(rec) %>%
        add_model(mod) %>%
        fit(datt))
    out.y = pluck(wf, "fit", "fit", "fit")
    
    # mediate is dependent measure    
    (rec =
        str_c("mediator") %>%
        str_c(" ~ .") %>%
        as.formula() %>%
        recipe(datt %>% dplyr::select(-gene_set)))
    
    (mod = linear_reg() %>%
        set_engine("lm"))
    
    (wf =
        workflow() %>%
        add_recipe(rec) %>%
        add_model(mod) %>%
        fit(datt %>% dplyr::select(-gene_set)))
    
    out.m = pluck(wf, "fit", "fit", "fit")
    
  } else if(datt$mediator %>% table() %>% length() == 2){
    
    datt$mediator = datt$mediator %>% as.character() %>%  as.numeric()
    keep = datt %>% complete.cases()
    datt_keep = datt[keep,]
    
    formula_y = str_c("gene_set", " ~ .") %>% as.formula()
    out.y = lm(formula_y, data = datt_keep)
    
    formula_m = str_c("mediator", " ~ .") %>% as.formula()
    out.m = glm(formula_m, family = binomial(link = "probit"), data = datt_keep %>% dplyr::select(-gene_set))
    
  } else if(datt$mediator %>% table() %>% length() > 2){
    
    datt$mediator = datt$mediator %>% as.factor()
    keep = datt %>% complete.cases()
    datt_keep = datt[keep,]
    
    formula_y = str_c("gene_set", " ~ .") %>% as.formula()
    out.y = lm(formula_y, data = datt_keep)
    # has met some numerical problem with the starting value for optimization
    # therefore not stable, might have again for other mediators outcome combination.
    # https://stackoverflow.com/questions/28916377/r-error-with-polr-initial-value-in-vmmin-is-not-finite
    formula_m = str_c("mediator", " ~ .") %>% as.formula()
    out.m = MASS::polr(formula_m, data = datt_keep %>% dplyr::select(-gene_set), Hess=TRUE)
    
  }
  
  
  med.out <- mediation::mediate(out.m, out.y, treat = "treatment", mediator = "mediator")
  med.sens <- mediation::medsens(med.out, rho.by = 0.1, sims = 1000,
                                 eps = 0.01, effect.type = "indirect")
  # out = NULL
  out$med = med.out
  out$sens = med.sens
  out
}


extract_m99 = function(m, out = NULL){
  
  extract_p = function(x) {
    m$med %>% 
      summary %>% 
      capture.output() %>% 
      str_subset(x) %>% 
      str_split(" ", simplify = TRUE) %>%
      str_subset("[0-9]") %>% 
      pluck(4)
    
  }
  out$detail = "too big" #m
  out$p = extract_p("ACME")
  out$other$sens = m$sens
  return(out = out)
}

