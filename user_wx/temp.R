mediate_pca = function(gene_set, pca){
  datt_m = 
    prepro(gene_set, treatment, c(controls, mediator)) %>%
    rename(treatment = treatment,
           mediator = mediator) %>% 
    remove_diseased_subjects_from_datt(gene_set_name, controls)
  
  ok = complete.cases(datt_m$mediator)
  datt_m = datt_m[ok, ]
  fit_m99(datt_m, gene_set) %>% extract_m99()
  
  
  outcome = colnames(pca$scores)
  args_m = 
    crossing(mediators, outcome) %>% 
    rename(mediator = mediators) %>% 
    mutate(names = str_c(mediator, outcome, sep = "_"))
  map2(args_m$mediator, args_m$outcome, safely(mediate)) %>% 
    set_names(args_m$names)
  
  
}

mediate = function(mediator, gene_set){
  
  datt_m = 
    prepro(gene_set, treatment, c(controls, mediator)) %>%
    rename(treatment = treatment,
           mediator = mediator) %>% 
    remove_diseased_subjects_from_datt(gene_set_name, controls)
  
  ok = complete.cases(datt_m$mediator)
  datt_m = datt_m[ok, ]
  fit_m99(datt_m, gene_set) %>% extract_m99()  
}

mediate_multiple = function(gene_set){ 
  
  args_m97 = 
    crossing(mediators, gene_set) %>% 
    rename(mediator = mediators) %>% 
    mutate(names = str_c(mediator, gene_set, sep = "_"))
  map2(args_m97$mediator, args_m97$gene_set, safely(mediate)) %>% 
    set_names(args_m97$names)
  
}