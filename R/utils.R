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

mediate_pca = function(mediator, gene_set, rotate, pca_out){
  datt_m = 
    prepro(gene_set, treatment, c(controls, mediator)) %>%
    rename(treatment = treatment
           # mediator = mediator
           ) %>% 
    remove_diseased_subjects_from_datt(gene_set_name, controls)
  
  pca_rotated = fit_pca_util(datt_m, gene_set, rotate) 
  
  datt_m = dplyr::select(datt_m, -gene_set) 
  outcome = colnames(pca_rotated$scores) %>% set_names()
  
  convert_outcome = function(outcome, datt_m){
    datt_pca = bind_cols(datt_m, !!outcome := pca_rotated$scores[, outcome])
  }
  
  datt_pca = outcome %>% map(convert_outcome, datt = datt_m)
  # only do mediation for significant PCs, p values are corrected by bonferonni
  # actually should divided by 4 the treatment numbers, but Mike sometimes keep some
  # signatures when their p value is close and a little bigger to the cut off, so here
  # choose correcting pca component numbers only so mediational analysis will be done to less cases, save time.
  threshold = 0.05 / length(pca_out$p)
  sig = pca_out$p < threshold
  print("please be aware: p value correction is bonferonni for now")
 # map2(datt_pca[sig], outcome[sig], fit_m99) %>% map(extract_m99)
 pmap(list(datt_pca[sig], outcome[sig], mediator), fit_m99) %>% map(extract_m99)
}

get_table1 = function(example){ 
  
  tab1a = 
    example %>% 
    hoist(result, !!!funcs)  %>% 
    unnest(!!funcs) %>% 
    # hoist(m1, pm1 = "p") %>% 
    # hoist(m2, pm2 = "p") %>% 
    # hoist(m3, pm3 = "p") %>% 
    # hoist(m5, pm5 = "p") %>%
    # hoist(m5b, pm5b = "p") %>%
    hoist(m8_fwer, m8_fwer_p = "p") %>% 
    hoist(m8_fdr, m8_fdr_p = "p") %>% 
    discard(is.list)
  
  if(0){ 
    
    # mediation
    tab1b = 
      example %>% 
      hoist(result, !!!funcs)  %>% 
      unnest(m99) %>%
      unnest_wider(m99) %>% 
      hoist(w5bmi, w5bmi_p = c("result", "p"))  %>% 
      hoist(bingedrink, bingedrink_p = c("result", "p"))  %>% 
      hoist(currentsmoke, currentsmoke_p = c("result", "p"))  %>% 
      hoist(phys_activ_ff5, phys_activ_ff5_p = c("result", "p"))  %>%   
      discard(is.list)
    
    tab1a %>% left_join(tab1b)
    
  } else {
    
    tab1a
    
  }
}

detach_package <- function(abcpkg, character.only = FALSE)
  # https://stackoverflow.com/questions/6979917/how-to-unload-a-package-without-restarting-r
{
  if(!character.only)
  {
    abcpkg <- deparse(substitute(abcpkg))
  }
  search_item <- paste("package", abcpkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}


remove_errors = function(example){
  
  # ERRORS:
  example %>%
    hoist(out, "error") %>%
    mutate(error = map(error, as.character)) %>%
    unnest(error) %>%
    group_by(error) %>%
    slice(1) %>% 
    print()
  
  # WHAT CAUSES ERROR? RELATE NA TO ARGS OF model_fit()
  example %>%
    hoist(out, p = list("result", "m5", 1, "p")) %>%
    with(table(gene_set_name, is.na(p))) 
  
  # REMOVE MODELS THAT ERR
  example = example %>% hoist(out, "result") %>% drop_na()
}