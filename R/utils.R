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

get_table1 = function(example){ 
  
  tab1a = 
    example %>% 
    # hoist(out, "result") %>% 
    hoist(result, !!!funcs)  %>% 
    unnest(!!funcs) %>% 
    hoist(m1, pm1 = "p") %>% 
    hoist(m2, pm2 = "p") %>% 
    hoist(m3, pm3 = "p") %>% 
    hoist(m5, pm5 = "p") %>%
    hoist(m5b, pm5b = "p") %>%
    hoist(m8_fwer, m8_fwer_p = "p") %>% 
    hoist(m8_fdr, m8_fdr_p = "p") %>% 
    discard(is.list)
  
  # mediation
  tab1b = 
    example %>% 
    # hoist(out, "result") %>% 
    hoist(result, !!!funcs)  %>% 
    unnest(m99) %>%
    unnest_wider(m99) %>% 
    hoist(w5bmi, w5bmi_p = c("result", "p"))  %>% 
    hoist(bingedrink, bingedrink_p = c("result", "p"))  %>% 
    hoist(currentsmoke, currentsmoke_p = c("result", "p"))  %>% 
    hoist(phys_activ_ff5, phys_activ_ff5_p = c("result", "p"))  %>%   
    discard(is.list)
  
  (
    tab1a %>% left_join(tab1b)
  )
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