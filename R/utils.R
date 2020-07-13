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