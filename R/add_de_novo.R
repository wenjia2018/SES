add_de_novo = function(signatures, de_novo_treatments, control = "all", remove_inflam = TRUE){
  
  ttT_table = 
    args %>%
    filter(gene_set_name == "whole_genome_and_tfbm",
           treatment %in%  de_novo_treatment, 
           names(controls) == control) %>% 
    mutate(out = pmap(., safely(model_fit), funcs),
           controls = control) %>% 
    hoist(out, "result")
  
  
  de_novo_signatures = 
    ttT_table %>%
    pluck("result") %>%
    set_names(ttT_table$treatment %>% str_c("_","de_novo")) %>% 
    map(pluck("ttT")) %>% 
    map(~ filter(.x,adj.P.Val<=0.05) %>%
          pull(gene))
 
   if(remove_inflam == TRUE){
     
    sig = Reduce(union, signatures$outcome_set[table1[1:11]])
    de_novo_signatures = map(de_novo_signatures, ~ setdiff(.x, sig))
    
   }
  
  signatures$outcome_set = c(signatures$outcome_set, de_novo_signatures)

  return(signatures)
}

