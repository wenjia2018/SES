get_sig_PCs = function(example, m7_model, threshold){ 
  
  ( 
    tabPCA =
      example %>% 
      # hoist(out, result = list("result" )) %>%
      unnest(result) %>% 
      unnest(matches("m7")) %>%
      filter(names(m7_nn) == "p") %>%
      unnest_wider(m7_nn, names_sep = "_") %>% 
      unnest_wider(m7_vx, names_sep = "_") %>% 
      unnest_wider(m7_ob, names_sep = "_")
  )
  
  ############################################################
  # FOR EACH SIGNIFICANT PC, DO GENE ENRICHMENT ANALYSIS OF THE IMPORTANT GENES THEREIN.
  ############################################################
  ############################################################
  # 1. WHICH PCs HAVE BONFERONNI CORRECTED SIGNIFICANT RELATION TO TREATMENT?
  ############################################################
  
  # HOW MANY SIGNIFICANT PCs
  
  # n_args =   dim(tabPCA)[1]  # MORE CONSERVATIVE: number of treatment/outcome/control on RHS
  # n_args =   length(table1) # LESS CONSEVRATIVE: the number of signature sets
  # n_pc_dims = 9 # number of dimensions on LHS
  # bonferonni_threshold = 0.05 / (n_args * n_pc_dims) 
  # bonferonni_threshold = 0.05 / n_pc_dims
  bonferonni_threshold = threshold
  tabPCA %>% 
    rowwise(treatment, gene_set_name, controls) %>%
    summarize(sum = sum(c_across(matches(m7_model)) < bonferonni_threshold)) %>% 
    arrange(-sum) %>% 
    print(n = Inf) 
  
  # WHERE ARE THE SIGNIFICANT PCS?
  x = 
    tabPCA %>% 
    rowwise() %>%
    transmute(across(matches(m7_model), ~ . < bonferonni_threshold)) %>% 
    group_split() %>% 
    map(as_vector) %>% 
    map(unname)
  # WHICH GENES ARE WELL LOADED?
  y = 
    example %>% 
    unnest(result) %>% 
    unnest_wider(m7_model) %>% 
    hoist(other,"well_loaded") %>% 
    pluck("well_loaded")
  
  # TAKE WELL-LOADED GENES OF SIGNIFICANT PCs, THEN DO ENRICHMENT
  example = 
    example %>% 
    mutate(well_loaded_genes_on_significant_PCs = pmap(list(x = x, y = y), function(x,y)  y[x])) 
  
}