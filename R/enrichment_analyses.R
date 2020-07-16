
DE_enrichplot = function(ttT){
  
  if (reproducible <- FALSE){
    # load our whole gene list 
    dt_batches1_2_recon <- readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_25.06.2020.rds")
    dt_batches1_2_unrecon <- readRDS("/home/share/preprocessed_two_batches/wx/dt_batches1_2_steve_waves_21042020.rds")
    
    Glist= union(featureNames(dt_batches1_2_recon), featureNames(dt_batches1_2_unrecon)) 
    # convert our whole gene list from hugo name to entrez
    mart = useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl')
    attributes = c('entrezgene_accession','entrezgene_id','hgnc_symbol')
    filters = 'entrezgene_accession'
    results = getBM(mart=mart,
                    attributes = attributes,
                    filters = filters,
                    values = Glist)
    # save genelist with both names for future lookup
    results %>% saveRDS("/home/share/preprocessed_two_batches/entrezgeneid.rds")
  }
  
  
  # load our whole gene name list
  results =  readRDS("/home/share/preprocessed_two_batches/entrezgeneid.rds")
  DE_list = ttT %>% filter(adj.P.Val <= 0.05)
  
  # coverting to entrezgene id
  de = DE_list %>%
    left_join(results, by = c("gene" = "hgnc_symbol")) %>% 
    dplyr::pull(entrezgene_id)
  
  # pass entrezgene_id to enrich function, then plot
  
  # edo_DOSE <- DOSE::enrichDGN(de)
  # fig_DOSE = barplot(edo_DOSE, showCategory=20, title = "DOSE")
  
  edo_reactome <- ReactomePA::enrichPathway(de)
  fig_reactome = barplot(edo_reactome, showCategory=20, title = "Reactome")
  
  edo_kegg <- clusterProfiler::enrichKEGG(de)
  fig_kegg = barplot(edo_kegg, showCategory=20, title = "Kegg")
  
  list("clusterProfiler", "DO.db", "ReactomePA", "reactome.db", "DOSE", "graphite", "enrichplot",  "GO.db", "GOSemSim", "AnnotationDbi") %>% map(detach_package)
  
  fig = list(reactome = fig_reactome, kegg = fig_kegg)
  return (fig)
}


get_well_loaded_genes = 
  function(datt, gene_set, rotate, loading_threshold = 0.1){
    # "loading_threshold" appears to the default threshold for psych
    # for each PC dimension, get well-loaded genes (i.e. factor loading > 0.1)
    
    pca_rotated = fit_pca_util(datt, gene_set, rotate)
    non_sparse = abs(pca_rotated$loadings[,]) > loading_threshold
    well_loaded_genes = non_sparse %>% as_tibble() %>% map(~rownames(non_sparse)[.x])
  }



my_vis = function(DE_list, p_val_threshold = 0.05){ 
  
  # see also DE_enrichplot()
  print(str_c(runif(1), "  please wait: calculating enrichment for the well-loaded genes of significant PCs ..."))
  results = readRDS("/home/share/preprocessed_two_batches/entrezgeneid.rds") 
  de = 
    enframe(DE_list, value = "gene") %>%
    left_join(results, by = c("gene" = "hgnc_symbol")) %>% 
    dplyr::pull(entrezgene_id)
  
  edo_reactome <- ReactomePA::enrichPathway(de)
  fig_reactome = barplot(edo_reactome, showCategory=20, title = "Reactome")
  enriched_physiology_reactome = edo_reactome@result %>% filter(p.adjust <= p_val_threshold) %>% pull(Description)  
  
  edo_kegg <- clusterProfiler::enrichKEGG(de)
  fig_kegg = barplot(edo_kegg, showCategory=20, title = "Kegg")
  enriched_physiology_kegg = edo_kegg@result %>% filter(p.adjust <= p_val_threshold) %>% pull(Description)  
  
  fig = list(reactome = fig_reactome, kegg = fig_kegg)
  enriched_physiology = list(reactome = enriched_physiology_reactome, kegg = enriched_physiology_kegg)
  
  list("clusterProfiler", "DO.db", "ReactomePA", "reactome.db", "DOSE", "graphite", "enrichplot",  "GO.db", "GOSemSim") %>% map(detach_package)
  
 
  list(out = list(fig = fig, enriched_physiology = enriched_physiology))
}



get_sig_PCs_and_sig_enrichment_on_those_PCs = function(example, m7_model){ 
  
  ( 
    tabPCA =
      example %>% 
      # hoist(out, result = list("result" )) %>%
      unnest(result) %>% 
      unnest(matches("m6|m7")) %>%
      filter(names(m7_nn) == "p") %>%
      unnest(matches("m6")) %>%
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
  
  n_args =   dim(tabPCA)[1]  # MORE CONSERVATIVE: number of treatment/outcome/control on RHS
  n_args =   length(table1) # LESS CONSEVRATIVE: the number of signature sets
  n_pc_dims = 9 # number of dimensions on LHS
  bonferonni_threshold = 0.05 / (n_args * n_pc_dims) 
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
    pluck("other") 
  
  # TAKE WELL-LOADED GENES OF SIGNIFICANT PCs, THEN DO ENRICHMENT
  example = 
    example %>% 
    mutate(well_loaded_genes_on_significant_PCs = pmap(list(x = x, y = y), function(x,y)  y[x]),
           enrichment_of_well_loaded_genes_on_significant_PCs = map_depth(well_loaded_genes_on_significant_PCs, 2, my_vis)) 
  
}