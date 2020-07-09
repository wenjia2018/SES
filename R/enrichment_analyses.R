
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
  
  # pass entrezgene_id to enrich function
  edo <- DOSE::enrichDGN(de)
  # plot
  barplot(edo, showCategory=20)
}


get_well_loaded_genes = 
  function(datt, gene_set, rotate, threshold = 0.1){
    # "threshold" appears to the default threshold for psych
    # for each PC dimension, get well-loaded genes (i.e. factor loading > 0.1)
    
    pca_rotated = fit_pca_util(datt, gene_set, rotate)
    non_sparse = abs(pca_rotated$loadings[,]) > threshold
    well_loaded_genes = non_sparse %>% as_tibble() %>% map(~rownames(non_sparse)[.x])
  }


get_well_loaded_genes_on_significant_PCs = function(example, m7_model, tabPCA){ 
  # get_well_loaded_genes_significant_genes(m7_model = "m7_ob")
  ids = c("treatment", "gene_set_name", "controls")
  binarize = . %>% 
    rowwise() %>% 
    mutate(across(matches("m6|m7"), ~ as.numeric(.x < 0.05))) 
  # select(tabPCA, all_of(ids), matches(m7_model)) 
  
  y = 
    example %>% 
    hoist(out, result = list("result" )) %>%
    unnest(result) %>% 
    unnest(matches("m7")) %>% 
    filter(names(m7_nn) == "other") %>% 
    pluck(m7_model) 
  x = 
    select(tabPCA, all_of(ids), matches(m7_model)) %>% 
    binarize() %>% 
    select(matches(m7_model)) %>%
    rowwise() %>%
    group_split() %>%
    map(unlist) %>%
    map(as.logical) 
  
  pmap(list(x = x, y = y), function(x,y)  y[x]) 
  
}



my_vis = function(DE_list){ 
  
  #  plaguerized from DE_enrichplot()
  
  results = readRDS("/home/share/preprocessed_two_batches/entrezgeneid.rds") 
  de = 
    enframe( unlist(DE_list), value = "gene") %>%
    left_join(results, by = c("gene" = "hgnc_symbol")) %>% 
    dplyr::pull(entrezgene_id)
  
  edo <- DOSE::enrichDGN(de)
  barplot(edo, showCategory=20)
  
}

enrichment_of_well_loaded_genes_on_significant_PCs =
  function(example, m7_model, tabPCA) {
  get_well_loaded_genes_on_significant_PCs(example, m7_model, tabPCA) %>% 
    discard(~length(.x) == 0) %>%
    map(my_vis) 
}
