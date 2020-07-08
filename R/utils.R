mediate = function(mediator){
  
  datt_m = 
    prepro(gene_set, treatment, c(controls, mediator)) %>%
    rename(treatment = treatment,
           mediator = mediator) %>% 
    remove_diseased_subjects_from_datt(gene_set_name, controls)
  
  ok = complete.cases(datt_m$mediator)
  datt_m = datt_m[ok, ]
  fit_m99(datt_m, gene_set) %>% extract_m99()  
}

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

