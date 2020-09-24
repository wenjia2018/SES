netplot = function(DE){
  results =  readRDS("/home/share/preprocessed_two_batches/entrezgeneid.rds")
  # coverting to entrezgene id
  DE = data.frame("hgnc_symbol" = DE)
  DE = DE %>% 
    left_join(results, by =  "hgnc_symbol") %>% 
    dplyr::pull(entrezgene_id)
  x = ReactomePA::enrichPathway(DE)
  edox = DOSE::setReadable(x, 'org.Hs.eg.db', 'ENTREZID')
  enrichplot::cnetplot(edox, categorySize="pvalue", foldChange = DE)
}
