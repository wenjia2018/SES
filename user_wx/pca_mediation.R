mediate_m7 = function(datt, gene_set, rotate){
  
  pca_rotated = fit_pca_util(datt, gene_set, rotate) 
  datt = dplyr::select(datt, -gene_set) 
  gene_set = colnames(pca_rotated$scores)
  
  outcome = gene_set %>% set_names()
      
  datt_pca = bind_cols(datt, !!outcome := pca_rotated$scores[, outcome])
      
}