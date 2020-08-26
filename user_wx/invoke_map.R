
method = list(DOSE = DOSE::enrichDGN,Reactome =ReactomePA::enrichPathway, Kegg = clusterProfiler::enrichKEGG)
pathway = function(de, method){
  edo = map(de, method)
  barplot(edo, showCategory=20)
  pathway = edo@result %>% filter(p.adjust<=0.05)
}

method = list(DOSE = DOSE::enrichDGN,Reactome =ReactomePA::enrichPathway, Kegg = clusterProfiler::enrichKEGG)
edo = method %>% invoke_map(gene=de)
pathway = edo %>% map(~ .x@result %>% filter(p.adjust<=0.05))
map2(edo, names(edo), barplot(.x), showCategory=20, title = .y)
edo %>% map(barplot, showCategory=20)