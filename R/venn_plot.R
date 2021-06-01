venn_plot = function(de_tfbm_result) {
  
  temp = de_tfbm_result %>%
    hoist(out, "result")
  
  result = temp %>%
    pluck("result") %>%
    set_names(temp$treatment)
  
  DE = result %>% 
    map(pluck("ttT")) %>% 
    map(~ filter(.x,adj.P.Val<=0.05) %>% pull(gene))
  # if redefine the name for each component
  # names(DE_removesig) = c("A", "B", "C")
  venn::venn(DE, zcolor = "style",
             ellipse = FALSE, opacity = 0.15, ilcs = 1,
             box = FALSE,
             sncs = 1)
}


gsea_webgestalt = function(treatment, ttT, file_output, enrichMethod = "GSEA", enrichDatabase="pathway_Reactome"){
  rankFile = str_c(file_output,"/", treatment,".rnk")
  
  ttT %>%
    dplyr::select(gene, logFC) %>% 
    write.table(file = rankFile, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  
  enrichResult = WebGestaltR::WebGestaltR(interestGeneFile = rankFile,
                                          interestGeneType = "genesymbol",
                                          enrichMethod = enrichMethod,
                                          organism = "hsapiens",
                                          enrichDatabase = enrichDatabase,
                                          fdrThr = 1,
                                          minNum=5,
                                          perNum = 1000,
                                          outputDirectory = file_output)
}
