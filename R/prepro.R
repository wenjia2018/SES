prepro = function(gene_set, treatment, controls){ 
  library(Biobase)
  library(tidyverse)
  
  genes = t(exprs(dat[gene_set, ]))
  # genes = genes %>% set_names( genes %>% names() %>% str_replace_all("-", "_"))
  pheno = pData(dat) %>% select(AID = AID, all_of(controls), all_of(treatment)) 
  datt <-pheno %>% 
    cbind(genes) %>% 
    as_tibble() %>% 
    select(-AID)
  
}