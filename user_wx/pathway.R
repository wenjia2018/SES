#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(tidyverse)
library(here)
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase)
library(enrichplot)
library(gridExtra)
library(grid)
library(ggpubr)
library(dbr) # my package
m7_ob = readRDS("/home/xu/ses-1/user_wx/m7_ob.rds")
a = m7_ob %>%
  dplyr::select(well_loaded_genes_on_significant_PCs) %>%
  pull(well_loaded_genes_on_significant_PCs) %>%
  setNames(m7_ob$gene_set_name)

all_pathway = m7_ob %>%
  unnest_longer(enrichment_of_well_loaded_genes_on_significant_PCs) %>%
  hoist(enrichment_of_well_loaded_genes_on_significant_PCs, pathway = list("out","enriched_physiology","reactome")) 

temp_with = all_pathway$pathway %>%
  map(function(x) as_tibble(x) %>% dplyr::select(Description, p.adjust) %>% dplyr::top_n(-4)) %>%
  set_names(all_pathway$gene_set_name) %>% 
  .[c(1,2,5,6,14)] %>% 
  map_df(~as.data.frame(.x), .id="Signature") %>%
  as_tibble() %>% 
  unique() %>% 
  .[c(1,2,3,4,5,7,8,9,10,11,12,14,15,16),]

temp_with %>% kableExtra::kable() %>% kableExtra::kable_styling()

#'`rmarkdown::render("/home/xu/ses-1/user_wx/pathway.R")`

