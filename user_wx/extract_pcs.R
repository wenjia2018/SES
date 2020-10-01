

set.seed(123)
library(here)
library(tidyverse)
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase)
library(enrichplot)
library(dbr) # my package

walk(dir(path = here("R"),full.names = TRUE), source)
fit_m4 = partial(fit_m4, n_perm = 1000) # specify n_perm
library(tidyverse)
example0 =  readRDS("/home/share/scratch/example0_without_1KI_de_novo2.rds")

source("./user_wx/get_sig_PCs.R")
m7_ob = example0 %>% 
  remove_errors() %>% 
  get_sig_PCs("m7_ob", threshold = 0.05/10/5) %>% 
  filter(map_lgl(well_loaded_genes_on_significant_PCs, ~ length(.x)!=0))

income = m7_ob %>% filter(gene_set_name=="income_unique_de_novo")
ses4 = m7_ob %>% filter(gene_set_name=="ses4_unique_de_novo")
ses4income = m7_ob %>% filter(gene_set_name=="ses4income_de_novo")

income_wellloaded_gene = income$well_loaded_genes_on_significant_PCs[1] %>% unlist(recursive=FALSE)
ses4_wellloaded_gene = ses4$well_loaded_genes_on_significant_PCs[5] %>% unlist(recursive=FALSE)
ses4income_wellloaded_gene = ses4income$well_loaded_genes_on_significant_PCs[2] %>% unlist(recursive=FALSE)

income_wellloaded_gene  %>% openxlsx::write.xlsx("./user_wx/income_unique_denovo_wellloaded_gene.xlsx")
ses4_wellloaded_gene %>% openxlsx::write.xlsx("./user_wx/ses4_unique_denovo_wellloaded_gene.xlsx")
ses4income_wellloaded_gene  %>% openxlsx::write.xlsx("./user_wx/income_ses4_denovo_wellloaded_gene.xlsx")


