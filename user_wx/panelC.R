
#' ---
#' title: Figure 1 Panel C
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE


#' #### Panle C
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
walk(dir(path = here("R"),full.names = TRUE), source)

#' ## pathways with 1KI

#+ echo=F, eval=T, warning=FALSE, message=FALSE 


m7_ob = readRDS("/home/xu/ses-1/user_wx_RESTORED/m7_ob.rds")
# a = m7_ob %>%
#   dplyr::select(well_loaded_genes_on_significant_PCs) %>%
#   pull(well_loaded_genes_on_significant_PCs) %>%
#   setNames(m7_ob$gene_set_name)
#   # not correct for 2 pc sig cases
# a[c(1,2,5,6,12)] %>% openxlsx::write.xlsx("./user_wx/well_loaded_gene_withinflam.xlsx")

all_pathway = m7_ob %>%
  unnest_longer(enrichment_of_well_loaded_genes_on_significant_PCs) %>%
  hoist(enrichment_of_well_loaded_genes_on_significant_PCs, pathway = list("out","enriched_physiology","reactome")) 

names = all_pathway$gene_set_name

temp_with = all_pathway$pathway %>%
  map(function(x) as_tibble(x) %>%
        dplyr::select(Description, p.adjust) #%>%
      # dplyr::top_n(-4)
  ) #%>%
# set_names(all_pathway$gene_set_name) %>% 
# .[c(1,2,5,6,14)] %>% 
# map_df(~as.data.frame(.x), .id="Signature") %>%
# as_tibble() %>% 
# unique() %>% 
# .[c(1,2,3,4,5,7,8,9,10,11,12,14,15,16),]
# depression
# depr=temp_with[c(2,4,7,10,16)] %>% bind_rows() %>% dplyr::distinct(Description, .keep_all = TRUE)
depr = temp_with[c(2,4,7,10,16)] %>%
  bind_rows() %>% 
  arrange(Description, p.adjust) %>% 
  dplyr::distinct(Description, .keep_all = TRUE) %>% 
  arrange(p.adjust) %>% 
  mutate_at(.vars =c("p.adjust"),
            .funs = list(~ .x %>% format(digits = 3, scientific =T))) 
#' ### depression 1 pc

temp_with[c(2)] %>% kableExtra::kable() %>% kableExtra::kable_styling()
# copd

#+ echo=F, eval=T, warning=FALSE, message=FALSE
copd = temp_with[c(1,3,8,9,15)] %>%
  bind_rows() %>%
  arrange(Description, p.adjust) %>% 
  dplyr::distinct(Description, .keep_all = TRUE) %>% 
  arrange(p.adjust) %>% 
  mutate_at(.vars =c("p.adjust"),
            .funs = list(~ .x %>% format(digits = 3, scientific =T))) 

#' ### copd two PCs
#' 
temp_with[c(9,15)] %>% kableExtra::kable() %>% kableExtra::kable_styling()
# inflam1k

#+ echo=F, eval=T, warning=FALSE, message=FALSE

i1k=temp_with[c(5,11,12)] %>%
  bind_rows() %>%
  arrange(Description, p.adjust) %>% 
  dplyr::distinct(Description, .keep_all = TRUE) %>% 
  arrange(p.adjust) %>% 
  mutate_at(.vars =c("p.adjust"),
            .funs = list(~ .x %>% format(digits = 3, scientific =T))) 

#' ### 1KI two pcs
temp_with[c(5,11)] %>% kableExtra::kable() %>% kableExtra::kable_styling()


#' ### Asthma and Arthritis has no pathways with adjusted p value under 0.05
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
m7_ob = readRDS("/home/xu/ses-1/user_wx/m7_ob_fullpathways.rds")
all_pathway = m7_ob %>%
  unnest_longer(enrichment_of_well_loaded_genes_on_significant_PCs) %>%
  hoist(enrichment_of_well_loaded_genes_on_significant_PCs, pathway = list("out","enriched_physiology","reactome")) 

names = all_pathway$gene_set_name

temp_with = all_pathway$pathway %>%
  map(function(x) as_tibble(x) %>%
        dplyr::select(Description, p.adjust)
      ) 

#' ### Asthma
temp_with[14][[1]] %>%  head(10) %>% kableExtra::kable() %>% kableExtra::kable_styling()

#' ### Arthritis
temp_with[6][[1]] %>%  head(10) %>% kableExtra::kable() %>% kableExtra::kable_styling()

# extract genes
gene_list = all_pathway %>%
  dplyr::select(well_loaded_genes_on_significant_PCs) %>%
  pull(well_loaded_genes_on_significant_PCs) 
  # not correct for 2 pc sig cases
list(Asthma = gene_list[14], Arthritis = gene_list[6]) %>%
  openxlsx::write.xlsx("./user_wx/gene_nonsig_pathway_withinflam.xlsx")




#' ## pathways without 1KI
#' 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m7_ob = readRDS("/home/xu/ses-1/user_wx/m7_ob_without1KI.rds")
b = m7_ob %>%
  dplyr::select(well_loaded_genes_on_significant_PCs) %>%
  pull(well_loaded_genes_on_significant_PCs) %>%
  setNames(m7_ob$gene_set_name)

# b[c(1,3,4,5,8)] %>% openxlsx::write.xlsx("./user_wx/well_loaded_gene_noinflam.xlsx")

all_pathway = m7_ob %>%
  unnest_longer(enrichment_of_well_loaded_genes_on_significant_PCs) %>%
  hoist(enrichment_of_well_loaded_genes_on_significant_PCs, pathway = list("out","enriched_physiology","reactome")) 
names = all_pathway$gene_set_name
temp_without = all_pathway$pathway %>%
  map(function(x) as_tibble(x) %>% dplyr::select(Description, p.adjust))
#     %>% dplyr::top_n(-2)) %>%
# set_names(all_pathway$gene_set_name) %>% 
# .[c(1,3,4,5,8)] %>% 
# map_df(~as.data.frame(.x), .id="Signature") %>%
# as_tibble() %>% 
# unique()

# copd

#+ echo=F, eval=T, warning=FALSE, message=FALSE
copd = temp_without[c(1,2,6,11)] %>%
  bind_rows() %>%
  arrange(Description, p.adjust) %>% 
  dplyr::distinct(Description, .keep_all = TRUE) %>% 
  arrange(p.adjust) %>% 
  mutate_at(.vars =c("p.adjust"),
            .funs = list(~ .x %>% format(digits = 3, scientific =T))) 

#' ### copd
temp_without[c(1)] %>% kableExtra::kable() %>% kableExtra::kable_styling()
# arthritis

#+ echo=F, eval=T, warning=FALSE, message=FALSE
arthritis = temp_without[c(4,9,13)] %>%
  bind_rows() %>%
  arrange(Description, p.adjust) %>% 
  dplyr::distinct(Description, .keep_all = TRUE) %>% 
  arrange(p.adjust) %>% 
  mutate_at(.vars =c("p.adjust"),
            .funs = list(~ .x %>% format(digits = 3, scientific =T))) 

#' ### arthritis
temp_without[c(4)]  %>% kableExtra::kable() %>% kableExtra::kable_styling()


hypertension = temp_without[c(8)] %>%
  bind_rows() %>%
  arrange(Description, p.adjust) %>% 
  dplyr::distinct(Description, .keep_all = TRUE) %>% 
  arrange(p.adjust) %>% 
  mutate_at(.vars =c("p.adjust"),
            .funs = list(~ .x %>% format(digits = 3, scientific =T))) 
#' ### hypertension
#' 
temp_without[c(8)] %>% kableExtra::kable() %>% kableExtra::kable_styling()


#' ### depression and Asthma has no pathways with adjusted p value under 0.05 

#+ echo=F, eval=T, warning=FALSE, message=FALSE
m7_ob = readRDS("/home/xu/ses-1/user_wx/m7_ob_without1KI_fullpathways.rds")
all_pathway = m7_ob %>%
  unnest_longer(enrichment_of_well_loaded_genes_on_significant_PCs) %>%
  hoist(enrichment_of_well_loaded_genes_on_significant_PCs, pathway = list("out","enriched_physiology","reactome")) 

names = all_pathway$gene_set_name

temp_without = all_pathway$pathway %>%
  map(function(x) as_tibble(x) %>%
        dplyr::select(Description, p.adjust)
  ) 

#' ### Asthma
temp_without[5][[1]] %>%
  head(10) %>%
  kableExtra::kable() %>% kableExtra::kable_styling()

#' ### Depression
temp_without[3][[1]] %>%  head(10) %>% kableExtra::kable() %>% kableExtra::kable_styling()


# extract genes
gene_list = all_pathway %>%
  dplyr::select(well_loaded_genes_on_significant_PCs) %>%
  pull(well_loaded_genes_on_significant_PCs) 
# not correct for 2 pc sig cases
list(Asthma = gene_list[5], Depression = gene_list[3],
     Hypertension = gene_list[8], Arthritis = gene_list[4]) %>%
  openxlsx::write.xlsx("./user_wx/gene_nonsig_pathway_noinflam.xlsx")

list(COPD = gene_list[1], Depression = gene_list[3], Arthritis = gene_list[4],
     Asthma = gene_list[5], 
     Hypertension = gene_list[8]) %>%
  openxlsx::write.xlsx("./user_wx/well_loaded_gene_noinflam.xlsx")

#'`rmarkdown::render("/home/xu/ses-1/user_wx/panelC.R")`