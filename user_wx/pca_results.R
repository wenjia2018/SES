#+ echo=F, eval=T
library(tidyverse)
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase)
library(ggformula)
# results for original 9 pc results
#+ echo=F, eval=T
example0 = readRDS("/home/share/scratch/example0.rds")
threshold = 0.0019
#' ## non rotation PCA
#+ echo=F, eval=F

# pc_nn <- example0 %>%
#   hoist(out, p = list("result", "m7_nn", 1, "p")) %>% 
#   unnest_longer(p) %>% 
#   dplyr::select(treatment, gene_set_name, p, p_id) %>% 
#   filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
#   filter(p<threshold)
# 
# pc_nn %>% 
#   kableExtra::kable() %>%
#   kableExtra::kable_styling()

#' \newpage
#' ## varimax rotation PCA
#+ echo=F, eval=F
# pc_vx <- example0 %>%
#   hoist(out, p = list("result", "m7_vx", 1, "p")) %>% 
#   unnest_longer(p) %>% 
#   dplyr::select(treatment, gene_set_name, p, p_id) %>% 
#   filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
#   filter(p<threshold)
# 
# pc_vx %>% 
#   kableExtra::kable() %>%
#   kableExtra::kable_styling()
#' 
# \newpage
#' ## oblique rotation PCA
#+ echo=F, eval=T
pc_ob <- example0 %>%
  hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
  unnest_longer(p) %>% 
  dplyr::select(treatment, gene_set_name, p, p_id) %>% 
  filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  filter(p< threshold)

pc_ob %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()


#' ## oblique rotation PCA well loaded gene pathway reactome results
#+ echo=F, eval=T

example0_m7_ob = readRDS("/home/xu/ses-1/user_wx/ob_rotation_pca9.rds")

m7_ob = example0_m7_ob %>%
  filter(map_lgl(well_loaded_genes_on_significant_PCs, ~ length(.x)!=0)) %>% 
  filter(treatment != "ses_composite_ff5")
 
all_pathway = m7_ob %>%
  unnest_longer(enrichment_of_well_loaded_genes_on_significant_PCs) %>%
  hoist(enrichment_of_well_loaded_genes_on_significant_PCs, pathway = list("out","enriched_physiology","reactome")) 

 all_pathway$pathway %>%
  map(function(x) as_tibble(x) %>% dplyr::select(Description, p.adjust)
      # %>% dplyr::top_n(-5)
      ) %>%
  set_names(all_pathway$gene_set_name) %>% 
  map_df(~as.data.frame(.x), .id="tag") %>%
  as_tibble() %>%
  unique() %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()


#'  `rmarkdown::render("/home/xu/ses-1/user_wx/pca_results.R")`
