#' ---
#' title: results of removing inflam1k signature from other disease signatures, pca 6
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE
#+ echo=F, eval=T, warning=FALSE, message=FALSE
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
library(here)
walk(dir(path = here("R"),full.names = TRUE), source)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
example0 = readRDS("/home/share/scratch/xu/example0_w5bmi_removeinflam_6pc_trial2.rds")
threshold = 0.05/4/6
#' ## panel A results: FDR corrected whole genome
#+ echo=F, eval=T, warning=FALSE, message=FALSE
  
example0 %>%
  hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
  filter(controls=="all") %>% 
  dplyr::select(treatment, gene_set_name, p) %>% 
  dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>%  
  filter(p< 0.05)%>% 
    kableExtra::kable() %>%
    kableExtra::kable_styling()


#' ## panel B results: oblique rotation PCA-6pcs 
#' threshold to pick up significant PCs is 0.05/4/6
#+ echo=F, eval=T, warning=FALSE, message=FALSE
pc_ob <- example0 %>%
  hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
  unnest_longer(p) %>% 
  dplyr::select(treatment, gene_set_name, p, p_id) %>% 
  filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  filter(p< threshold)

pc_ob %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()


#' ## panel C results: oblique rotation PCA well loaded gene pathway reactome results
#+ echo=F, eval=T, warning=FALSE, message=FALSE
example0_noerror = remove_errors(example0)

example0_m7_ob = example0_noerror %>% get_sig_PCs_and_sig_enrichment_on_those_PCs("m7_ob", threshold = threshold)

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
  map_df(~as.data.frame(.x), .id="Signature") %>%
  as_tibble() %>%
  unique() %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()


#'  `rmarkdown::render("/home/xu/ses-1/user_wx/figure1.R")`
