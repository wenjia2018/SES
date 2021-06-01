
#' ---
#' title: aging signature regression on color related snps
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
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
library(stringi)
walk(dir(path = here("R"),full.names = TRUE), source)
source("/home/xu/ses-1/user_wx/extract_v2.R")
load_data(reconciled = FALSE, remove_inflam = FALSE)
source("/home/xu/ses-1/user_wx/plotting_utils.R")
control = "ancestryPC_ses"
p_eqtl = 0.05
# skincolor_eqtl005_aging_composite_ancestry_11.05.2021.rds m7 m8 without genowide mediation
snp_NonHwhite <- readRDS("~/ses-1/user_wx/snp_NonHwhite_strata_eqtl005_aging_composite_ancestry_19.05.2021.rds")
snp <- readRDS("~/ses-1/user_wx/snp_eqtl005_aging_composite_ancestry_19.05.2021.rds")

#' ## among whites

#' ### 1 omnibus with aging composite ancestry controls
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=10, fig.height=8
outomni = p_eqtl %>% map(m8_present, control, snp_NonHwhite)

outomni %>%
  bind_rows() %>%
  dplyr::filter(p_omnibus < 0.05) %>%
  select(-gene_sig, -logFC) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_omnibus) %>%
  arrange(gene_set_name) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_omnibus=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### 2 PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE

outm71 = p_eqtl %>% map(outm7pca_anova, control, snp_NonHwhite)


outm71 %>%
  bind_rows() %>%
  unnest_longer(p_pca) %>% 
  mutate(p= p_pca %>% p.adjust("fdr")) %>% 
  filter(p<0.05) %>% 
  select(-coef, -loadings, - well_loaded) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_pca) %>%
  arrange(gene_set_name) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_pca=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ## blood sample

#' ### 3 omnibus with aging composite ancestry controls
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=10, fig.height=8
outomni = p_eqtl %>% map(m8_present, control, snp)

outomni %>%
  bind_rows() %>%
  dplyr::filter(p_omnibus < 0.05) %>%
  select(-gene_sig, -logFC) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_omnibus) %>%
  arrange(gene_set_name) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_omnibus=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()


#' ### 4 PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE

outm71 = p_eqtl %>% map(outm7pca_anova, control, snp)


outm71 %>%
  bind_rows() %>%
  unnest_longer(p_pca) %>% 
  mutate(p= p_pca %>% p.adjust("fdr")) %>% 
  filter(p<0.05) %>% 
  select(-coef, -loadings, - well_loaded) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_pca) %>%
  arrange(gene_set_name) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_pca=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()
