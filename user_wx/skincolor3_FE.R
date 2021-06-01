
#' ---
#' title: skin color 3 levels with basic control
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE


#' ## Controls, treatment and threshold
#' 
#' ### controls:
#'   basic = 
#' c(
#'  "sex_interv", "famid_fullsib",
#'   "Plate", "AvgCorrelogram100" ,"age_w5",
#'   # over-representation 
#'  "NK.cells.activated",
#'  "T.cells.CD8",
#'  # under-representation
#'  "Macrophages.M0", 
#'  "Macrophages.M2",
#'  "B.cells.naive",
#'  "T.cells.CD4.memory.resting"
#' )
#' 

#' ### treatments:
#' 
#' * skin color 3 levels: white as reference, LightMed and DarkBlack as treatments

#' 
#' ## omnibus regression p values are corrected genowide

#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(here)
library(tidyverse)
example_race_with1KI <- readRDS("~/ses-1/user_wx/example_skincolor_FE_08.04.2021.rds")


example_race_with1KI %>% 
  hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
  hoist(out, ttT = list("result", "m8_fdr", 1, "other", "m")) %>% 
  filter(ttT != "NULL") %>% 
  mutate(gene_sig = map(ttT, ~ dplyr::filter(., adj.P.Val<0.05) %>% pull(gene))) %>% 
  dplyr::select(treatment, gene_set_name, controls, p, gene_sig) %>% 
  dplyr::filter(p<0.05) %>%
  # left_join(race_ge_tfbm, by = c("treatment", "gene_set_name", "controls")) %>% 
  rename(p_omnibus = p) %>% 
  filter(controls == "basic") %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ## PCA regression p values uncorrected, but only the cases where the p value is less than 0.05/10 are presented
#' (no p <0.005)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
threshold = 0.05/10
threshold_med = 0.05


example_race_with1KI %>%
  hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
  unnest_longer(p) %>%
  dplyr::select(treatment, gene_set_name, controls, p, p_id) %>% 
  dplyr::filter(p < threshold) %>% 
  kableExtra::kable() %>%
  # kableExtra::column_spec(column = 7, width = "16in", width_min="8in") %>% 
  kableExtra::kable_styling() 



