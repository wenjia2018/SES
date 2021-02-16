
#' ---
#' title: race pattern with basic control, ses composite control and ancestryPC control(customised by us)
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

#' ### 3 sets of controls:
#' * basic control as in SES project
#' * ses control: basic + ses
#' * ancestryPC control : basic + ses composite + ancestryPC
#'   + ancestryPC is calculated using snps associated (p<0.05) with genes in our disease sets(~ 1800 genes)
#'
#' ### two focal treatment:
#'
#' * white_black  
#' * white_hispanic
#' * Non hispanic white is the reference group
#' 
#' ###  For pca results, only the cases where the p value is less than 0.05/10 are presented


#' ## omnibus regression and self contained gene set test (mroast) results

#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(here)
library(tidyverse)
example_race_with1KI <- readRDS("~/ses-1/user_wx/example_race_dummy_26.11.rds")
race_ge_tfbm <- readRDS("~/ses-1/user_wx/race_detfbm_dummy_26.11.rds")

# example_race_with1KI <- readRDS("~/ses-1/user_wx/example_race_with1KI_1030_withoutses4.rds")
# race_ge_tfbm <- readRDS("~/ses-1/user_wx/example_race_detfbm_withoutses4.rds")

race_ge_tfbm = race_ge_tfbm %>%
  hoist(out, gsea = list("result", "gene_set_test")) %>% 
  dplyr::select(1,3,4) %>%
  unnest(gsea) 

example_race_with1KI %>% 
  hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
  dplyr::select(treatment, gene_set_name, controls, p) %>% 
  dplyr::filter(p<0.05) %>%
  left_join(race_ge_tfbm, by = c("treatment", "gene_set_name", "controls")) %>% 
  rename(p_omnibus = p) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ## PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
threshold = 0.05/10
threshold_med = 0.05
var = example_race_with1KI %>%
  hoist(out, var_explained = list("result", "m7_ob", 1, "other", "varexplained")) %>% 
  dplyr::select(1,2,3,4)

var$var_explained = var$var_explained %>% map(~ set_names(.x, str_c("d", 1:10)))

var = var %>% unnest_longer(var_explained)


gene_list = example_race_with1KI %>%
  hoist(out, well_loaded = list("result", "m7_ob", 1, "other", "well_loaded")) %>% 
  dplyr::select(1,2,3,4)

gene_list$well_loaded = gene_list$well_loaded %>% map(~ set_names(.x, str_c("d", 1:10)))

gene_list = gene_list %>% unnest_longer(well_loaded)

example_race_with1KI %>%
  hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
  unnest_longer(p) %>%
  dplyr::select(treatment, gene_set_name, controls, p, p_id) %>% 
  dplyr::filter(p < threshold) %>%
  left_join(var, by = c("treatment", "gene_set_name", "controls", "p_id"= "var_explained_id")) %>% 
  left_join(gene_list, by = c("treatment", "gene_set_name","controls", "p_id"= "well_loaded_id")) %>% 
  dplyr::select(1:6) %>% 
  kableExtra::kable() %>%
  # kableExtra::column_spec(column = 7, width = "16in", width_min="8in") %>% 
  kableExtra::kable_styling() 


