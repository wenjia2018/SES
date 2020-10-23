
#' ---
#' title: genowide regression and oblique rotation PCA analysis
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE


#' ## full control as in SES project: basic plus biological controls

#' ## genowide regression with 1KI

#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(here)
library(tidyverse)
example_race_with1KI <- readRDS("~/ses-1/user_wx/example_race_with1KI.rds")

example_race_with1KI %>%  
  hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
  filter(controls=="all") %>% 
  dplyr::select(treatment, gene_set_name, p) %>% 
  filter(p<0.05) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ## PCA regression with 1KI
#+ echo=F, eval=T, warning=FALSE, message=FALSE
threshold = 0.05/10
example_race_with1KI %>% 
  hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
  unnest_longer(p) %>% 
  dplyr::select(treatment, gene_set_name, p, p_id) %>% 
  dplyr::filter(p < threshold) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling() 

#' ## genowide regression without 1KI
#+ echo=F, eval=T, warning=FALSE, message=FALSE
example_race_without1KI <- readRDS("~/ses-1/user_wx/example_race_without1KI.rds")

example_race_without1KI %>%  
  hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
  filter(controls=="all") %>% 
  dplyr::select(treatment, gene_set_name, p) %>% 
  filter(p<0.05) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ## PCA regression without 1KI
#+ echo=F, eval=T, warning=FALSE, message=FALSE
threshold = 0.05/10
example_race_without1KI %>% 
  hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
  unnest_longer(p) %>% 
  dplyr::select(treatment, gene_set_name, p, p_id) %>% 
  dplyr::filter(p < threshold) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()
#+ echo=F, eval=T, warning=FALSE, message=FALSE
 # `rmarkdown::render("/home/xu/ses-1/user_wx/race_extraction_m7m8.R")`
