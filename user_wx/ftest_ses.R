
#' ---
#' title: Ftest for ses
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE, include=FALSE, comment=NA

#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
library(tidyverse)

#' ###  with 1KI
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
example_tmm_m12_withinflame_Ftest <- readRDS("~/ses-1/user_wx/example_tmm_m12_withinflame_Ftest.rds")
example_tmm_m12_noinflame_Ftest <- readRDS("~/ses-1/user_wx/example_tmm_m12_noinflame_Ftest.rds")

with =
  example_tmm_m12_withinflame_Ftest %>% 
  hoist(out, f = list("result", "m12_fdr", 1, "F_pval")) %>%
  hoist(out, siggene = list("result", "m12_fdr", 1, "sig_genes")) %>% 
  mutate(prop = f %>% map_dbl(~ sum(.x<0.05)/length(.x)),
         no_siggene = siggene %>% map_int(length)) %>% 
  select(-out)
with %>% select(-f) %>% kableExtra::kable() %>% kableExtra::kable_styling()
#' ###  without 1KI
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
without =
  example_tmm_m12_noinflame_Ftest %>% 
  hoist(out, f = list("result", "m12_fdr", 1, "F_pval")) %>%
  hoist(out, siggene = list("result", "m12_fdr", 1, "sig_genes")) %>% 
  mutate(prop = f %>% map_dbl(~ sum(.x<0.05)/length(.x)),
         no_siggene = siggene %>% map_int(length)) %>% 
  select(-out)
without %>% select(-f) %>% kableExtra::kable() %>% kableExtra::kable_styling()
