#' ---
#' title: tfbm for race and skincolor
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

#' #### tfbm: telis method and regression method
#' * telis over and telis under
#' * regression uni (logFC ~ tfbm1) and regression cov (logFC ~ tfbm1 + ... + tfbmN)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(tidyverse)
library(here)
source("/home/xu/ses-1/user_wx/extract_v2.R")
p_eqtl = c(0.01, 0.05)
control = "ancestryPC_ses"

race_bespoke_DE_22.03.2021 <- readRDS("~/ses-1/user_wx/race_bespoke_DE_22.03.2021.rds")
skincolor3_bespoke_DE_22.03.2021 <- readRDS("~/ses-1/user_wx/skincolor3_bespoke_DE_22.03.2021.rds")
data = race_bespoke_DE_22.03.2021 %>% 
  rbind(skincolor3_bespoke_DE_22.03.2021)
tfbm_race_skincolor = p_eqtl %>% map_df(~ outtfbm(.x, control, data)) %>% 
  # with fdr correction, no sig tfbm
  mutate_at(.vars = vars(c("tfbm")),
            .funs = list(~ map(., ~ mutate_if(., is.numeric, p.adjust, method = "fdr")))) %>%
  mutate(`tellis_p_under < 0.05` = tfbm %>% map(~ filter(.x, p_under < 0.05) %>% pull("tfbm")),
         `tellis_p_over < 0.05` = tfbm %>% map(~ filter(.x, p_over < 0.05) %>% pull("tfbm")),
         `regression p uni < 0.05` = tfbm %>% map(~ filter(.x, m_uni.p.value< 0.05) %>% pull("tfbm")),
         `regression p cov < 0.05` = tfbm %>% map(~ filter(.x, m_cov.p.value< 0.05) %>% pull("tfbm"))) %>% 
  select(-controls, -tfbm)
tfbm_race_skincolor %>% kableExtra::kable() %>% kableExtra::kable_styling()
