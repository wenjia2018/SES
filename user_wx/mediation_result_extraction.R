
#' ---
#' title: mediational analysis for 3 pca rotations 
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE


#' We do PCA for each disease signature using 3 rotations and then do mediational analysis
#'  for significant PCs using w5bmi as mediator
#'  
#' p: p value for indirected effect / mediated effect;
#' result_id: which pc is significant 
library(here)
library(tidyverse)
example4_1 = readRDS("/home/share/scratch/xu/example0_w5bmi.rds")

#' ## none rotation
m7_nn = example4_1 %>%
  hoist(out, "result") %>%
  hoist(result, "m7_nn","m7_vx","m7_ob") %>% 
  unnest(matches("^m7")) %>% 
  hoist(m7_nn, result = list("mediation","w5bmi", "result")) %>% 
  filter(result!="NULL") %>% 
  unnest_longer(result) %>% 
  hoist(result,p = "p") %>% 
  select(1:4,6)

m7_nn %>% kableExtra::kable() %>% kableExtra::kable_styling()

#' ## varimax rotation
m7_vx = example4_1 %>%
  hoist(out, "result") %>%
  hoist(result, "m7_vx","m7_nn","m7_ob") %>% 
  unnest(matches("^m7")) %>% 
  hoist(m7_vx, result = list("mediation","w5bmi", "result")) %>% 
  filter(result!="NULL") %>% 
  unnest_longer(result) %>% 
  hoist(result,p = "p")%>% 
  select(1:4,6)

m7_vx %>% kableExtra::kable() %>% kableExtra::kable_styling()
#' ## oblique rotation
m7_ob = example4_1 %>%
  hoist(out, "result") %>%
  hoist(result, "m7_ob","m7_vx","m7_nn") %>% 
  unnest(matches("^m7")) %>% 
  hoist(m7_ob, result = list("mediation","w5bmi", "result")) %>% 
  filter(result!="NULL") %>% 
  unnest_longer(result) %>% 
  hoist(result,p = "p")%>% 
  select(1:4,6)

m7_ob %>% kableExtra::kable() %>% kableExtra::kable_styling()

#'  `rmarkdown::render("/home/xu/ses-1/user_wx/mediation_result_extraction.R")`
