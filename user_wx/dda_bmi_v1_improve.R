
#' ---
#' title: Results of direction dependency
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
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase)
library(ggformula)
library(ggpubr)
library(here)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(tidyverse)
library(stringi)
walk(dir(path = here("R"),full.names = TRUE), source)


example_bmi_m7_DDA <- readRDS("~/ses-1/user_wx/example_bmi_m7_DDA_30.07.2021.rds")
example_bmi_m7_DDA_10.08.2021_missingfrom3007 <- readRDS("~/ses-1/user_wx/example_bmi_m7_DDA_10.08.2021_missingfrom3007.rds")
example_bmi_m7_DDA_11.08.2021_missingfrom1008 <- readRDS("~/ses-1/user_wx/example_bmi_m7_DDA_11.08.2021_missingfrom1008.rds")

example_bmi_m7_DDA =
  example_bmi_m7_DDA %>%
  mutate(out2 = hoist(., out, keep = list("result", "m7_ob", 1, "dda", "w5bmi_lm", "result"))$keep) %>%
  filter(!map_lgl(.$out2, ~ is.null(.x))) %>% 
  dplyr::select(-out2) %>%
  rbind(example_bmi_m7_DDA_10.08.2021_missingfrom3007) %>%
  mutate(out2 = hoist(., out, keep = list("result", "m7_ob", 1, "dda", "w5bmi_lm", "result"))$keep) %>%
  filter(!map_lgl(.$out2, ~ is.null(.x))) %>% 
  dplyr::select(-out2) %>%
  rbind(example_bmi_m7_DDA_11.08.2021_missingfrom1008)


old =
  example_bmi_m7_DDA %>%
  hoist(out, estimate = list("result", "m7_ob", 1, "detail", "t")) %>%
  unnest_longer(estimate) %>%
  hoist(estimate, p = "p.value") %>%
  hoist(estimate, coef = "estimate") %>%
  dplyr::select(treatment, gene_set_name, p, coef, estimate_id) %>% 
  filter(p < 0.05/50)
#  old %>% kableExtra::kable() %>% kableExtra::kable_styling()
temp = 
  example_bmi_m7_DDA %>% 
  hoist(out, dda = list("result", "m7_ob", 1, "dda", "w5bmi_lm", "result")) %>%
  dplyr::select(treatment, gene_set_name, dda) %>% 
  filter(!map_lgl(.$dda, ~ is.null(.x))) %>%
  unnest_longer(dda) 


dda = temp %>% inner_join(old %>% select(1,2,5), by = c("treatment"="treatment","gene_set_name" = "gene_set_name","dda_id" = "estimate_id"))


# temp %>% filter(treatment %>% str_detect("edu") & gene_set_name %>% str_detect("COPD")) %>% `[[`("dda") %>% `[[`(1) %>% `[[`("d3")


for (i in 1 : dim(dda)[1]) {
  
  cat(" ############################################################","\n",
      "treatment is", dda[i, ]$treatment, "gene_set_name is", dda[i, ]$gene_set_name,"\n",
      "############################################################","\n")
  # print(dda[i,]$dda)
  print(dda[i, ] %>% hoist(dda, "vardist") %>% pluck("vardist"))
  print(dda[i, ] %>% hoist(dda, "resdist") %>% pluck("resdist"))
  print(dda[i, ] %>% hoist(dda, "indep4") %>% pluck("indep4"))
}
