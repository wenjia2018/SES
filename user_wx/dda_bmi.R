
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

#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
#' ### perform DDAmodel selection
#' * evaluate distributional properties of observed variables adjusting for covariates
#' * evaluate distributional properties of residuals
#' * evaluate independence properties of predictors and residuals

example_bmi_m7_DDA_30.07.2021 <- readRDS("~/ses-1/user_wx/example_bmi_m7_DDA_30.07.2021.rds")
example_bmi_m7_DDA_10.08.2021_missingfrom3007 <- readRDS("~/ses-1/user_wx/example_bmi_m7_DDA_10.08.2021_missingfrom3007.rds")
example_bmi_m7_DDA_11.08.2021_missingfrom1008 <- readRDS("~/ses-1/user_wx/example_bmi_m7_DDA_11.08.2021_missingfrom1008.rds")


part3 = 
  example_bmi_m7_DDA_11.08.2021_missingfrom1008 %>%
  mutate(id = str_c(treatment, gene_set_name))

part2 = 
  example_bmi_m7_DDA_10.08.2021_missingfrom3007 %>% 
  mutate(id = str_c(treatment, gene_set_name)) %>% 
  filter(!(id %in% part3$id)) %>% 
  rbind(part3)

part1 = 
  example_bmi_m7_DDA_30.07.2021 %>% 
  mutate(id = str_c(treatment, gene_set_name)) %>% 
  filter(!(id %in% part2$id)) %>% 
  rbind(part2)


old =
  example_bmi_m7_DDA_30.07.2021 %>%
  hoist(out, estimate = list("result", "m7_ob", 1, "detail", "t")) %>%
  unnest_longer(estimate) %>%
  hoist(estimate, p = "p.value") %>%
  hoist(estimate, coef = "estimate") %>%
  filter(p<0.005) %>%
  dplyr::select(treatment, gene_set_name, p, coef, estimate_id) %>% 
  filter(p < 0.05/50)
#  old %>% kableExtra::kable() %>% kableExtra::kable_styling()
temp = 
  part1 %>%
  hoist(out, dda = list("result", "m7_ob", 1, "dda", "w5bmi_lm", "result")) %>%
  unnest_longer("dda") %>%
  filter(!is.na(dda_id)) %>% 
  dplyr::select(treatment, gene_set_name, dda, dda_id) 

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
