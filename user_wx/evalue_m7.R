
#' ---
#' title: Evalue for PCA
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE, include=FALSE, comment=NA

#' ### Evalue
#' *  The E-value is defined as the minimum strength of association, on the risk ratio scale,
#'  that an unmeasured confounder would need to have with both the treatment and the outcome to fully explain away a specific
#'  treatment–outcome association, conditional on the measured covariates.
#' * Input required to calculate Evalue: OLS estimate and standard deviation, and standard deviation of the model (outcome)
 

#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
library(tidyverse)
library(rlang)
library(skimr)
library(furrr)
library(limma)
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

example0_with1k <- readRDS("~/ses-1/user_wx/example_tmm_m7_withinflame.rds")
example0_without1k <- readRDS("~/ses-1/user_wx/example_tmm_m7_noinflame0209.rds")

threshold = 0.05
nfactors = 10

temp = 
  example0_with1k %>% 
  hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
  unnest_longer(p) %>% 
  mutate(p = p.adjust(p, method = "fdr"),
         # gene_set_name = str_c(gene_set_name,"_",p_id),
         "1KI Genes" = "With 1KI Genes") %>% 
dplyr::filter(p < threshold, p_id %in% c(str_c("d",1:nfactors))) %>% 
  hoist(out, evals = list("result", "m7_ob", 1, "evalue")) %>% 
  mutate(evalue = map2(.$p_id, .$evals, ~ .y[[.x]])) %>% 
  dplyr::select(treatment, gene_set_name, p, p_id, evalue)

evaluetable =
  temp %>% 
  mutate(`E-values` = map_dbl(.$evalue, ~ .[2,1])) %>% 
  mutate(pcmin = p_id %>% str_remove("d") %>% as.numeric()) %>%
  group_by(treatment, gene_set_name) %>%
  mutate(p_no = n()) %>% 
  slice(which.min(pcmin)) %>% 
  ungroup %>% 
  select(treatment, gene_set_name, `E-values`) %>% 
  mutate(
    treatment= case_when(treatment == "edu_p" ~ "Parental Education",
                         treatment =="income_pp1_log" ~  "Parental Income" ,
                         treatment =="SEI_max_p_w12" ~ "Parental SEI",
                         treatment =="ses_composite_pp1" ~ "Parental SES Composite",
                         treatment =="work_collar_rm_f12" ~ "Mother's Occupation",
                         treatment =="work_collar_rf_f12" ~ "Father's Occupation" ,
                         treatment =="work_collar_ff5" ~ "Occupation Work Collar",
                         treatment =="edu_max" ~ "Education" ,
                         treatment =="income_hh_ff5" ~ "Income"     ,
                         treatment =="SEI_ff5" ~ "Occupation"      ,
                         treatment =="ses_sss_composite" ~ "SES Composite"  ,
                         treatment =="sss_5" ~ "Subjective Social Status",
                         treatment =="ses_composite_ff5"  ~ "SES Composite 3")) 
  
evaluetable %>%
  pivot_wider(names_from = treatment, values_from = `E-values`) %>% 
  mutate(across(2:6, ~ .x %>% as.numeric() %>% format(digits = 4))) %>%
  mutate(across(everything(), ~ ifelse(.x=="   NA", "", .x))) %>% 
  # mutate(across(2:6, ~ .x %>% as.numeric() %>% format(digits = 4))) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()

temp2 = 
  example0_without1k %>% 
  hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
  unnest_longer(p) %>% 
  mutate(p = p.adjust(p, method = "fdr"),
         # gene_set_name = str_c(gene_set_name,"_",p_id),
         "1KI Genes" = "With 1KI Genes") %>% 
  dplyr::filter(p < threshold, p_id %in% c(str_c("d",1:nfactors))) %>% 
  hoist(out, evals = list("result", "m7_ob", 1, "evalue")) %>% 
  mutate(evalue = map2(.$p_id, .$evals, ~ .y[[.x]])) %>% 
  dplyr::select(treatment, gene_set_name, p, p_id, evalue)

evaluetable2=
  temp2 %>% 
  mutate(`E-values` = map_dbl(.$evalue, ~ .[2,1])) %>% 
  mutate(pcmin = p_id %>% str_remove("d") %>% as.numeric()) %>%
  group_by(treatment, gene_set_name) %>%
  mutate(p_no = n()) %>% 
  slice(which.min(pcmin)) %>% 
  ungroup %>% 
  select(treatment, gene_set_name, `E-values`) %>% 
  mutate(
         treatment= case_when(treatment == "edu_p" ~ "Parental Education",
                              treatment =="income_pp1_log" ~  "Parental Income" ,
                              treatment =="SEI_max_p_w12" ~ "Parental SEI",
                              treatment =="ses_composite_pp1" ~ "Parental SES Composite",
                              treatment =="work_collar_rm_f12" ~ "Mother's Occupation",
                              treatment =="work_collar_rf_f12" ~ "Father's Occupation" ,
                              treatment =="work_collar_ff5" ~ "Occupation Work Collar",
                              treatment =="edu_max" ~ "Education" ,
                              treatment =="income_hh_ff5" ~ "Income"     ,
                              treatment =="SEI_ff5" ~ "Occupation"      ,
                              treatment =="ses_sss_composite" ~ "SES Composite"  ,
                              treatment =="sss_5" ~ "Subjective Social Status",
                              treatment =="ses_composite_ff5"  ~ "SES Composite 3")) 
list(with1KI = evaluetable, without1KI = evaluetable2) %>% 
  openxlsx::write.xlsx("./user_wx/evalues_m7.xlsx")
evaluetable2 %>%
  pivot_wider(names_from = treatment, values_from = `E-values`) %>% 
  mutate(across(2:6, ~ .x %>% as.numeric() %>% format(digits = 4))) %>%
  mutate(across(everything(), ~ ifelse(.x=="   NA", "", .x))) %>% 
  # mutate(across(2:6, ~ .x %>% as.numeric() %>% format(digits = 4))) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()
#' ### Evalue for total effect with 1KI
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
for (i in 1 : dim(temp)[1]) {
  
  cat(" ############################################################","\n",
      "treatment is", temp[i, ]$treatment, "gene_set_name is", temp[i, ]$gene_set_name,"\n",
      "############################################################","\n")
  cat("total effect with 1KI","\n","\n")
  
  print(temp[i,]$evalue)

  cat("\n","\n")
  
}

#' ### Evalue for total effect without 1KI
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
for (i in 1 : dim(temp2)[1]) {
  
  cat(" ############################################################","\n",
      "treatment is", temp2[i, ]$treatment, "gene_set_name is", temp2[i, ]$gene_set_name,"\n",
      "############################################################","\n")
  cat("total effect without 1KI","\n","\n")
  
  print(temp2[i,]$evalue)
  cat("\n","\n")
  
}
