
#' ---
#' title: race with aging cluster complement outcome
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

#' ### treatment--raceethnicity: 
#' * "NonHwhite" as reference
#' * "NonHblack",
#' * "Hispanic"

control = "ancestryPC_ses"
threshold = 0.05/10
threshold_med = 0.05
#' ### controls: 
#' 
#' * ancestryPC_ses : basic + ses + bespoke ancestryPC
#' 
#' ### outcome:
#'c(
#'  "ctra",
#'  "inflamation1k",
#' "aging"
#')
#'
#' ### analysis:
#' * omnibus P value: association between each disease set and color related snps, 
#' the smallest, whole-genome FDR corrected p-value via whole genome regressions
#'
#' * PCA p value: p value associated between each PC of the disease set and color related snps
#' 
#' * p_eqtl is a sequence of p (0.05, 0.01, 1e-3,...,1e-10) used to choose asscoiated snps
#' 
#' * omnibus p values are genowide corrected p values
#' 
#' * pca regression p value and mediation p values are unadjusted 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(here)
library(tidyverse)
source("/home/xu/ses-1/user_wx/extract_v2.R")
# bespoke <- readRDS("/home/xu/ses-1/user_wx/race_bespoke_12.02.2021.rds")
bespoke <- readRDS("/home/xu/ses-1/user_wx/race_bespoke_15.03.2021.rds")

mediators = 
  c(
    "stress_perceived_lm",
    "bills_binary",
    "currentsmoke_binary",
    "w5bmi_lm",
    "insurance_lack_binary",
    "lowbirthweight_binary",
    "high_lowbirth_binary",
    
    "totdiscrim2_binary",  # binary
    "discrim2_binary",   # binary
    "totdiscrim1_category",   # categorical of 4
    # special treatment in mediation
    "totdiscrim1_gamma", #exponential fit glm for mediation
    "totdiscrim2_gamma", #exponential fit glm for mediation
    "countdiscrimwhy_gamma",#exponential fit glm for mediation
    
    "totdiscrim1_pois", #poisson fit glm for mediation
    "totdiscrim2_pois", #poisson fit glm for mediation
    "countdiscrimwhy_pois",#poisson fit glm for mediation
    
    
    "totdiscrim2_category", #ordered logistic fit glm for mediation
    "countdiscrimwhy_category"#ordered logistic fit glm for mediation
  )


#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
p_eqtl = c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)

outomni = p_eqtl %>% map(m8_present, control, bespoke)

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

#' ### mean of all the omnibus significant genes mediation analysis
#+ echo=F, eval=T, warning=FALSE, message=FALSE
outomnimed = p_eqtl %>% map(outm8med, control, bespoke)
outomnimed %>% 
  bind_rows() %>%
  dplyr::filter(p_med < 0.05) %>%
  select(-table1, -controls, -med, -control_set) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_med) %>%
  arrange(gene_set_name) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_omnibus=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()




#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE

# outm7pca = function(p, control){
#   m7_present(f0(p), control) %>% mutate(p_eqtl = p)
# }

outm71 = p_eqtl %>% map(outm7pca, control, bespoke)


outm71 %>%
  bind_rows() %>%
  dplyr::filter(p_pca < threshold) %>%
  select(-coef, -loadings, - well_loaded) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_pca) %>%
  arrange(gene_set_name) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_pca=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()
#' ### mediation 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
# outm7med = function(p, mediators, control){
#   out = mediators %>% 
#     map_df(med_extact, f0(p), control) %>% 
#     filter(med!="NULL") %>%
#     unnest_longer(med) %>% 
#     hoist(med, p = "p") %>%
#     filter(p < threshold_med) %>% 
#     rename(p_med= p)
# }

outm72 = p_eqtl %>% map_df(outm7med, mediators, control, bespoke)

outm72 %>% 
  select(treatment, gene_set_name, med_id, mediator,p_med, p_eqtl) %>% 
  dplyr::filter(p_med < threshold_med) %>%
  mutate_at(.vars = vars(c("p_med")), .funs = funs(. %>%  str_remove("<") %>% as.numeric %>% format(scientific =TRUE))) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_med) %>%
  arrange(gene_set_name) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_med=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

