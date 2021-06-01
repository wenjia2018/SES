
#' ---
#' title: skin color 3 within Race strata with ancestryPC control
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

#' ### treatments: 
#' * "color_byinterviewer3_White",
#' * "color_byinterviewer3_LightMed"
#' * "color_byinterviewer3_DarkBlack"

#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(tidyverse)
control = "ancestryPC_ses"
threshold = 0.05/10
threshold_med = 0.05
p_eqtl = c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)


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
source("/home/xu/ses-1/user_wx/extract_v2.R")

#' ## Strata: Nonhispanic White(white as reference)
#+ echo=F, eval=T, warning=FALSE, message=FALSE


# bespoke <- readRDS("/home/xu/ses-1/user_wx/color3_bespoke_NonHwhite_strata_25.02.2021.rds")
bespoke <- readRDS("/home/xu/ses-1/user_wx/color3_bespoke_NonHwhite_strata_16.03.2021.rds")
#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
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



#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE


outm71 = p_eqtl %>% map(outm7pca, control, bespoke)


outm71 %>%
  bind_rows() %>%
  dplyr::filter(p_pca < threshold) %>%
  select(-coef, -well_loaded, -loadings) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_pca) %>%
  arrange(gene_set_name) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_pca=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### PCA regression coef
#+ echo=F, eval=T, warning=FALSE, message=FALSE
outm71 %>%
  bind_rows() %>%
  dplyr::filter(p_pca < threshold) %>%
  select(-p_pca, -well_loaded, -loadings) %>% 
  pivot_wider(names_from = p_eqtl, values_from = coef) %>%
  arrange(treatment) %>% 
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = FALSE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("coef=", .x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=", .x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()


#' ### PCA mediation 
#+ echo=F, eval=T, warning=FALSE, message=FALSE

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

#' ## Strata: Nonhispanic black(darkblack as reference)
#+ echo=F, eval=T, warning=FALSE, message=FALSE


# bespoke <- readRDS("/home/xu/ses-1/user_wx/color3_bespoke_NonHblack_strata_25.02.2021.rds")
bespoke <- readRDS("/home/xu/ses-1/user_wx/color3_bespoke_NonHblack_strata_16.03.2021.rds")
outomni = p_eqtl %>% map(m8_present, control, bespoke)
#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE

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



#' ### omnibus regression coefficient
#+ echo=F, eval=T, warning=FALSE, message=FALSE

outomni %>%
  bind_rows() %>%
  dplyr::filter(p_omnibus < 0.05) %>%
  select(-p_omnibus) %>% 
  pivot_wider(names_from = p_eqtl,values_from = logFC) %>%
  arrange(treatment) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = FALSE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("logFC=", .x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=", .x))) %>% 
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


outm71 = p_eqtl %>% map(outm7pca, control, bespoke)


outm71 %>%
  bind_rows() %>%
  dplyr::filter(p_pca < threshold) %>%
  select(-coef, -well_loaded, -loadings) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_pca) %>%
  arrange(gene_set_name) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_pca=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### PCA regression coef
#+ echo=F, eval=T, warning=FALSE, message=FALSE
outm71 %>%
  bind_rows() %>%
  dplyr::filter(p_pca < threshold) %>%
  select(-p_pca, -well_loaded, -loadings) %>% 
  pivot_wider(names_from = p_eqtl, values_from = coef) %>%
  arrange(treatment) %>% 
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = FALSE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("coef=", .x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=", .x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### PCA mediation 
#+ echo=F, eval=T, warning=FALSE, message=FALSE

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

#' ## Strata: Hispanic (white as reference)
#+ echo=F, eval=T, warning=FALSE, message=FALSE


# bespoke <- readRDS("/home/xu/ses-1/user_wx/color3_bespoke_Hispanic_strata_25.02.2021.rds")
bespoke <- readRDS("/home/xu/ses-1/user_wx/color3_bespoke_Hispanic_strata_16.03.2021.rds")
#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
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



#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE


outm71 = p_eqtl %>% map(outm7pca, control, bespoke)


outm71 %>%
  bind_rows() %>%
  dplyr::filter(p_pca < threshold) %>%
  select(-coef, -well_loaded, -loadings) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_pca) %>%
  arrange(gene_set_name) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_pca=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()


#' ### PCA regression coef
#+ echo=F, eval=T, warning=FALSE, message=FALSE
outm71 %>%
  bind_rows() %>%
  dplyr::filter(p_pca < threshold) %>%
  select(-p_pca, -well_loaded, -loadings) %>% 
  pivot_wider(names_from = p_eqtl, values_from = coef) %>%
  arrange(treatment) %>% 
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = FALSE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("coef=", .x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=", .x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### PCA mediation 
#+ echo=F, eval=T, warning=FALSE, message=FALSE

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

