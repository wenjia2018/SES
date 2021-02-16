#' ### stratification--raceethnicity: 
#' * "NonHwhite"
#' * "NonHblack",
#' * "Hispanic"
#' 
#' 
#' ### treatments: 
#' * "color_byinterviewer3_White"
#' * "color_byinterviewer3_DarkBlack",
#' * "color_byinterviewer3_LightMed"



control = "ancestryPC_ses"
threshold = 0.05/10
threshold_med = 0.05
p_eqtl = c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
#' ### controls: 
#' 
#' * ancestryPC_ses : basic + ses + bespoke ancestryPC
#' 
#' ### outcome:
#'c(
#'  "ctra_mRNA",
#'  "inflame_mRNA",
#'  "interferon_mRNA",
#'  "AntBIntF_mRNA", 
#' "aging_mRNA",
#' "aging_up_mRNA",
#'  "aging_down_mRNA"
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


f0 = function(p){
  example0 <- bespoke %>%
    hoist(out, out1 = list("result","example0")) %>%
    select(1:3) %>%
    filter(p_eqtl == p) %>% 
    rename(out=out1) %>%
    unnest(out) 
  # %>% 
  # filter(control_set =="ancestryPC")
  return(example0)
}

f1 = function(p){
  example1 <- bespoke %>%
    hoist(out, out1 = list("result","example1")) %>%
    select(1:3) %>%
    filter(p_eqtl == p) %>% 
    rename(out=out1) %>% 
    unnest(out) 
  # %>% 
  #   filter(control_set =="ancestryPC") 
  return(example1)
}

mroast_present = function(example1){
  example1 %>%
    hoist(out, gsea = list("result", "gene_set_test")) %>% 
    dplyr::select(1,3,4) %>%
    unnest(gsea) %>% 
    as_tibble() %>% 
    filter(control_set == control) %>% 
    dplyr::filter(FDR.Mixed<0.05) %>% 
    kableExtra::kable() %>%
    kableExtra::kable_styling() 
}

m8_present = function(example0, control){
  
  out = example0 %>% 
    hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
    hoist(out, logFC = list("result", "m8_fdr", 1, "detail", "logFC")) %>% 
    dplyr::select(treatment, gene_set_name, p, logFC, control_set) %>% 
    dplyr::filter(p<0.05) %>%
    rename(p_omnibus = p) %>% 
    filter(control_set == control)
  return(out)
}

outm8 = function(p, control){
  m8_present(f0(p), control) %>% mutate(p_eqtl = p)
}

m7_present = function(example0, control){

  var = example0 %>%
    hoist(out, var_explained = list("result", "m7_ob", 1, "other", "varexplained")) %>%
    dplyr::select(treatment, gene_set_name, var_explained, control_set) %>% 
    filter(var_explained!="NULL")
  # filter out the non estimable results
  var$var_explained = var$var_explained %>% map(~ set_names(.x, str_c("d", 1:10)))

  var = var %>% unnest_longer(var_explained)
  
  
  gene_list = example0 %>%
    hoist(out, well_loaded = list("result", "m7_ob", 1, "other", "well_loaded")) %>%
    dplyr::select(treatment, gene_set_name, well_loaded, control_set) %>% 
    filter(well_loaded!="NULL")
  gene_list$well_loaded = gene_list$well_loaded %>% map(~ set_names(.x, str_c("d", 1:10)))

  gene_list = gene_list %>% unnest_longer(well_loaded)
  # 
  out = example0 %>%
    hoist(out, estimate = list("result", "m7_ob", 1, "detail", "t")) %>% 
    unnest_longer(estimate) %>% 
    hoist(estimate, p = "p.value") %>% 
    hoist(estimate, coef = "estimate") %>% 
    rename(p_id = estimate_id) %>% 
    dplyr::select(treatment, gene_set_name, p, coef, p_id, control_set) %>% 
    dplyr::filter(p < threshold) %>%
    left_join(var, by = c("treatment", "gene_set_name", "p_id"= "var_explained_id", "control_set")) %>%
    left_join(gene_list, by = c("treatment", "gene_set_name", "p_id"= "well_loaded_id", "control_set")) %>%
    filter(control_set == control) %>%
    dplyr::select(1:7) %>%
    rename(p_pca = p) 
  return(out)
  
}

outm7pca = function(p, control){
  m7_present(f0(p), control) %>% mutate(p_eqtl = p)
}

med_extact = function (focal, example0, control){
  out = example0 %>%
    hoist(out, med=list("result", "m7_ob", 1, "mediation", focal, "result")) %>% 
    mutate(mediator = focal) %>% 
    filter(control_set=="ancestryPC") %>% 
    select(-out)
}

outm7med = function(p, mediators, control){
  out = mediators %>% 
    map_df(med_extact, f0(p), control) %>% 
    filter(med!="NULL") %>%
    unnest_longer(med) %>% 
    hoist(med, p = "p") %>%
    filter(p < threshold_med) %>% 
    rename(p_med= p)
}


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

#' ## Strata: Nonhispanic White(white as reference)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE

bespoke <- readRDS("~/ses-1/user_wx/color_nonHwhite_strata_bespoke.rds")

outomni = p_eqtl %>% map(outm8,control)

outomni %>%
  bind_rows() %>%
  select(-logFC) %>% 
  pivot_wider(names_from = p_eqtl,values_from = p_omnibus) %>%
  arrange(treatment) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_omnibus=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()


#' ### omnibus regression coefficient
#+ echo=F, eval=T, warning=FALSE, message=FALSE

outomni %>%
  bind_rows() %>%
  select(-p_omnibus) %>% 
  pivot_wider(names_from = p_eqtl,values_from = logFC) %>%
  arrange(treatment) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = FALSE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("logFC=", .x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=", .x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()



#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE


outm71 = p_eqtl %>% map(outm7pca, control)


outm71 %>%
  bind_rows() %>%
  select(-coef) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_pca) %>%
  arrange(treatment) %>% 
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_pca=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### PCA regression coef
#+ echo=F, eval=T, warning=FALSE, message=FALSE
outm71 %>%
  bind_rows() %>%
  select(-p_pca) %>% 
  pivot_wider(names_from = p_eqtl, values_from = coef) %>%
  arrange(treatment) %>% 
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = FALSE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("coef=", .x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=", .x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### mediation
#+ echo=F, eval=T, warning=FALSE, message=FALSE

outm72 = p_eqtl %>% map_df(outm7med, mediators, control)

outm72 %>%
  select(treatment, gene_set_name, med_id, mediator,p_med, p_eqtl) %>%
  mutate_at(.vars = vars(c("p_med")), .funs = funs(. %>%  str_remove("<") %>% as.numeric %>% format(scientific =TRUE))) %>%
  pivot_wider(names_from = p_eqtl, values_from = p_med) %>%
  arrange(treatment) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_med=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()


#' ## Strata: Nonhispanic black(black as reference, level white has only one observation, could not get the estimation for this level)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
bespoke <- readRDS("~/ses-1/user_wx/color_NonHblack_strata_bespoke.rds")

#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE

outomni = p_eqtl %>% map(outm8,control)

outomni %>%
  bind_rows() %>%
  select(-logFC) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_omnibus) %>%
  arrange(treatment) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_omnibus=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()


#' ### omnibus regression coefficient
#+ echo=F, eval=T, warning=FALSE, message=FALSE

outomni %>%
  bind_rows() %>%
  select(-p_omnibus) %>% 
  pivot_wider(names_from = p_eqtl,values_from = logFC) %>%
  arrange(treatment) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = FALSE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("logFC=", .x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=", .x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()
#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE


outm71 = p_eqtl %>% map(outm7pca, control)


outm71 %>%
  bind_rows() %>%
  select(-coef) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_pca) %>%
  arrange(treatment) %>% 
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_pca=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### PCA regression coef
#+ echo=F, eval=T, warning=FALSE, message=FALSE
outm71 %>%
  bind_rows() %>%
  select(-p_pca) %>% 
  pivot_wider(names_from = p_eqtl, values_from = coef) %>%
  arrange(treatment) %>% 
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = FALSE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("coef=", .x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=", .x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()



#' ### mediation (not estimable)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

# outm72 = p_eqtl %>% map_df(outm7med, mediators, control)
#
# outm72 %>%
#   select(treatment, gene_set_name, med_id, mediator,p_med, p_eqtl) %>%
#   mutate_at(.vars = vars(c("p_med")), .funs = funs(. %>%  str_remove("<") %>% as.numeric %>% format(scientific =TRUE))) %>%
#   pivot_wider(names_from = p_eqtl, values_from = p_med) %>%
#   arrange(treatment) %>%
#   mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_med=",.x))) %>%
#   rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>%
#   kableExtra::kable() %>%
#   kableExtra::kable_styling()

#' ## Strata: Hispanic(white as reference, removing pregnant and tube control variables)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
bespoke <- readRDS("~/ses-1/user_wx/color_Hispanic_strata_bespoke_removetubepregnant.rds")

#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE

outomni = p_eqtl %>% map(outm8,control)

outomni %>%
  bind_rows() %>%
  select(-logFC) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_omnibus) %>%
  arrange(treatment) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_omnibus=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()


#' ### omnibus regression coefficient
#+ echo=F, eval=T, warning=FALSE, message=FALSE

outomni %>%
  bind_rows() %>%
  select(-p_omnibus) %>% 
  pivot_wider(names_from = p_eqtl,values_from = logFC) %>%
  arrange(treatment) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = FALSE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("logFC=", .x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=", .x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()
#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE


outm71 = p_eqtl %>% map(outm7pca, control)


outm71 %>%
  bind_rows() %>%
  select(-coef) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_pca) %>%
  arrange(treatment) %>% 
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_pca=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### PCA regression coef
#+ echo=F, eval=T, warning=FALSE, message=FALSE
outm71 %>%
  bind_rows() %>%
  select(-p_pca) %>% 
  pivot_wider(names_from = p_eqtl, values_from = coef) %>%
  arrange(treatment) %>% 
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = FALSE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("coef=", .x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=", .x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### mediation (not estimable)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

# outm72 = p_eqtl %>% map_df(outm7med, mediators, control)
# 
# outm72 %>% 
#   select(treatment, gene_set_name, med_id, mediator,p_med, p_eqtl) %>% 
#   mutate_at(.vars = vars(c("p_med")), .funs = funs(. %>%  str_remove("<") %>% as.numeric %>% format(scientific =TRUE))) %>% 
#   pivot_wider(names_from = p_eqtl, values_from = p_med) %>%
#   arrange(treatment) %>% 
#   mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_med=",.x))) %>%
#   rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
#   kableExtra::kable() %>%
#   kableExtra::kable_styling()
