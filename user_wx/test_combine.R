
#' ### treatments: black and hispanic

control = "ancestryPC"
threshold = 0.05/10
threshold_med = 0.05
#' ### controls: 3 sets controls: 
#' * basic : controls in ses paper
#' * ses : basic + ses
#' * ancestryPC : basic + ses + bespoke ancestryPC
#' 
#' 
#' ### bespoke ancestryPC:
#' 
#' * eqtl ancestry: p_eqtl is a sequence of p (0.05, 0.01, 1e-3,...,1e-10) used to choose asscoiated snps
#' * spatial ancestry: 2000 upstream and down stream region of gene in our outcome set
#' 
#' ### outcome:
#'c(
#'  "ctra_mRNA",
#'  "inflame_mRNA",
#'  "interferon_mRNA",
#'  "AntBIntF_mRNA", 
#'  "inflam1k_mRNA",
#' "aging_mRNA",
#' "aging_up_mRNA",
#'  "aging_down_mRNA"
#')
#'
#' ### analysis: eqtl ancestry and spatial ancestry 
#' * omnibus P value: association between each disease set and color related snps, 
#' the smallest, whole-genome FDR corrected p-value via whole genome regressions
#'
#' * PCA p value: p value associated between each PC of the disease set and color related snps
#' 
#' 
#' * omnibus p values are genowide corrected p values
#' 
#' * pca regression p value and mediation p values are unadjusted 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(here)
library(tidyverse)
library(naniar)


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
    dplyr::select(treatment, gene_set_name, p, control_set) %>% 
    dplyr::filter(p<0.05) %>%
    rename(p_omnibus = p) %>% 
    filter(control_set == control)
  return(out)
}

m7_present = function(example0, control){
  threshold = 0.05/10
  var = example0 %>%
    hoist(out, var_explained = list("result", "m7_ob", 1, "other", "varexplained")) %>% 
    dplyr::select(treatment, gene_set_name, var_explained, control_set)
  
  var$var_explained = var$var_explained %>% map(~ set_names(.x, str_c("d", 1:10)))
  
  var = var %>% unnest_longer(var_explained)
  
  
  gene_list = example0 %>%
    hoist(out, well_loaded = list("result", "m7_ob", 1, "other", "well_loaded")) %>% 
    dplyr::select(treatment, gene_set_name, well_loaded, control_set)
  
  gene_list$well_loaded = gene_list$well_loaded %>% map(~ set_names(.x, str_c("d", 1:10)))
  
  gene_list = gene_list %>% unnest_longer(well_loaded)
  
  out = example0 %>%
    hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
    unnest_longer(p) %>%
    dplyr::select(treatment, gene_set_name, p, p_id, control_set) %>% 
    dplyr::filter(p < threshold) %>%
    left_join(var, by = c("treatment", "gene_set_name", "p_id"= "var_explained_id", "control_set")) %>% 
    left_join(gene_list, by = c("treatment", "gene_set_name", "p_id"= "well_loaded_id", "control_set")) %>% 
    filter(control_set == control) %>%
    dplyr::select(1:6) %>% 
    rename(p_pca = p) 
  return(out)
  
}


med_extact = function (focal, example0, control){
  out = example0 %>%
    hoist(out, med=list("result", "m7_ob", 1, "mediation", focal, "result")) %>% 
    mutate(mediator = focal) %>% 
    filter(control_set=="ancestryPC") %>% 
    select(-out)
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
    
    "totdiscrim1_category", #ordered logistic fit glm for mediation
    "totdiscrim2_category", #ordered logistic fit glm for mediation
    "countdiscrimwhy_category"#ordered logistic fit glm for mediation
  )



#+ echo=F, eval=T, warning=FALSE, message=FALSE
bespoke_g2s <- readRDS("~/ses-1/user_wx/bespoke_g2s.rds")

example0_g2s <- bespoke_g2s %>%
  hoist(out, out1 = list("result","example0")) %>%
  select(1:2) %>%
  rename(out=out1) %>%
  unnest(out) 

g2s_m8 = m8_present(example0_g2s, control) %>%
  mutate_at(vars("p_omnibus"), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars("p_omnibus"), .funs = list( ~ str_c("p_omnibus=",.x))) %>%
  rename_at(vars("p_omnibus"), .funs = list( ~ "spatial_2000")) 

#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE

g2s_m71 = m7_present(example0_g2s, control) %>%
  mutate_at(vars("p_pca"), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars("p_pca"), .funs = list( ~ str_c("p_pca=",.x))) %>%
  rename_at(vars("p_pca"), .funs = list( ~ "spatial_2000")) 

#' ### mediation 
#+ echo=F, eval=T, warning=FALSE, message=FALSE

g2s_m72 = mediators %>% 
  map_df(med_extact, example0, control) %>% 
  filter(med!="NULL") %>%
  unnest_longer(med) %>% 
  hoist(med, p = "p") %>%
  filter(p < threshold_med) %>% 
  rename(p_med= p) %>% 
  select(2,3,5,7,8,9) %>% 
  mutate_at(vars("p_med"), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars("p_med"), .funs = list( ~ str_c("p_med=",.x))) %>%
  rename_at(vars("p_med"), .funs = list( ~ "spatial_2000")) 





#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
bespoke <- readRDS("~/ses-1/user_wx/bespoke_v4.rds") 
p_eqtl = c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
outm8 = function(p, control){
  m8_present(f0(p), control) %>% mutate(p_eqtl = p)
}

outomni = p_eqtl %>% map(outm8,control)

outomni %>%
  bind_rows() %>%
  pivot_wider(names_from = p_eqtl, values_from = p_omnibus) %>%
  arrange(treatment) %>%
  mutate_at(vars(c(4:13)), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(c(4:13)), .funs = list( ~ str_c("p_omnibus=",.x))) %>%
  rename_at(vars(c(4:13)), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  left_join(g2s_m8) %>% 
  relocate("spatial_2000", .after = "control_set") %>% 
  replace_with_na_all(condition = ~ .x=="p_omnibus=     NA") %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()
#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE

outm7pca = function(p, control){
  m7_present(f0(p), control) %>% mutate(p_eqtl = p)
}

outm71 = p_eqtl %>% map(outm7pca,control)


outm71 %>%
  bind_rows() %>%
  pivot_wider(names_from = p_eqtl, values_from = p_pca) %>%
  arrange(treatment) %>% 
  mutate_at(vars(c(6:15)), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(c(6:15)), .funs = list( ~ str_c("p_pca=",.x))) %>%
  rename_at(vars(c(6:15)), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  left_join(g2s_m71) %>%
  relocate("spatial_2000", .after = "var_explained") %>% 
  replace_with_na_all(condition = ~ .x=="p_pca=     NA") %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()
#' ### mediation 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
outm7med = function(p, mediators, control){
  out = mediators %>% 
    map_df(med_extact, f0(p), control) %>% 
    filter(med!="NULL") %>%
    unnest_longer(med) %>% 
    hoist(med, p = "p") %>%
    filter(p < threshold_med) %>% 
    rename(p_med= p)
}

outm72 = p_eqtl %>% map_df(outm7med, mediators, control)

outm72 %>% 
  select(treatment, gene_set_name, med_id, mediator,p_med, p_eqtl) %>% 
  mutate_at(.vars = vars(c("p_med")), .funs = funs(. %>%  str_remove("<") %>% as.numeric %>% format(scientific =TRUE))) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_med) %>%
  arrange(treatment) %>% 
  mutate_at(vars(c(5:14)), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(c(5:14)), .funs = list( ~ str_c("p_med=",.x))) %>%
  rename_at(vars(c(5:14)), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  left_join(g2s_m72) %>% 
  relocate("spatial_2000", .after = "mediator") %>% 
  replace_with_na_all(condition = ~ str_detect(.x, "=NA")) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()



