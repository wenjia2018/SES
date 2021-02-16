bespoke = readRDS("~/ses-1/user_wx/xxx.rds")
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
#' * omnibus P value: association between each disease set and via whole genome regressions
#' * regresssion results of whole-genome FDR corrected p-value <0.05 
#' * mediational results for the adj.p<0.05 genes.

#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(tidyverse)
f0 = function(p){
  example0 <- bespoke %>%
    hoist(out, out1 = list("result","example0")) %>%
    dplyr::select(1:3) %>%
    filter(p_eqtl == p) %>% 
    dplyr::rename(out=out1) %>%
    unnest(out) 
  # %>% 
  # filter(control_set =="ancestryPC")
  return(example0)
}


m8_present = function(example0, control){
  
  out = example0 %>% 
    hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
    hoist(out, logFC = list("result", "m8_fdr", 1, "detail", "logFC")) %>% 
    dplyr::select(treatment, gene_set_name, p, logFC, control_set) %>% 
    # dplyr::filter(p<0.05) %>%
    dplyr::rename(p_omnibus = p) %>% 
    filter(control_set == control)
  return(out)
}
outm8_allsig = function(p, control){
  out = f0(p) %>% 
    hoist(out, ttT = list("result", "m8_fdr", 1, "other", "m")) %>% 
    filter(ttT != "NULL") %>% 
    mutate(gene_sig = map(ttT, ~ dplyr::filter(., adj.P.Val<0.05) %>% pull(gene))) %>% 
    # select(-out) %>% 
    filter(control_set == control) %>%
    mutate(p_eqtl = p)
}

med_extact = function (focal, example0, control){
  out = example0 %>%
    hoist(out, med=list("result", "m8_fdr", 1, "mediation")) %>% 
    mutate(mediator = focal) %>% 
    filter(control_set== control) %>% 
    dplyr::select(-out)
}

outm8med = function(p, mediators, control){
  out = mediators %>% 
    map_df(med_extact, f0(p), control) %>% 
    filter(med!="NULL") %>%
    unnest_longer(med) %>% 
    hoist(med, p = list("result","p")) %>%
    filter(!is.na(p)) %>% 
    # filter(p < threshold_med) %>% 
    dplyr::rename(p_med= p)
}




outomni = p_eqtl %>% map(outm8_allsig, control =control)

temp = 
  example %>% 
  hoist(out, ttT = list("result", "m8_fdr", 1, "other", "m")) %>% 
  filter(ttT != "NULL") %>% 
  mutate(gene_sig = map(ttT, ~ dplyr::filter(., adj.P.Val<0.05) %>% pull(gene))) %>% 
  filter(control_set == control) %>%
  filter(gene_sig %>% map_dfc(~ length(.))>0) %>% 
  unnest_longer(gene_sig)
  
  
  temp %>% 
  hoist(out, p_omnibus = list("result", "m8_fdr", 1, "p")) %>%
  hoist(out, logFC = list("result", "m8_fdr", 1, "detail", "logFC")) %>%  
  hoist(out, med_single=list("result", "m8_fdr", 1, "mediation_signle")) %>% 
  unnest_longer(med_single) %>% 
  hoist(med_single, med_single_p = list("result", "p")) %>% 
  filter(!is.na(med_single_p)) %>% 
  hoist(med_single, med_single_prop = list("result", "other", "med_prop")) %>% 
  hoist(med_single, med_single_beta = list("result", "other", "med_beta")) %>% 
  dplyr::select(-controls, -ttT, -out, -med_single, -table1) %>%
  mutate(p_eqtl = p_eqtl %>% format(scientific =T)) %>%
  filter(med_single_p<0.05) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()

  temp %>% 
    hoist(out, p_omnibus = list("result", "m8_fdr", 1, "p")) %>%
    hoist(out, logFC = list("result", "m8_fdr", 1, "detail", "logFC")) %>%  
    hoist(out, med_single=list("result", "m8_fdr", 1, "mediation_single")) %>% 
    unnest_longer(med_single) %>% 
    hoist(med_single, med_single_p = list("result", "p")) %>% 
    filter(!is.na(med_single_p)) %>% 
    hoist(med_single, med_single_prop = list("result", "other", "med_prop")) %>% 
    hoist(med_single, med_single_beta = list("result", "other", "med_beta")) %>% 
    dplyr::select(-controls, -ttT, -out, -med_single, -table1) %>%
    mutate(p_eqtl = p_eqtl %>% format(scientific =T)) %>%
    filter(med_single_p<0.05) %>% 
    kableExtra::kable() %>% 
    kableExtra::kable_styling()
  
  
  
temp %>% 
  hoist(out, med_mean=list("result", "m8_fdr", 1, "mediation_mean")) %>% 
  unnest_longer(med_mean) %>% 
  hoist(med_mean, med_mean_p = list("result", "p")) %>% 
  filter(!is.na(med_mean_p)) %>% 
  hoist(med_mean, med_mean_prop = list("result", "other", "med_prop")) %>% 
  hoist(med_mean, med_mean_beta = list("result", "other", "med_beta")) %>% 
  dplyr::select(-controls, -ttT, -out, -med_mean, -table1) %>%
  mutate(p_eqtl = p_eqtl %>% format(scientific =T)) %>%
  filter(med_mean_p<0.05) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()


temp %>% 
  hoist(out, med=list("result", "m8_fdr", 1, "mediation")) %>% 
  unnest_longer(med) %>% 
  hoist(med, med_p = list("result", "p")) %>% 
  filter(!is.na(med_p)) %>% 
  hoist(med, med_prop = list("result", "other", "med_prop")) %>% 
  hoist(med, med_beta = list("result", "other", "med_beta")) %>% 
  dplyr::select(-controls, -ttT, -out, -med, -table1) %>%
  mutate(p_eqtl = p_eqtl %>% format(scientific =T)) %>%
  filter(med_p<0.05) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()

