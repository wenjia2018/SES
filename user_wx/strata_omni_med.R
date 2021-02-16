#' ### stratification--raceethnicity: 
#' * "NonHwhite"
#' * "NonHblack",
#' * "Hispanic"
#' 
#' 
#' ### treatments: 
#' * "color_byinterviewer3_White"
#' * "color_byinterviewer3_Light",
#' * "color_byinterviewer3_Medium"
#' * "color_byinterviewer3_Dark",
#' * "color_byinterviewer3_Black"


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

# outm8med = p_eqtl %>% map_df(outm8med, mediators, control)


#' ## Strata: Nonhispanic White(white as reference)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

bespoke = readRDS("~/ses-1/user_wx/color5_NonHwhite_strata_bespoke_withtubepregnant_m8_v2.rds")

outomni = p_eqtl %>% map(outm8_allsig, control)


temp = outomni %>%
  bind_rows() %>% 
  filter(gene_sig %>% map_dfc(~ length(.))>0) %>% 
  unnest_longer(gene_sig)

temp %>% 
  hoist(out, p_omnibus = list("result", "m8_fdr", 1, "p")) %>%
  hoist(out, logFC = list("result", "m8_fdr", 1, "detail", "logFC")) %>%  
  hoist(out, med=list("result", "m8_fdr", 1, "mediation")) %>% 
  unnest_longer(med) %>% 
  hoist(med, med_p = list("result", "p")) %>% 
  filter(!is.na(med_p)) %>% 
  hoist(med, med_prop = list("result", "other", "med_prop")) %>% 
  hoist(med, med_beta = list("result", "other", "med_beta")) %>% 
  select(-controls, -ttT, -out, -med, -table1) %>%
  mutate(p_eqtl = p_eqtl %>% format(scientific =T)) %>%
  kableExtra::kable() %>% 
  kableExtra::kable_styling()
  
#' ## Strata: Nonhispanic black(Medium as reference)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

bespoke = readRDS("~/ses-1/user_wx/color5_NonHblack_strata_bespoke_withtubepregnant_m8_v2.rds")

outomni = p_eqtl %>% map(outm8_allsig, control)


temp = outomni %>%
  bind_rows() %>% 
  filter(gene_sig %>% map_dfc(~ length(.))>0) %>% 
  unnest_longer(gene_sig)

temp %>% 
  hoist(out, p_omnibus = list("result", "m8_fdr", 1, "p")) %>%
  hoist(out, logFC = list("result", "m8_fdr", 1, "detail", "logFC")) %>%  
  hoist(out, med=list("result", "m8_fdr", 1, "mediation")) %>% 
  unnest_longer(med) %>% 
  hoist(med, med_p = list("result", "p")) %>% 
  filter(!is.na(med_p)) %>% 
  hoist(med, med_prop = list("result", "other", "med_prop")) %>% 
  hoist(med, med_beta = list("result", "other", "med_beta")) %>% 
  dplyr::select(-controls, -ttT, -out, -med, -table1) %>%
  mutate(p_eqtl = p_eqtl %>% format(scientific =T)) %>%
  # filter(med_p<0.05) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()


