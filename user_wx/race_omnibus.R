#' ### treatment--raceethnicity: 
#' * "NonHwhite" as reference
#' * "NonHblack",
#' * "Hispanic"

focal = c("inflam1k_mRNA",
   "aging_down_cl1_mRNA",
   "aging_down_cl1a_mRNA",            
   "aging_down_cl1b_mRNA",
   "aging_down_cl1c_mRNA",
   "aging_down_cl2_mRNA",
   "aging_down_cl3_mRNA",
   "aging_up_cl1_mRNA",
   "aging_up_cl2_mRNA",
   "aging_up_cl3_mRNA", 
   "aging_up_cl4_mRNA" )
bespoke <- readRDS("/home/xu/ses-1/user_wx/race_bespoke_12.02.2021.rds")
p_eqtl <- c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
control = "ancestryPC_ses"
#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(tidyverse)
# functions to extract data
source("/home/xu/ses-1/user_wx/extract_utils.R")

#' ### minimum sig p from omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE

outomni_min = p_eqtl %>% map(m8_present, control)
outomni_min %>%
  bind_rows() %>%
  filter(p_omnibus<0.05) %>% 
  mutate(p_eqtl = p_eqtl %>% format(scientific =T)) %>%
  filter(gene_set_name %in% focal) %>% 
  arrange(gene_set_name) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### all sig p from omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE

outomni = p_eqtl %>% map(outm8_allsig, control)
temp = outomni %>%
  bind_rows() %>%
  filter(gene_sig %>% map_dfc(~ length(.))>0) %>%
  unnest_longer(gene_sig)

temp %>%
  hoist(out, p_omnibus = list("result", "m8_fdr", 1, "p")) %>%
  hoist(out, logFC = list("result", "m8_fdr", 1, "detail", "logFC")) %>% 
  select(-out, -ttT, -controls, - table1) %>% 
  arrange(gene_set_name) %>% 
  filter(gene_set_name %in% focal) %>% 
  mutate(p_eqtl = p_eqtl %>% format(scientific =T)) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### mediation for mean of all sig genes from omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
temp %>%
  hoist(out, med_mean=list("result", "m8_fdr", 1, "mediation_mean")) %>%
  unnest_longer(med_mean) %>%
  hoist(med_mean, med_mean_p = list("result", "p")) %>%
  filter(!is.na(med_mean_p)) %>%
  hoist(med_mean, med_mean_prop = list("result", "other", "med_prop")) %>%
  hoist(med_mean, med_mean_beta = list("result", "other", "med_beta")) %>%
  dplyr::select(-controls, -ttT, -out, -med_mean, -table1, - gene_sig) %>%
  distinct(.keep_all = TRUE) %>% 
  mutate(p_eqtl = p_eqtl %>% format(scientific =T)) %>%
  arrange(gene_set_name) %>% 
  filter(gene_set_name %in% focal) %>% 
  filter(med_mean_p<0.05) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()
