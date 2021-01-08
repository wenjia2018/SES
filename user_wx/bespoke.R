# for first version when multiple results are generated, need to filter out results by gene_name
#' ### treatments: black and hispanic

#' ### controls: Bespoke ancestryPC
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
#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(here)
library(tidyverse)

mroast_present = function(race_ge_tfbm){
  race_ge_tfbm %>%
    hoist(out, gsea = list("result", "gene_set_test")) %>% 
    dplyr::select(1,3,4) %>%
    unnest(gsea) %>% 
    as_tibble() %>% 
    filter(controls == "ancestryPC") %>% 
    dplyr::filter(FDR.Mixed<0.05) %>% 
    kableExtra::kable() %>%
    kableExtra::kable_styling() 
}

m8_present = function(example_race_with1KI){
  
  example_race_with1KI %>% 
    hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
    dplyr::select(treatment, gene_set_name, controls, p) %>% 
    dplyr::filter(p<0.05) %>%
    rename(p_omnibus = p) %>% 
    # filter(controls == "ancestryPC") %>% 
    kableExtra::kable() %>%
    kableExtra::kable_styling()
}

m7_present = function(race_ge_tfbm, example_race_with1KI){
  threshold = 0.05/10
  threshold_med = 0.05
  var = example_race_with1KI %>%
    hoist(out, var_explained = list("result", "m7_ob", 1, "other", "varexplained")) %>% 
    dplyr::select(1,2,3,5)
  
  var$var_explained = var$var_explained %>% map(~ set_names(.x, str_c("d", 1:10)))
  
  var = var %>% unnest_longer(var_explained)
  
  
  gene_list = example_race_with1KI %>%
    hoist(out, well_loaded = list("result", "m7_ob", 1, "other", "well_loaded")) %>% 
    dplyr::select(1,2,3,5)
  
  gene_list$well_loaded = gene_list$well_loaded %>% map(~ set_names(.x, str_c("d", 1:10)))
  
  gene_list = gene_list %>% unnest_longer(well_loaded)
  
  example_race_with1KI %>%
    hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
    unnest_longer(p) %>%
    dplyr::select(treatment, gene_set_name, controls, p, p_id) %>% 
    dplyr::filter(p < threshold) %>%
    left_join(var, by = c("treatment", "gene_set_name", "p_id"= "var_explained_id")) %>% 
    left_join(gene_list, by = c("treatment", "gene_set_name", "p_id"= "well_loaded_id")) %>% 
    # filter(controls == "ancestryPC") %>% 
    dplyr::select(1:6) %>% 
    kableExtra::kable() %>%
    kableExtra::kable_styling() 
  
}

#' ## p=0.05
#+ echo=F, eval=T, warning=FALSE, message=FALSE
bespoke <- readRDS("~/ses-1/user_wx/bespoke_all.rds") 
example_race_with1KI <- bespoke %>%
  hoist(out, out1 = list("result","example0")) %>%
  select(1:3) %>%
  rename(out=out1, gene = gene_set_name) %>% 
  unnest(out) %>% 
  filter(gene ==gene_set_name) %>% 
  filter(p_eqtl == 0.05, control_set =="ancestryPC") %>% 
  select(-gene)

race_ge_tfbm <- bespoke %>%
  hoist(out, out1 = list("result","example1")) %>%
  select(1:3) %>%
  rename(out=out1, gene = gene_set_name) %>% 
  unnest(out) %>% 
  filter(p_eqtl == 0.05, control_set =="ancestryPC") %>% 
  unique()

# mroast_present(race_ge_tfbm)
#' ## omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m8_present(example_race_with1KI)

#' ## PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m7_present(race_ge_tfbm, example_race_with1KI)
