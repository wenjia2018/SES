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
  filter(controls == "ancestryPC") %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()
}

m7_present = function(race_ge_tfbm, example_race_with1KI){
threshold = 0.05/10
threshold_med = 0.05
var = example_race_with1KI %>%
  hoist(out, var_explained = list("result", "m7_ob", 1, "other", "varexplained")) %>% 
  dplyr::select(1,2,3,4)

var$var_explained = var$var_explained %>% map(~ set_names(.x, str_c("d", 1:10)))

var = var %>% unnest_longer(var_explained)


gene_list = example_race_with1KI %>%
  hoist(out, well_loaded = list("result", "m7_ob", 1, "other", "well_loaded")) %>% 
  dplyr::select(1,2,3,4)

gene_list$well_loaded = gene_list$well_loaded %>% map(~ set_names(.x, str_c("d", 1:10)))

gene_list = gene_list %>% unnest_longer(well_loaded)

example_race_with1KI %>%
  hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
  unnest_longer(p) %>%
  dplyr::select(treatment, gene_set_name, controls, p, p_id) %>% 
  dplyr::filter(p < threshold) %>%
  left_join(var, by = c("treatment", "gene_set_name", "controls", "p_id"= "var_explained_id")) %>% 
  left_join(gene_list, by = c("treatment", "gene_set_name","controls", "p_id"= "well_loaded_id")) %>% 
  filter(controls == "ancestryPC") %>% 
  dplyr::select(1:6) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling() 

}

#' ## p=0.05
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example_race_with1KI <- readRDS("~/ses-1/user_wx/example_race_dummy_26.11.rds")
race_ge_tfbm <- readRDS("~/ses-1/user_wx/race_detfbm_dummy_26.11.rds")
mroast_present(race_ge_tfbm)
m8_present(example_race_with1KI)
m7_present(race_ge_tfbm, example_race_with1KI)


#' ## p=1e-3
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example_race_with1KI <- readRDS("~/ses-1/user_wx/example_race_dummy1e-03.rds")
race_ge_tfbm <- readRDS("~/ses-1/user_wx/race_detfbm_dummy1e-03.rds")
mroast_present(race_ge_tfbm)
m8_present(example_race_with1KI)
m7_present(race_ge_tfbm, example_race_with1KI)



#' ## p=1e-4
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example_race_with1KI <- readRDS("~/ses-1/user_wx/example_race_dummy1e-04.rds")
race_ge_tfbm <- readRDS("~/ses-1/user_wx/race_detfbm_dummy1e-04.rds")
mroast_present(race_ge_tfbm)
m8_present(example_race_with1KI)
m7_present(race_ge_tfbm, example_race_with1KI)


#' ## p=1e-5
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example_race_with1KI <- readRDS("~/ses-1/user_wx/example_race_dummy1e-05.rds")
race_ge_tfbm <- readRDS("~/ses-1/user_wx/race_detfbm_dummy1e-05.rds")
mroast_present(race_ge_tfbm)
m8_present(example_race_with1KI)
m7_present(race_ge_tfbm, example_race_with1KI)


#' ## p=1e-6
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example_race_with1KI <- readRDS("~/ses-1/user_wx/example_race_dummy1e-06.rds")
race_ge_tfbm <- readRDS("~/ses-1/user_wx/race_detfbm_dummy1e-06.rds")
mroast_present(race_ge_tfbm)
m8_present(example_race_with1KI)
m7_present(race_ge_tfbm, example_race_with1KI)


#' ## p=1e-7
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example_race_with1KI <- readRDS("~/ses-1/user_wx/example_race_dummy1e-07.rds")
race_ge_tfbm <- readRDS("~/ses-1/user_wx/race_detfbm_dummy1e-07.rds")
mroast_present(race_ge_tfbm)
m8_present(example_race_with1KI)
m7_present(race_ge_tfbm, example_race_with1KI)


#' ## p=1e-8
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example_race_with1KI <- readRDS("~/ses-1/user_wx/example_race_dummy1e-08.rds")
race_ge_tfbm <- readRDS("~/ses-1/user_wx/race_detfbm_dummy1e-08.rds")
mroast_present(race_ge_tfbm)
m8_present(example_race_with1KI)
m7_present(race_ge_tfbm, example_race_with1KI)



#' ## p=1e-9
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example_race_with1KI <- readRDS("~/ses-1/user_wx/example_race_dummy1e-09.rds")
race_ge_tfbm <- readRDS("~/ses-1/user_wx/race_detfbm_dummy1e-09.rds")
mroast_present(race_ge_tfbm)
m8_present(example_race_with1KI)
m7_present(race_ge_tfbm, example_race_with1KI)

#' ## p=1e-10
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example_race_with1KI <- readRDS("~/ses-1/user_wx/example_race_dummy1e-10.rds")
race_ge_tfbm <- readRDS("~/ses-1/user_wx/race_detfbm_dummy1e-10.rds")
mroast_present(race_ge_tfbm)
m8_present(example_race_with1KI)
m7_present(race_ge_tfbm, example_race_with1KI)
