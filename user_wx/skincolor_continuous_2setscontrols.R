
#' ---
#' title: skin color as a continuous variable with basic control, and race ethnicity control
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE


#' ## Controls, treatment and threshold
#' 
#' ### 2 sets of controls:
#' 
#' * basic control as in SES project
#' * race control: basic + raceenthnicity_NonHblack + raceethnicity_Hispanic
#' 
#' ### one focal treatment:
#' 
#' * white is 0, black is 5
#' 
#' ###  For pca results, only the cases where the p value is less than 0.05/10 are presented
#' 
#' ## omnibus regression and self contained gene set test (mroast) results

#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(here)
library(tidyverse)
example_race_with1KI <- readRDS("~/ses-1/user_wx/example_skincolor_continuous_13.11.rds")
race_ge_tfbm <- readRDS("~/ses-1/user_wx/skincolor_detfbm_continuous_13.11.rds")

# example_race_with1KI <- readRDS("~/ses-1/user_wx/example_race_with1KI_1030_withoutses4.rds")
# race_ge_tfbm <- readRDS("~/ses-1/user_wx/example_race_detfbm_withoutses4.rds")

race_ge_tfbm = race_ge_tfbm %>%
  hoist(out, gsea = list("result", "gene_set_test")) %>% 
  dplyr::select(1,3,4) %>%
  unnest(gsea) 

example_race_with1KI %>% 
  hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
  dplyr::select(treatment, gene_set_name, controls, p) %>% 
  dplyr::filter(p<0.05) %>%
  left_join(race_ge_tfbm, by = c("treatment", "gene_set_name", "controls")) %>% 
  rename(p_omnibus = p) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ## PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
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
  dplyr::select(1:6) %>% 
  kableExtra::kable() %>%
  # kableExtra::column_spec(column = 7, width = "16in", width_min="8in") %>% 
  kableExtra::kable_styling() 

#' ## only p <0.05 are presentd for mediation analysis 

#' ### BMI mediation
#+ echo=F, eval=T, warning=FALSE, message=FALSE

med_present = function (focal, data){
  out =  data %>%
    hoist(out, "result") %>%
    hoist(result, "m7_ob") %>% 
    unnest(matches("^m7")) %>% 
    hoist(m7_ob, result = list("mediation", focal , "result")) %>% 
    filter(result!="NULL") %>% 
    unnest_longer(result) %>% 
    hoist(result,p = "p")%>% 
    dplyr::select(1:4,6)  %>%
    filter(p < threshold_med)
  return(out)
}

med_present(focal = "w5bmi_lm", data = example_race_with1KI)  %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### stress_perceived mediation
#+ echo=F, eval=T, warning=FALSE, message=FALSE

med_present(focal = "stress_perceived_lm", data = example_race_with1KI)  %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### bills mediation
#+ echo=F, eval=T, warning=FALSE, message=FALSE

med_present(focal = "bills_binary", data = example_race_with1KI)  %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### current smoke mediation
#+ echo=F, eval=T, warning=FALSE, message=FALSE

med_present(focal = "currentsmoke_binary", data = example_race_with1KI)  %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### insurance lack mediation
#+ echo=F, eval=T, warning=FALSE, message=FALSE

med_present(focal = "insurance_lack_binary", data = example_race_with1KI)  %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### lowbirthweight mediation
#+ echo=F, eval=T, warning=FALSE, message=FALSE

med_present(focal = "lowbirthweight_binary", data = example_race_with1KI)  %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### "high_lowbirth_binary"  mediation
#+ echo=F, eval=T, warning=FALSE, message=FALSE

med_present(focal = "high_lowbirth_binary" , data = example_race_with1KI)  %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()
#' ### "totdiscrim2_binary" mediation
#+ echo=F, eval=T, warning=FALSE, message=FALSE

med_present(focal = "totdiscrim2_binary", data = example_race_with1KI)  %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()
#' ### "discrim2_binary" mediation
#+ echo=F, eval=T, warning=FALSE, message=FALSE

med_present(focal = "discrim2_binary", data = example_race_with1KI)  %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()
#' ### "totdiscrim1_category" mediation
#+ echo=F, eval=T, warning=FALSE, message=FALSE

med_present(focal = "totdiscrim1_category", data = example_race_with1KI)  %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

