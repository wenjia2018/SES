
#' ---
#' title: race pattern with controlling for ses composite
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

#' ### Full control as in SES project: basic plus biological controls plus ses composite

#' ### two treatments:
#' * white_black (nonhispanic white and nonhispanic black) 
#' * white_hispanic(nonhispanic white and hispanic).
#' * Non hispanic white is the reference group
#' 
#' ###  For pca results, only the cases where the p value is less than 0.05/10 are presented



#' ## omnibus regression and self contained gene set test (mroast) results

#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(here)
library(tidyverse)
example_race_with1KI <- readRDS("~/ses-1/user_wx/example_race_with1KI_1030_withses4.rds")
race_ge_tfbm <- readRDS("~/ses-1/user_wx/example_race_detfbm_1030_withses4.rds")

# example_race_with1KI <- readRDS("~/ses-1/user_wx/example_race_with1KI_1030_withoutses4.rds")
# race_ge_tfbm <- readRDS("~/ses-1/user_wx/example_race_detfbm_withoutses4.rds")

race_ge_tfbm = race_ge_tfbm %>%
  filter(treatment=="raceethnicity") %>% 
  hoist(out, gsea = list("result", "gene_set_test")) %>% 
  select(1,4) %>%
  unnest(gsea) 

example_race_with1KI %>% 
  filter(treatment=="raceethnicity") %>% 
  hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
  filter(controls=="all") %>% 
  dplyr::select(treatment, gene_set_name, p) %>% 
  filter(p<0.05) %>%
  left_join(race_ge_tfbm, by = c("treatment", "gene_set_name")) %>% 
  rename(p_omnibus = p) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ## PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
threshold = 0.05/10
threshold_med = 0.05
var = example_race_with1KI %>%
  filter(treatment=="raceethnicity") %>% 
  hoist(out, var_explained = list("result", "m7_ob", 1, "other", "varexplained")) %>% 
  select(1,2,4)

var$var_explained = var$var_explained %>% map(~ set_names(.x, str_c("d", 1:10)))

var = var %>% unnest_longer(var_explained)


gene_list = example_race_with1KI %>%
  filter(treatment=="raceethnicity") %>% 
  hoist(out, well_loaded = list("result", "m7_ob", 1, "other", "well_loaded")) %>% 
  select(1,2,4)

gene_list$well_loaded = gene_list$well_loaded %>% map(~ set_names(.x, str_c("d", 1:10)))

gene_list = gene_list %>% unnest_longer(well_loaded)

example_race_with1KI %>%
  filter(treatment=="raceethnicity") %>% 
  hoist(out, p = list("result", "m7_ob", 1, "detail", "t")) %>%
  unnest_longer(p) %>%
  unnest(p) %>% 
  dplyr::select(1:9) %>% 
  dplyr::filter(p.value < threshold) %>%
  left_join(var, by = c("treatment", "gene_set_name", "p_id"= "var_explained_id")) %>% 
  left_join(gene_list, by = c("treatment", "gene_set_name", "p_id"= "well_loaded_id")) %>% 
  kableExtra::kable() %>%
  kableExtra::column_spec(column = 6, width = "16in", width_min="8in") %>% 
  kableExtra::kable_styling() 

#' ## BMI mediation
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example_race_with1KI %>%
  filter(treatment!="raceethnicity") %>%
  hoist(out, "result") %>%
  hoist(result, "m7_ob") %>% 
  unnest(matches("^m7")) %>% 
  hoist(m7_ob, result = list("mediation","w5bmi", "result")) %>% 
  filter(result!="NULL") %>% 
  unnest_longer(result) %>% 
  hoist(result,p = "p")%>% 
  select(1:4,6)  %>%
  filter(p < threshold_med) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ## stress_perceived mediation
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example_race_with1KI %>%
  filter(treatment!="raceethnicity") %>%
  hoist(out, "result") %>%
  hoist(result, "m7_ob") %>% 
  unnest(matches("^m7")) %>% 
  hoist(m7_ob, result = list("mediation","stress_perceived", "result")) %>% 
  filter(result!="NULL") %>% 
  unnest_longer(result) %>% 
  hoist(result,p = "p")%>% 
  select(1:4,6)  %>%
  filter(p < threshold_med) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ## bills mediation
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example_race_with1KI %>%
  filter(treatment!="raceethnicity") %>%
  hoist(out, "result") %>%
  hoist(result, "m7_ob") %>% 
  unnest(matches("^m7")) %>% 
  hoist(m7_ob, result = list("mediation","bills", "result")) %>% 
  filter(result!="NULL") %>% 
  unnest_longer(result) %>% 
  hoist(result,p = "p")%>% 
  select(1:4,6) %>%
  filter(p < threshold_med) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ## current smoke mediation
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example_race_with1KI %>%
  filter(treatment!="raceethnicity") %>%
  hoist(out, "result") %>%
  hoist(result, "m7_ob") %>% 
  unnest(matches("^m7")) %>% 
  hoist(m7_ob, result = list("mediation","currentsmoke", "result")) %>% 
  filter(result!="NULL") %>% 
  unnest_longer(result) %>% 
  hoist(result,p = "p")%>% 
  select(1:4,6) %>%
  filter(p < threshold_med) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ## insurance lack mediation
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example_race_with1KI %>%
  filter(treatment!="raceethnicity") %>%
  hoist(out, "result") %>%
  hoist(result, "m7_ob") %>% 
  unnest(matches("^m7")) %>% 
  hoist(m7_ob, result = list("mediation","insurance_lack", "result")) %>% 
  filter(result!="NULL") %>% 
  unnest_longer(result) %>% 
  hoist(result,p = "p")%>% 
  select(1:4,6) %>%
  filter(p < threshold_med) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()



#+ echo=F, eval=T, warning=FALSE, message=FALSE
# `rmarkdown::render("/home/xu/ses-1/user_wx/race_extraction_m7m8.R")`
