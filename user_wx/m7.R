
#' ---
#' title: oblique rotation PCA mediational analysis
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE


#' We do PCA for each disease signature using 3 rotations and then do mediational analysis
#'  for significant PCs using w5bmi as mediator
#'  
#' p: p value for indirected effect / mediated effect;
#' result_id: which pc is significant 
library(here)
library(tidyverse)
example0_with1k = readRDS("/home/share/scratch/xu/example0_new2signature_withinflam.rds")

#' ## keep inflamation signatures
#' 
threshold_with1k = 0.05/11/4
threshold_med = 0.05/11/5
fig1panelB <- example0_with1k %>%
  hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
  unnest_longer(p) %>% 
  dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  dplyr::filter(p < threshold_with1k) %>% 
  dplyr::select(treatment, gene_set_name, p_id) %>% 
  mutate(id = str_c(treatment,"_",gene_set_name,"_",p_id))

#' ### bmi at wave 5

#+ echo=F, eval=T, warning=FALSE, message=FALSE
w5bmi = example0_with1k %>%
  hoist(out, "result") %>%
  hoist(result, "m7_ob") %>% 
  unnest(matches("^m7")) %>% 
  hoist(m7_ob, result = list("mediation","w5bmi", "result")) %>% 
  filter(result!="NULL") %>% 
  unnest_longer(result) %>% 
  hoist(result,p = "p")%>% 
  select(1:4,6) %>% 
  filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  mutate(id = str_c(treatment,"_",gene_set_name,"_",result_id)) %>% 
  right_join(fig1panelB %>% select(id))

w5bmi %>% kableExtra::kable() %>% kableExtra::kable_styling()

#' ### perceived stress

#+ echo=F, eval=T, warning=FALSE, message=FALSE
stress_perceived = example0_with1k %>%
  hoist(out, "result") %>%
  hoist(result, "m7_ob") %>% 
  unnest(matches("^m7")) %>% 
  hoist(m7_ob, result = list("mediation","stress_perceived", "result")) %>% 
  filter(result!="NULL") %>% 
  unnest_longer(result) %>% 
  hoist(result,p = "p")%>% 
  select(1:4,6) %>% 
  filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  mutate(id = str_c(treatment,"_",gene_set_name,"_",result_id)) %>% 
  right_join(fig1panelB %>% select(id))

stress_perceived %>% kableExtra::kable() %>% kableExtra::kable_styling()

#' ### late for bills

#+ echo=F, eval=T, warning=FALSE, message=FALSE
bills = example0_with1k %>%
  hoist(out, "result") %>%
  hoist(result, "m7_ob") %>% 
  unnest(matches("^m7")) %>% 
  hoist(m7_ob, result = list("mediation","bills", "result")) %>% 
  filter(result!="NULL") %>% 
  unnest_longer(result) %>% 
  hoist(result,p = "p")%>% 
  select(1:4,6) %>% 
  filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  mutate(id = str_c(treatment,"_",gene_set_name,"_",result_id)) %>% 
  right_join(fig1panelB %>% select(id))

bills %>% kableExtra::kable() %>% kableExtra::kable_styling()

#' ### current smoke

#+ echo=F, eval=T, warning=FALSE, message=FALSE
currentsmoke = example0_with1k %>%
  hoist(out, "result") %>%
  hoist(result, "m7_ob") %>% 
  unnest(matches("^m7")) %>% 
  hoist(m7_ob, result = list("mediation","currentsmoke", "result")) %>% 
  filter(result!="NULL") %>% 
  unnest_longer(result) %>% 
  hoist(result,p = "p")%>% 
  select(1:4,6) %>% 
  filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  mutate(id = str_c(treatment,"_",gene_set_name,"_",result_id)) %>% 
  right_join(fig1panelB %>% select(id))

currentsmoke %>% kableExtra::kable() %>% kableExtra::kable_styling()

#' ### no insurance because of no money

#+ echo=F, eval=T, warning=FALSE, message=FALSE
insurance_lack = example0_with1k %>%
  hoist(out, "result") %>%
  hoist(result, "m7_ob") %>% 
  unnest(matches("^m7")) %>% 
  hoist(m7_ob, result = list("mediation","insurance_lack", "result")) %>% 
  filter(result!="NULL") %>% 
  unnest_longer(result) %>% 
  hoist(result,p = "p")%>% 
  select(1:4,6) %>% 
  filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  mutate(id = str_c(treatment,"_",gene_set_name,"_",result_id)) %>% 
  right_join(fig1panelB %>% select(id))

insurance_lack %>% kableExtra::kable() %>% kableExtra::kable_styling()

#' ## remove inflamation signatures


example0_no1k = readRDS("/home/share/scratch/xu/example0_new2signature_noinflam.rds")

threshold_with1k = 0.05/11/4
threshold_med = 0.05/11/5
fig1panelB <- example0_no1k %>%
  hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
  unnest_longer(p) %>% 
  dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  dplyr::filter(p < threshold_with1k) %>% 
  dplyr::select(treatment, gene_set_name, p_id) %>% 
  mutate(id = str_c(treatment,"_",gene_set_name,"_",p_id))

#' ### bmi at wave 5

#+ echo=F, eval=T, warning=FALSE, message=FALSE
w5bmi = example0_no1k %>%
  hoist(out, "result") %>%
  hoist(result, "m7_ob") %>% 
  unnest(matches("^m7")) %>% 
  hoist(m7_ob, result = list("mediation","w5bmi", "result")) %>% 
  filter(result!="NULL") %>% 
  unnest_longer(result) %>% 
  hoist(result,p = "p")%>% 
  select(1:4,6) %>% 
  filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  mutate(id = str_c(treatment,"_",gene_set_name,"_",result_id)) %>% 
  right_join(fig1panelB %>% select(id))

w5bmi %>% kableExtra::kable() %>% kableExtra::kable_styling()

#' ### perceived stress

#+ echo=F, eval=T, warning=FALSE, message=FALSE
stress_perceived = example0_no1k %>%
  hoist(out, "result") %>%
  hoist(result, "m7_ob") %>% 
  unnest(matches("^m7")) %>% 
  hoist(m7_ob, result = list("mediation","stress_perceived", "result")) %>% 
  filter(result!="NULL") %>% 
  unnest_longer(result) %>% 
  hoist(result,p = "p")%>% 
  select(1:4,6) %>% 
  filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  mutate(id = str_c(treatment,"_",gene_set_name,"_",result_id)) %>% 
  right_join(fig1panelB %>% select(id))

stress_perceived %>% kableExtra::kable() %>% kableExtra::kable_styling()

#' ### late for bills

#+ echo=F, eval=T, warning=FALSE, message=FALSE
bills = example0_no1k %>%
  hoist(out, "result") %>%
  hoist(result, "m7_ob") %>% 
  unnest(matches("^m7")) %>% 
  hoist(m7_ob, result = list("mediation","bills", "result")) %>% 
  filter(result!="NULL") %>% 
  unnest_longer(result) %>% 
  hoist(result,p = "p")%>% 
  select(1:4,6) %>% 
  filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  mutate(id = str_c(treatment,"_",gene_set_name,"_",result_id)) %>% 
  right_join(fig1panelB %>% select(id))

bills %>% kableExtra::kable() %>% kableExtra::kable_styling()

#' ### current smoke

#+ echo=F, eval=T, warning=FALSE, message=FALSE
currentsmoke = example0_no1k %>%
  hoist(out, "result") %>%
  hoist(result, "m7_ob") %>% 
  unnest(matches("^m7")) %>% 
  hoist(m7_ob, result = list("mediation","currentsmoke", "result")) %>% 
  filter(result!="NULL") %>% 
  unnest_longer(result) %>% 
  hoist(result,p = "p")%>% 
  select(1:4,6) %>% 
  filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  mutate(id = str_c(treatment,"_",gene_set_name,"_",result_id)) %>% 
  right_join(fig1panelB %>% select(id))

currentsmoke %>% kableExtra::kable() %>% kableExtra::kable_styling()

#' ### no insurance because of no money

#+ echo=F, eval=T, warning=FALSE, message=FALSE
insurance_lack = example0_no1k %>%
  hoist(out, "result") %>%
  hoist(result, "m7_ob") %>% 
  unnest(matches("^m7")) %>% 
  hoist(m7_ob, result = list("mediation","insurance_lack", "result")) %>% 
  filter(result!="NULL") %>% 
  unnest_longer(result) %>% 
  hoist(result,p = "p")%>% 
  select(1:4,6) %>% 
  filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  mutate(id = str_c(treatment,"_",gene_set_name,"_",result_id)) %>% 
  right_join(fig1panelB %>% select(id))

insurance_lack %>% kableExtra::kable() %>% kableExtra::kable_styling()


#'  `rmarkdown::render("/home/xu/ses-1/user_wx/m7.R")`
