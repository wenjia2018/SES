
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


#' We do PCA for two new disease signature, which are constructed by DE ses, using oblique rotation and then do mediational analysis
#'  for significant PCs using 5 mediators
#'  
#' p: p value for indirected effect / mediated effect;
#' result_id: which pc is significant;
#' adjP is adjusted for 13 disease signatures and 5 mediators, ie. adjP = p /65;
#' results are presented for adjP<0.05 
library(here)
library(tidyverse)
# example0_with1k = readRDS("/home/share/scratch/xu/example0_new2signature_withinflam.rds")

example0_with1k <- readRDS("/home/share/scratch/m7_desig_withinflam.rds")
#' ## keep inflamation signatures
#' 
threshold_with1k = 0.05/10/5
threshold_med = 13*5
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
  right_join(fig1panelB %>% select(id))  %>% 
  mutate(p = p %>% str_remove("<") %>% as.numeric() %>% format(digits = 3, scientific =T),
         adjP = ((p %>% str_remove("<") %>% as.numeric()) * threshold_med),
         adjP = ifelse(adjP>1, 1, adjP)) %>% 
  filter(adjP<0.05) %>% 
  mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric
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
  right_join(fig1panelB %>% select(id)) %>% 
  mutate(p = p %>% str_remove("<") %>% as.numeric() %>% format(digits = 3, scientific =T),
         adjP = ((p %>% str_remove("<") %>% as.numeric()) * threshold_med),
         adjP = ifelse(adjP>1, 1, adjP)) %>% 
  filter(adjP<0.05) %>% 
  mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric

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
  right_join(fig1panelB %>% select(id))  %>% 
  mutate(p = p %>% str_remove("<") %>% as.numeric() %>% format(digits = 3, scientific =T),
         adjP = ((p %>% str_remove("<") %>% as.numeric()) * threshold_med),
         adjP = ifelse(adjP>1, 1, adjP)) %>% 
  filter(adjP<0.05) %>% 
  mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric

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
  right_join(fig1panelB %>% select(id))  %>% 
  mutate(p = p %>% str_remove("<") %>% as.numeric() %>% format(digits = 3, scientific =T),
         adjP = ((p %>% str_remove("<") %>% as.numeric()) * threshold_med),
         adjP = ifelse(adjP>1, 1, adjP)) %>% 
  filter(adjP<0.05) %>% 
  mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric

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
  right_join(fig1panelB %>% select(id))  %>% 
  mutate(p = p %>% str_remove("<") %>% as.numeric() %>% format(digits = 3, scientific =T),
         adjP = ((p %>% str_remove("<") %>% as.numeric()) * threshold_med),
         adjP = ifelse(adjP>1, 1, adjP)) %>% 
  filter(adjP<0.05) %>% 
  mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric

insurance_lack %>% kableExtra::kable() %>% kableExtra::kable_styling()

#' ## remove inflamation signatures


example0_no1k = readRDS("/home/share/scratch/xu/example0_new2signature_noinflam.rds")
example0_no1k = readRDS("/home/share/scratch/m7_desig_noinflam.rds")

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
  right_join(fig1panelB %>% select(id))  %>% 
  mutate(p = p %>% str_remove("<") %>% as.numeric() %>% format(digits = 3, scientific =T),
         adjP = ((p %>% str_remove("<") %>% as.numeric()) * threshold_med),
         adjP = ifelse(adjP>1, 1, adjP)) %>% 
  filter(adjP<0.05) %>% 
  mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric

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
  right_join(fig1panelB %>% select(id))  %>% 
  mutate(p = p %>% str_remove("<") %>% as.numeric() %>% format(digits = 3, scientific =T),
         adjP = ((p %>% str_remove("<") %>% as.numeric()) * threshold_med),
         adjP = ifelse(adjP>1, 1, adjP)) %>% 
  filter(adjP<0.05) %>% 
  mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric

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
  right_join(fig1panelB %>% select(id))  %>% 
  mutate(p = p %>% str_remove("<") %>% as.numeric() %>% format(digits = 3, scientific =T),
         adjP = ((p %>% str_remove("<") %>% as.numeric()) * threshold_med),
         adjP = ifelse(adjP>1, 1, adjP)) %>% 
  filter(adjP<0.05) %>% 
  mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric
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
  right_join(fig1panelB %>% select(id)) %>% 
  mutate(p = p %>% str_remove("<") %>% as.numeric() %>% format(digits = 3, scientific =T),
         adjP = ((p %>% str_remove("<") %>% as.numeric()) * threshold_med),
         adjP = ifelse(adjP>1, 1, adjP)) %>% 
  filter(adjP<0.05) %>% 
  mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric

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
  right_join(fig1panelB %>% select(id)) %>% 
  mutate(p = p %>% str_remove("<") %>% as.numeric() %>% format(digits = 3, scientific =T),
         adjP = ((p %>% str_remove("<") %>% as.numeric()) * threshold_med),
         adjP = ifelse(adjP>1, 1, adjP)) %>% 
  filter(adjP<0.05) %>% 
  mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric

insurance_lack %>% kableExtra::kable() %>% kableExtra::kable_styling()


#'  `rmarkdown::render("/home/xu/ses-1/user_wx/m7.R")`
