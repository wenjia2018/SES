
#' ---
#' title: mediational Evalue for PCA
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE, include=FALSE, comment=NA

#' ### Evalue
#' *  The E-value is defined as the minimum strength of association, on the risk ratio scale,
#'  that an unmeasured confounder would need to have with both the treatment and the outcome to fully explain away a specific
#'  treatment–outcome association, conditional on the measured covariates.
#' * Input required to calculate Evalue: OLS estimate and standard deviation, and standard deviation of the model (outcome)
#'
#' ### approximation or speical treatment made to calculate evalue for causal mediated effect
#' * estimate from mediation package is based on bootstrap. use bootstrap estimate and sd as OLS estimate and sd.  

#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
library(tidyverse)
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(Biobase)
library(ggformula)
library(ggpubr)
library(here)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(tidyverse)
library(stringi)
library(kableExtra)
walk(dir(path = here("R"),full.names = TRUE), source)
example0_with1k <- readRDS("~/ses-1/user_wx/example_tmm_m7_withinflame.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/example_tmm_m7_noinflame0209.rds")

threshold = 0.05
mediators = 
  c(
    "stress_perceived_lm",
    "bills_binary",
    "currentsmoke_binary",
    "w5bmi_lm",
    "insurance_lack_binary"
  )
threshold = 0.05
extract_pp = function(m, fig1panelB, mediator){
  m %>%
    hoist(out, "result") %>%
    hoist(result, "m7_ob") %>% 
    unnest(matches("^m7")) %>% 
    hoist(m7_ob, result = list("mediation", mediator, "result")) %>% 
    filter(result!="NULL") %>% 
    unnest_longer(result) %>% 
    hoist(result, p = "p")%>% 
    # hoist(result, detail = "other") %>% 
    # unnest_wider("detail") %>% 
    dplyr::select(-out, -m7_ob) %>%
    # filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
    mutate(id = str_c(treatment,"_",gene_set_name,"_",result_id)) %>% 
    right_join(fig1panelB %>% dplyr::select(id)) %>%
    mutate(p = p %>% str_remove("<") %>% as.numeric() %>% format(digits = 3, scientific =T),
           mediator = mediator)
  # %>% 
  # adjP = ((p %>% str_remove("<") %>% as.numeric())) %>% p.adjust("fdr")) %>% 
  # filter(adjP<0.05) %>% 
  # mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric
}
med_out_all = function(example_data) {
  
  # we use fig1panleB to choose the mediation we wanted to do, that is only for significant total effect
  fig1panelB <- example_data  %>%
    hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
    unnest_longer(p) %>% 
    # dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
    mutate(p = p.adjust(p, method = "fdr")) %>% 
    dplyr::filter(p < threshold) %>% 
    dplyr::select(treatment, gene_set_name, p_id) %>% 
    mutate(id = str_c(treatment,"_",gene_set_name,"_",p_id))
  
  
  out = mediators %>% set_names() %>%  map(~ extract_pp(example_data, fig1panelB, .x)) %>% bind_rows() %>% 
    mutate(adjP = ((p %>% str_remove("<") %>% as.numeric())) %>% p.adjust("fdr")) %>% 
    filter(adjP<0.05) %>%
    mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric
  # out %>% select(-id, -controls) %>% kableExtra::kable() %>% kableExtra::kable_styling()
}


exB_data_with = 
  example0_with1k %>%
  med_out_all() 

with = 
  exB_data_with %>% 
  # mutate(pcmin = result_id %>% str_remove("d") %>% as.numeric()) %>% 
  # group_by(treatment, gene_set_name, controls, mediator) %>% 
  # slice(which.min(adjP)) %>% 
  # ungroup %>% 
  mutate(pcmin = result_id %>% str_remove("d") %>% as.numeric()) %>%
  # select(-c(6,7,8,9)) %>% 
  group_by(treatment, gene_set_name, controls, mediator) %>%
  mutate(p_no = n()) %>% 
  slice(which.min(pcmin)) %>% 
  ungroup
# m_evalue = function(ols_est, sd) {
#   est = toRR(OLS(ols_est, sd))
#   return(est + sqrt(est*(est-1)))
# }
temp2 = with %>%
  mutate(
    mediator = case_when(mediator =="w5bmi_lm" ~ "BMI",
                         mediator =="bills_binary" ~ "Financial Stress",
                         mediator =="currentsmoke_binary" ~ "Current Smoking Status",
                         mediator =="stress_perceived_lm" ~ "Perceived Stress",
                         mediator =="insurance_lack_binary" ~ "Lack of Insurance"),
    treatment = case_when(treatment == "edu_p" ~ "Parental Education",
                         treatment =="income_pp1_log" ~  "Parental Income" ,
                         treatment =="SEI_max_p_w12" ~ "Parental SEI",
                         treatment =="ses_composite_pp1" ~ "Parental SES Composite",
                         treatment =="work_collar_rm_f12" ~ "Mother's Occupation",
                         treatment =="work_collar_rf_f12" ~ "Father's Occupation" ,
                         treatment =="work_collar_ff5" ~ "Occupation Work Collar",
                         treatment =="edu_max" ~ "Education" ,
                         treatment =="income_hh_ff5" ~ "Income"     ,
                         treatment =="SEI_ff5" ~ "Occupation"      ,
                         treatment =="ses_sss_composite" ~ "SES Composite"  ,
                         treatment =="sss_5" ~ "Subjective Social Status",
                         treatment =="ses_composite_ff5"  ~ "SES Composite 3")) %>% 
  unnest_wider(result) %>% 
  mutate(`E-values` = evalue %>% map_dbl(~ .x[2,1]))
# %>%
  # unnest_wider(sims) %>% 
  # unnest_wider(other) %>% 
  # mutate(med_ACME = med_ACME%>% as.numeric()) %>% 
  # mutate(med_evalue = m_evalue(med_ACME, y_sd))

# temp2 %>%
#   select(treatment, gene_set_name, mediator,`E-values`) %>% 
#   pivot_wider(names_from = treatment, values_from = `E-values`) %>% 
#   mutate(across(3:7, ~ .x %>% as.numeric() %>% format(digits = 4))) %>%
#   mutate(across(everything(), ~ ifelse(.x=="   NA", "", .x))) %>% 
#   kableExtra::kbl() %>%
#   kableExtra::kable_classic() %>%
#   add_header_above(c(" " = 1, "Group 1" = 2, "Group 2" = 2, "Group 3" = 2))



temp2 %>%
  select(treatment, gene_set_name, mediator,`E-values`) %>% 
  pivot_wider(names_from = treatment, values_from = `E-values`) %>% 
  mutate(across(3:7, ~ .x %>% as.numeric() %>% format(digits = 4))) %>%
  mutate(across(everything(), ~ ifelse(.x=="   NA", "", .x))) %>% 
  arrange(mediator) %>% 
  select(-mediator) %>% 
  kbl() %>% 
  kable_paper("striped", full_width = F) %>%
  pack_rows("BMI", 1, 10,label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Current Smoking Status", 11, 19,label_row_css = "background-color: #666; color: #fff;") %>% 
  pack_rows("Financial Stress", 20, 27,label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Perceived Stress", 28, 31,label_row_css = "background-color: #666; color: #fff;") %>% 
  add_indent(c(1:31))
  # pivot_wider(names_from = "mediator",values_from = c(3:7)) %>% 
  # kableExtra::kbl() %>%
  # kableExtra::kable_classic() %>%
  # add_header_above(c(" " = 1, "Education" = 4, "Income" = 4,"Occupation" = 4, "SES Composite" = 4, "Subjective Social Status" = 4)) 
temp2 %>%
  select(treatment, gene_set_name, mediator,`E-values`) %>% 
  pivot_wider(names_from = treatment, values_from = `E-values`) %>% 
  mutate(across(3:7, ~ .x %>% as.numeric() %>% format(digits = 4))) %>%
  mutate(across(everything(), ~ ifelse(.x=="   NA", "", .x))) %>% 
  arrange(mediator) %>% 
  # select(-mediator) %>% 
  kbl() %>% 
  kable_paper("striped", full_width = F) %>%
  pack_rows("BMI", 1, 9,label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Current Smoking Status", 10, 17,label_row_css = "background-color: #666; color: #fff;") %>% 
  pack_rows("Financial Stress", 18, 24,label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Perceived Stress", 25, 27,label_row_css = "background-color: #666; color: #fff;") %>% 
  add_indent(c(1:27))
