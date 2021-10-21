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
example0_with1k <- readRDS("~/ses-1/user_wx/mediate_cluster_with1k_v2.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/mediate_cluster_without1k.rds")

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

extract_pp = function(m, mediator){
  m %>%
    hoist(out, result = list("result", mediator, "result")) %>% 
    filter(result!="NULL") %>% 
    unnest_longer(result) %>% 
    hoist(result, p = "p")%>% 
    hoist(result, evalue = "evalue") %>% 
    # hoist(result, detail = "other") %>% 
    # unnest_wider("detail") %>% 
    dplyr::select(-out, -result) %>% 
    mutate(p = p %>% str_remove("<") %>% as.numeric() %>% format(digits = 3, scientific =T),
           mediator = mediator)
  # %>% 
  # adjP = ((p %>% str_remove("<") %>% as.numeric())) %>% p.adjust("fdr")) %>% 
  # filter(adjP<0.05) %>% 
  # mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric
}
med_out_all = function(example_data) {
  out = mediators %>% set_names() %>%  map(~ extract_pp(example_data, .x)) %>% bind_rows() %>% 
    mutate(adjP = ((p %>% str_remove("<") %>% as.numeric())) %>% p.adjust("fdr")) %>% 
    filter(adjP<0.05) %>%
    mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric
  # out %>% select(-id, -controls) %>% kableExtra::kable() %>% kableExtra::kable_styling()
}


exB_data_with = 
  example0_with1k %>%
  med_out_all() %>% 
  mutate(eval = evalue %>% map_dbl(~ .x[2,1]),
         ) %>% 
  select(treatment, gene_set_name, result_id, mediator, eval)

with = 
  exB_data_with %>% 
  group_by(treatment, gene_set_name, mediator) %>% 
  mutate(`E-values` = mean(eval)) %>% 
  ungroup %>% 
  select(treatment, gene_set_name, mediator, `E-values`) %>% 
  unique()

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
                          treatment =="ses_composite_ff5"  ~ "SES Composite 3"))


temp2 %>%
  select(treatment, gene_set_name, mediator,`E-values`) %>% 
  pivot_wider(names_from = treatment, values_from = `E-values`) %>% 
  mutate(across(3:7, ~ .x %>% as.numeric() %>% format(digits = 4))) %>%
  mutate(across(everything(), ~ ifelse(.x=="   NA", "", .x))) %>% 
  arrange(mediator) %>% 
  select(-mediator) %>% 
  kbl() %>% 
  kable_paper("striped", full_width = F) %>%
  pack_rows("BMI", 1, 9,label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Current Smoking Status", 10, 17,label_row_css = "background-color: #666; color: #fff;") %>% 
  pack_rows("Financial Stress", 18, 24,label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Perceived Stress", 25, 28,label_row_css = "background-color: #666; color: #fff;") %>% 
  add_indent(c(1:28))
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
  select(-mediator) %>%
  select(gene_set_name, Education, Income, Occupation, `SES Composite`, `Subjective Social Status`) %>% 
  kbl() %>% 
  kable_paper("striped", full_width = F) %>%
  pack_rows("BMI", 1, 7,label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Current Smoking Status", 8, 12,label_row_css = "background-color: #666; color: #fff;") %>% 
  pack_rows("Financial Stress", 13, 16,label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Perceived Stress", 17, 18,label_row_css = "background-color: #666; color: #fff;") %>% 
  add_indent(c(1:18))
