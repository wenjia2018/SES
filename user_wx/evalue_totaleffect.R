# evalue total effect
with <- readRDS("~/ses-1/user_wx/Evalue_cluster_with1k.rds")
without <- readRDS("~/ses-1/user_wx/Evalue_cluster_without1k_v2.rds")

extract_total_eval = function(data){
  data %>% 
    hoist(out, eval = list("result")) %>% 
    unnest(eval) %>% 
    mutate(Evalue = eval %>% map_dbl(~ .x[2,1])) %>% 
    select(treatment, gene_set_name, Evalue) %>% 
    group_by(treatment, gene_set_name) %>% 
    mutate(`E-values` = mean(Evalue)) %>% 
    ungroup %>% 
    select(treatment, gene_set_name, `E-values`) %>% 
    unique() %>% 
    mutate(
      treatment= case_when(treatment == "edu_p" ~ "Parental Education",
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
}
temp_with = with %>% extract_total_eval()
temp_with %>%   
  pivot_wider(names_from = treatment, values_from = `E-values`) %>% 
  mutate(across(2:6, ~ .x %>% as.numeric() %>% format(digits = 4))) %>%
  mutate(across(everything(), ~ ifelse(.x=="   NA", "", .x))) %>% 
  # mutate(across(2:6, ~ .x %>% as.numeric() %>% format(digits = 4))) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()

temp_without = without %>% extract_total_eval()
temp_without %>% 
  pivot_wider(names_from = treatment, values_from = `E-values`) %>% 
  mutate(across(2:6, ~ .x %>% as.numeric() %>% format(digits = 4))) %>%
  mutate(across(3:6, ~ ifelse(.x=="   NA", "", .x))) %>% 
  mutate(across(everything(), ~ ifelse(.x=="  NA", "", .x))) %>% 
  # mutate(across(2:6, ~ .x %>% as.numeric() %>% format(digits = 4))) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()
