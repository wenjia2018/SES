
#' ###  with1k
#+ echo=F, eval=T, warning=FALSE, message=FALSE
example <- readRDS("~/ses-1/user_wx/mediate_cluster_multiple_with1k.rds")

 with = 
   example %>% 
  hoist(out, proporp = list ("result")) %>% 
  unnest_longer(proporp) %>% 
  filter(!is.na(proporp_id)) %>%
  dplyr::mutate(multimed = map_dbl(.$proporp, ~ .x[2,2]) %>% abs()) %>% 
  group_by(treatment, gene_set_name) %>% 
  dplyr::mutate(proportion = mean(multimed)) %>% 
  ungroup %>% 
  dplyr::mutate(names = str_c(treatment,gene_set_name, proporp_id)) %>% 
  select(treatment, gene_set_name, proportion) %>% 
  mutate(proportion = proportion %>% scales::percent(accuracy = 0.1),
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
  unique 
 
 with %>% 
  kableExtra::kable() %>% kableExtra::kable_styling()

#' ###  without1k
#+ echo=F, eval=T, warning=FALSE, message=FALSE
example <- readRDS("~/ses-1/user_wx/mediate_cluster_multiple_without1k.rds")

without = 
  example %>% 
  hoist(out, proporp = list ("result")) %>% 
  unnest_longer(proporp) %>% 
  filter(!is.na(proporp_id)) %>%
  dplyr::mutate(multimed = map_dbl(.$proporp, ~ .x[2,2]) %>% abs()) %>% 
  group_by(treatment, gene_set_name) %>% 
  dplyr::mutate(proportion = mean(multimed)) %>% 
  ungroup %>% 
  dplyr::mutate(names = str_c(treatment,gene_set_name, proporp_id)) %>% 
  select(treatment, gene_set_name, proportion) %>% 
  mutate(proportion = proportion %>% scales::percent(accuracy = 0.1),
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
  unique 

without %>% 
  kableExtra::kable() %>% kableExtra::kable_styling()



example <- readRDS("~/ses-1/user_wx/mediate_cluster_multiple_DeNovo.rds")

denovo = 
  example %>% 
  hoist(out, proporp = list ("result")) %>% 
  unnest_longer(proporp) %>% 
  filter(!is.na(proporp_id)) %>%
  dplyr::mutate(multimed = map_dbl(.$proporp, ~ .x[2,2]) %>% abs()) %>% 
  group_by(treatment, gene_set_name) %>% 
  dplyr::mutate(proportion = mean(multimed)) %>% 
  ungroup %>% 
  dplyr::mutate(names = str_c(treatment,gene_set_name, proporp_id)) %>% 
  select(treatment, gene_set_name, proportion) %>% 
  mutate(proportion = proportion %>% scales::percent(accuracy = 0.1),
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
  unique 
# 0.6794162
denovo %>% 
  kableExtra::kable() %>% kableExtra::kable_styling()
#' #' ### multiple mediation results
#' #+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
#' for (i in 1 : dim(ses)[1]) {
#' 
#'   cat(" ############################################################","\n",
#'       "treatment is", ses[i, ]$treatment, "gene_set_name is", ses[i, ]$gene_set_name,"\n",
#'       "############################################################","\n")
#'   cat("multiple mediators","\n","\n")
#' 
#'   print(ses[i,]$multimed)
#'   cat("\n","\n")
#' 
#' }


