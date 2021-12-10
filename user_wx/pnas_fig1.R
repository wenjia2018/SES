# extract ses omnibus analysis results
# load results dataset
example0_with1k  <- readRDS("~/ses-1/user_wx/m12_withinflam_ses_3setscontrols.rds")
example0_without1k  <- readRDS("~/ses-1/user_wx/m12_withoutinflam_ses_3setscontrols.rds")
# specify which control set, "basic", "basic_less", "basic_less_ancestry"
control_set = "basic"
table1 =
  c(
    "CVD_mRNA",
    "diabetes_mRNA",
    "inflam1k_mRNA",
    # "breast_cancer_mRNA",
    # "Lupus_mRNA", "Colorectal_mRNA",
    "Rheumatoid_Arthritis_mRNA", "Alzheimers_mRNA",
    "Aortic_Aneurysm_mRNA", "COPD_mRNA",
    "Asthma_mRNA","Hypertension_mRNA",
    "Depression_mRNA",
    "CKD_mRNA"
    # "kidney_transplant_tolerance_mRNA"
  )

ex0_without1k <- example0_without1k %>%
  hoist(out, p = list("result", "m12_fdr", 1, "p")) %>% 
  filter(controls ==control_set) %>% 
  dplyr::select(treatment, gene_set_name, p) %>% 
  filter(p<0.05) %>% 
  mutate(pval=case_when(p<0.0001 ~ 0.0001,
                        p<0.001 ~ 0.001,
                        p<0.01 ~ 0.01,
                        p<0.05 ~ 0.05,
                        p>0.05 ~ 100),
         pval2=case_when(p<0.0001 ~ 100000,
                         p<0.001 ~ 25000,
                         p<0.01 ~ 15000,
                         p<0.05 ~ 10000,
                         p>0.05 ~ 0.0000001),
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
                              treatment =="ses_composite_ff5"  ~ "SES Composite 3"))  %>%
  mutate(treatment = factor(treatment, levels = c("Parental SES Composite", "Parental Education","Parental Income",
                                                  "Parental SEI", "Mother's Occupation", "Father's Occupation",
                                                  "SES Composite 3", "SES Composite", "Education", "Income",
                                                  "Occupation", "Subjective Social Status","Occupation Work Collar"
  ))) %>%
  mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
         gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
         gene_set_name = gene_set_name %>% replace(gene_set_name == "Cvd", "CVD"),
         gene_set_name = gene_set_name %>% replace(gene_set_name == "Copd", "COPD"),
         gene_set_name = gene_set_name %>% replace(gene_set_name == "Ckd", "CKD"),
         gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI"),
         "1KI Genes" = "Without 1KI Genes"%>% as.factor())


ex0_with1k <- example0_with1k %>%
  hoist(out, p = list("result", "m12_fdr", 1, "p")) %>% 
  filter(controls == control_set) %>% 
  dplyr::select(treatment, gene_set_name, p) %>% 
  filter(p<0.05) %>% 
  mutate(pval=case_when(p<0.0001 ~ 0.0001,
                        p<0.001 ~ 0.001,
                        p<0.01 ~ 0.01,
                        p<0.05 ~ 0.05,
                        p>0.05 ~ 100),
         pval2=case_when(p<0.0001 ~ 100000,
                         p<0.001 ~ 25000,
                         p<0.01 ~ 15000,
                         p<0.05 ~ 10000,
                         p>0.05 ~ 0.0000001),
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
                              treatment =="ses_composite_ff5"  ~ "SES Composite 3"))  %>%
  mutate(treatment = factor(treatment, levels = c("Parental SES Composite", "Parental Education","Parental Income",
                                                  "Parental SEI", "Mother's Occupation", "Father's Occupation",
                                                  "SES Composite 3", "SES Composite", "Education", "Income",
                                                  "Occupation", "Subjective Social Status","Occupation Work Collar"
  ))) %>%
  mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
         gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
         gene_set_name = gene_set_name %>% replace(gene_set_name == "Cvd", "CVD"),
         gene_set_name = gene_set_name %>% replace(gene_set_name == "Copd", "COPD"),
         gene_set_name = gene_set_name %>% replace(gene_set_name == "Ckd", "CKD"),
         gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI"),
         "1KI Genes" = "With 1KI Genes"%>% as.factor())


ex_0 = bind_rows(ex0_with1k, ex0_without1k)
# reorder y axis
ex_0$gene_set_name <- factor(ex_0$gene_set_name, levels = c("1KI",
                                                            "Alzheimers","Aortic Aneurysm", "Asthma", 
                                                            "CKD", "COPD", "CVD",
                                                            "Depression", "Diabetes",
                                                            "Hypertension",
                                                            "Rheumatoid Arthritis"
)) %>% fct_rev

ex_0 %>% saveRDS("./xxfig1.rds")
