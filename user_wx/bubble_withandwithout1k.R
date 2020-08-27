
#' #### Panle A
#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(tidyverse)

example0_with1k = readRDS("/home/share/scratch/xu/example0_new2signature_withinflam.rds")
example0_without1k = readRDS("/home/share/scratch/xu/example0_new2signature_noinflam.rds")


ex0_with1k <- example0_with1k %>%
  hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
  filter(controls=="all") %>% 
  dplyr::select(treatment, gene_set_name, p) %>% 
  dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
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
         inflam = "with inflamation",
         inflam = inflam %>% as.factor())

ex0_without1k <- example0_without1k %>%
  hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
  filter(controls=="all") %>% 
  dplyr::select(treatment, gene_set_name, p) %>% 
  dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
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
         inflam = "without inflamation",
         inflam = inflam %>% as.factor())

ex_0 = bind_rows(ex0_with1k, ex0_without1k)
# need to specify shape, because some shape does not have borders. 

  # ex_0  %>% 
  # ggplot(aes(x = treatment, y = gene_set_name, size = pval2)) + 
  # geom_point(data = filter(ex_0, inflam == "with inflamation"),shape = 21, alpha = 0.4, colour = "darkred") +
  #   
  # geom_point(data = filter(dat, inflam == "without inflamation"), shape = 21, alpha = 0.4, colour = "darkblue") +
  #   theme_bw() +
  #   labs(
  #     # title = "Figure 1. Associations between Indicators of Socioeconomic Status 
  #     #         and mRNA-Based Disease Signatures, Add Health 
  #     #         (p-values reported, FDR-corrected for whole genome)",
  #     y = "mRNA Signatures",
  #     x = "SES Indicators") +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   scale_size_continuous(name = "P-value",
  #                         range = c(0, 14),
  #                         limits = c(0.0000001, 100000), breaks = c(0.0000001, 10000, 15000, 25000, 100000),
  #                         labels = c("n.s.", "p<0.05", "p<0.01", "p<0.001", "p<0.0001"))+
  #   scale_alpha(guide = 'none')
  # 
# reorder y axis
ex_0$gene_set_name <- factor(ex_0$gene_set_name, levels = c("Inflam1k",
                                                            "Alzheimers","Aortic Aneurysm", "Asthma", 
                                                            "CKD", "COPD", "CVD",
                                                            "Depression", "Diabetes",
                                                            "Hypertension",
                                                            "Rheumatoid Arthritis"
                                                             )) %>% fct_rev
ggplot(ex_0, aes(treatment, gene_set_name, size = pval2, fill = inflam, colour = inflam)) +
  geom_point(stroke = 1.5, shape = 21, alpha = 0.4) +
  scale_fill_manual(values = c("red", "cornflowerblue")) +
  scale_color_manual(values = c( "darkred","darkblue")) +
  # geom_jitter(height = 0.00000025) +
  # gghighlight(class == "inflam") +
  theme_bw() +
  labs(
    # title = "Figure 1. Associations between Indicators of Socioeconomic Status 
    #         and mRNA-Based Disease Signatures, Add Health 
    #         (p-values reported, FDR-corrected for whole genome)",
    y = "mRNA Signatures",
    x = "SES Indicators") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_size_continuous(name = "P-value",
                        range = c(0, 14),
                        limits = c(0.0000001, 100000), breaks = c(0.0000001, 10000, 15000, 25000, 100000),
                        labels = c("n.s.", "p<0.05", "p<0.01", "p<0.001", "p<0.0001"))+
  scale_alpha(guide = 'none')

#' #### Panle A result to see the difference of removing and keeping inflam1k singatures, only a couple of difference
#+ echo=F, eval=F, warning=FALSE, message=FALSE

with1k = ex0_with1k %>% mutate(p_withinflam = p) %>% select(1,2,7)
without1k = ex0_without1k %>% mutate(p_withoutinflam = p) %>% select(1,2,7)
  with1k %>% left_join(without1k, by.x = treatment, by.y = gene_set_name) %>% 
    kableExtra::kable() %>% kableExtra::kable_styling()

#' #### Panle B, pca dimention is totally different things before and after removing inflam1k singatures, so I remove the dimension information
#' and only keep the signature names
  
#+ echo=F, eval=T, warning=FALSE, message=FALSE
  
  threshold_with1k = 0.05/11/4
  exB <- example0_with1k %>%
    hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
    unnest_longer(p) %>% 
    dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
    mutate(pval2=case_when(p<0.0001 ~ 100000,
                           p<0.001 ~ 25000,
                           p<0.01 ~ 15000,
                           p<0.05 ~ 10000,
                           p>0.05 ~ 0.0000001)
    ) %>% 
    dplyr::select(treatment, gene_set_name, p, p_id, pval2) %>% 
    dplyr::filter(p < threshold_with1k)
  
  exB_data_with1k = exB %>% 
    mutate(gene_set_name = str_c(gene_set_name,"_",p_id),
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
    mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA","")) %>% 
    mutate(inflam = "with inflamation") %>% 
    mutate(gene_set_name = gene_set_name %>% str_remove("_d.*$") %>% str_trim,
           gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Cvd", "CVD"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Ckd", "CKD"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Copd", "COPD"))
  
  threshold_without1k = threshold_with1k
  exB <- example0_without1k %>%
    hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
    unnest_longer(p) %>% 
    dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
    mutate(pval2=case_when(p<0.0001 ~ 100000,
                           p<0.001 ~ 25000,
                           p<0.01 ~ 15000,
                           p<0.05 ~ 10000,
                           p>0.05 ~ 0.0000001)
    ) %>% 
    dplyr::select(treatment, gene_set_name, p, p_id, pval2) %>% 
    dplyr::filter(p < threshold_without1k) 
  
  exB_data_without1k = exB %>% 
    mutate(gene_set_name = str_c(gene_set_name,"_",p_id),
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
    mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA","")) %>% 
    mutate(inflam = "without inflamation") %>% 
    mutate(gene_set_name = gene_set_name %>% str_remove("_d.*$") %>% str_trim,
           gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Cvd", "CVD"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Ckd", "CKD"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Copd", "COPD"))
  
  exB_data = full_join(exB_data_with1k %>% select(1,2,3,6), exB_data_without1k %>% select(1,2,3,6), by.x = treatment, by.y = "gene_set_name") %>% 
    mutate(pval2 = case_when(p<0.0001 ~ 100000,
                           p<0.001 ~ 25000,
                           p<0.01 ~ 15000,
                           p<0.05 ~ 10000,
                           p>0.05 ~ 0.0000001),
           gene_set_name = factor(gene_set_name, levels = c("Inflam1k",
                                                                       "Alzheimers","Aortic Aneurysm", "Asthma", 
                                                                       "CKD", "COPD", "CVD",
                                                                       "Depression", "Diabetes",
                                                                       "Hypertension",
                                                                       "Rheumatoid Arthritis"
           )) %>% fct_rev
    ) 
  
  # ggplot(exB_data, aes(treatment, gene_set_name, size = pval2, alpha = 0.4)) +
  #   geom_point(fill = "red", color="navy") +
  ggplot(exB_data, aes(treatment, gene_set_name, size = pval2, fill = inflam, colour = inflam)) +
    geom_point(stroke = 1.5, shape = 21, alpha = 0.4) +
    scale_fill_manual(values = c("red", "cornflowerblue")) +
    scale_color_manual(values = c( "darkred","darkblue")) +  
    theme_bw() +
    labs(
      # title = "Figure 1. Associations between Indicators of Socioeconomic Status 
      #         and mRNA-Based Disease Signatures, Add Health 
      #         (p-values reported, FDR-corrected for whole genome)",
      y = "mRNA Signatures",
      x = "SES Indicators") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_size_continuous(name = "P-value",
                          range = c(0, 14),
                          limits = c(0.0000001, 100000), breaks = c(15000, 25000),
                          labels = c("p<0.005", "p<0.001"))+
    scale_alpha(guide = 'none')
#+ echo=F, eval=F, warning=FALSE, message=FALSE 

  example0_noerror = remove_errors(example0_with1k)
  example0_m7_ob = example0_noerror %>% get_sig_PCs_and_sig_enrichment_on_those_PCs("m7_ob", threshold = threshold_with1k)
  m7_ob = example0_m7_ob %>%
    filter(map_lgl(well_loaded_genes_on_significant_PCs, ~ length(.x)!=0))%>% 
    filter(treatment != "ses_composite_ff5")
  
  all_pathway = m7_ob %>%
    unnest_longer(enrichment_of_well_loaded_genes_on_significant_PCs) %>%
    hoist(enrichment_of_well_loaded_genes_on_significant_PCs, pathway = list("out","enriched_physiology","reactome")) 
  
  temp = all_pathway$pathway %>%
    map(function(x) as_tibble(x) %>% dplyr::select(Description, p.adjust) %>% dplyr::top_n(-3)) %>%
    set_names(all_pathway$gene_set_name) %>% 
    map_df(~as.data.frame(.x), .id="Signature") %>%
    as_tibble() %>% 
    unique()
  
  rownames(temp) <- NULL
  
  figure = ggarrange(
    ggarrange(panelA, panelB, ncol = 2, labels = c("A", "B"), common.legend = TRUE, legend="right"), 
    panelC,
    nrow = 2,
    heights = c(2, 0.7),
    labels = c("","C")       # Label of the line plot
  )
  
  title <- expression(atop(bold("Figure 1:"),
                           scriptstyle("panel A:Associations between Indicators of Socioeconomic Status and mRNA-Based Disease Signatures, Add Health (p-values reported, FDR-corrected for whole genome); panel B: Significant PCs; panel C: pathways")))
  
  # title <- expression(atop(bold("Figure 1:"),
  # scriptstyle("panel A:Associations between Indicators of Socioeconomic Status and mRNA-Based Disease Signatures, Add Health (p-values reported, FDR-corrected for whole genome) \n panel B: Significant PCs \n panel C: pathways")))
  
  annotate_figure(figure,bottom = text_grob(title)
                  # bottom = text_grob("Figure 1: \n panel A: Associations between Indicators of Socioeconomic Status and mRNA-Based Disease Signatures, Add Health (p-values reported, FDR-corrected for whole genome \n panel B: Significant PCs \n panel C: pathways")
                  #           fig.lab = "Figure 1:
                  #           panel A - Associations between Indicators of Socioeconomic Status
                  #      and mRNA-Based Disease Signatures, Add Health
                  # (p-values reported, FDR-corrected for whole genome)",
                  #           fig.lab.face = "bold",
                  #           fig.lab.pos = "bottom.right"
  )
  
  
#'`rmarkdown::render("/home/xu/ses-1/user_wx/bubble_withandwithout1k.R")`