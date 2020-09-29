library(tidyverse)
example0 =  readRDS("/home/share/scratch/example0_without_1KI_de_novo.rds")

 panelA_data = example0 %>%
   hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
   filter(controls=="all") %>% 
   dplyr::select(treatment, gene_set_name, p) %>% 
   filter(p<0.05) %>% 
   dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5"))
 panelA_data = panelA_data %>%
   mutate(match = str_detect(gene_set_name,treatment)) %>% 
   filter(match ==TRUE)
 panelA_data = panelA_data %>% 
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
   )))
 
 threshold = 0.05/7
 panelA_data %>% dplyr::select(1:3) %>% kableExtra::kable() %>% kableExtra::kable_styling()
 
 exB <- example0 %>%
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
   dplyr::filter(p < threshold) %>% 
   mutate(match = str_detect(gene_set_name,treatment)) %>% 
   filter(match ==TRUE) %>% 
   mutate(denovo= case_when(gene_set_name=="income_hh_ff5_de_novo"  ~ "Income PC",
                            gene_set_name=="ses_sss_composite_de_novo"  ~ "SES Composite PC",
                            gene_set_name=="sss_5_de_novo"~ "Subjective Social Status PC")) %>%
   group_by(denovo) %>%
   mutate(cons_pc=1:n()) %>% 
   ungroup() %>% 
   mutate(denovo_pc= str_c(denovo, cons_pc))
 
 panelB_data = exB %>% 
   mutate(id = str_c(gene_set_name,"_",p_id),
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
   ))) 
 panelB_data %>% 
   dplyr::select(1:3, 9) %>%
   mutate(adj.p = p * 7 * 5) %>% #bonferroni correction
   kableExtra::kable() %>% kableExtra::kable_styling()

 
 panelA = ggplot(panelA_data, aes(treatment, gene_set_name, size = pval2)) +
   geom_point(stroke = 1.5, shape = 21, alpha = 0.4, fill = "darkblue", color ="darkblue") +
   theme_bw() +
   labs(
     # title = "Figure 1. Associations between Indicators of Socioeconomic Status 
     #         and mRNA-Based Disease Signatures, Add Health 
     #         (p-values reported, FDR-corrected for whole genome)",
     y = "mRNA Signatures",
     x = "SES Indicators") +
   theme(axis.text.x = element_text(angle = 30, hjust = 1),
         # plot.margin=unit(c(1, 1, 0.1, 1), "cm"),
         text = element_text(size=10, face = "bold")) +
   scale_size_continuous(name = "P-value",
                         range = c(0, 16),
                         limits = c(0.0000001, 100000), breaks = c(10000, 15000, 25000, 100000),
                         labels = c("p<0.05", "p<0.01", "p<0.001", "p<0.0001"))+
   scale_alpha(guide = 'none') +
   guides(shape = guide_legend(override.aes = list(size = 10)),
          fill = guide_legend(override.aes = list(size = 8)))
 
 # color order from bottom to top
 axiscolor = c("grey30", "grey30", "grey30", "darkseagreen1","darkseagreen4" ,"grey30","red", "darkred")
 panelB = ggplot(panelB_data, aes(treatment, denovo_pc, size = pval2)) +
   geom_point(stroke = 1.5, shape = 21, alpha = 0.4, fill = "darkblue", color ="darkblue") +
   theme_bw() +
   labs(
     y = "mRNA Signatures",
     x = "SES Indicators") +
   theme(axis.text.x = element_text(angle = 30, hjust = 1), 
         text = element_text(size=10, face = "bold")) +
   scale_size_continuous(name = "P-value",
                         range = c(0, 16),
                         limits = c(0.0000001, 100000), breaks = c(15000, 25000),
                         labels = c("p<0.005", "p<0.001"))+
   scale_alpha(guide = 'none')+
   guides(shape = guide_legend(override.aes = list(size = 10)),
          fill = guide_legend(override.aes = list(size = 8)))
 ggpubr::ggarrange(panelA, panelB, ncol = 2, labels = c("A", "B"), common.legend = TRUE, legend="right") 
 
 
 
 
