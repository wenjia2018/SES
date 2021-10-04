



example_tmm_m7_denovosesDE <- readRDS("~/ses-1/user_wx/example_tmm_m7_denovosesDE.rds") %>% 
  filter(gene_set_name == "ses4_de_novo") %>% 
  mutate(gene_set_name = case_when(gene_set_name =="ses4_de_novo" ~ "SES Composite de novo without all signatures",
                                   gene_set_name =="ses4_with1k_de_novo" ~ "SES Composite de novo"))
exB <- example_tmm_m7_denovosesDE %>%
  hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
  unnest_longer(p) %>% 
  # dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  mutate(p = p.adjust(p, method = "fdr"),
         pval2=case_when(p<0.0001 ~ 100000,
                         p<0.001 ~ 25000,
                         p<0.01 ~ 15000,
                         p<0.05 ~ 10000,
                         p>0.05 ~ 0.0000001)
  ) %>% 
  dplyr::select(treatment, gene_set_name, p, p_id, pval2) %>% 
  dplyr::filter(p < threshold) 
# %>% 
  # group_by(treatment, gene_set_name) %>% 
  # # slice(which.min(p)) %>% 
  # mutate(pcmin = p_id %>% str_remove("d") %>% as.numeric()) %>%
  # group_by(treatment, gene_set_name) %>%
  # mutate(p_no = n()) %>% 
  # slice(which.min(pcmin)) %>% 
  # ungroup

exB_data = exB%>% 
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
  ))) 



ggplot(exB_data, aes(treatment, gene_set_name, size = pval2)) +
  geom_point(stroke = 1.5, shape = 21, alpha = 0.4, fill = "darkblue" , colour = "darkblue") +
  scale_fill_manual(values = c("darkblue", "goldenrod3")) +
  scale_color_manual(values = c("darkblue", "goldenrod3")) +
  # scale_fill_manual(values = c("red", "cornflowerblue")) +
  # scale_color_manual(values = c( "darkred","darkblue")) +  
  theme_bw() +
  labs(
    # title = "Figure 1. Associations between Indicators of Socioeconomic Status 
    #         and mRNA-Based Disease Signatures, Add Health 
    #         (p-values reported, FDR-corrected for whole genome)",
    y = "mRNA Signatures",
    x = "SES Indicators") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), 
        # plot.margin=unit(c(1, 1, 0.1, 1), "cm"),
        # axis.text.y = element_text(color=axiscolor), 
        text = element_text(size=10, face = "bold")) +
  scale_size_continuous(name = "P-value",
                        range = c(0, 16),
                        # limits = c(1000, 100000), breaks = c(1000, 5000, 25000,100000),
                        limits = c(0.0000001, 100000), breaks = c(10000, 15000, 25000, 100000),
                        labels = c("p<0.05", "p<0.01","p<0.001","p<0.0001"))+
  scale_alpha(guide = 'none')+
  guides(shape = guide_legend(override.aes = list(size = 10)),
         fill = guide_legend(override.aes = list(size = 8)))
mediators = 
  c(
    "stress_perceived_lm",
    "bills_binary",
    "currentsmoke_binary",
    "w5bmi_lm",
    "insurance_lack_binary"
  )
threshold = 0.05
extract_pp = function(m,fig1panelB, mediator){
  m %>%
    hoist(out, "result") %>%
    hoist(result, "m7_ob") %>% 
    unnest(matches("^m7")) %>% 
    hoist(m7_ob, result = list("mediation", mediator, "result")) %>% 
    filter(result!="NULL") %>% 
    unnest_longer(result) %>% 
    hoist(result, p = "p")%>% 
    hoist(result, detail = "other") %>% 
    unnest_wider("detail") %>% 
    dplyr::select(-out, -result, -m7_ob) %>% 
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
  example_tmm_m7_denovosesDE %>%
  med_out_all() %>%
  mutate(proportion = scales::percent(as.numeric(med_prop), accuracy = 0.1))

with = 
  exB_data_with %>% 
  mutate(pcmin = result_id %>% str_remove("d") %>% as.numeric()) %>%
  dplyr::select(-c(5:8)) 
# %>% 
  # group_by(treatment, gene_set_name, controls, mediator) %>%
  # mutate(p_no = n()) %>% 
  # slice(which.min(pcmin)) %>% 
  # ungroup


exB_data = with %>% mutate(treatment= case_when(treatment == "edu_p" ~ "Parental Education",
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
         gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI")) %>% 
  mutate(
    # gene_set_name = gene_set_name %>% str_c(" ","PC"),
    adjP = adjP %>% as.numeric(),
    gene_set_name = gene_set_name %>% str_c(" ","PC_", result_id),
    pval2 = case_when(adjP<0.0001 ~ 100000,
                      adjP<0.001 ~ 25000,
                      adjP<0.01 ~ 15000,
                      adjP<0.05 ~ 10000,
                      adjP>0.05 ~ 0.0000001)
  ) 
if(0) {
  med = mediators[1]
  fig =
    exB_data %>% 
    filter(mediator == med) %>% 
    ggplot(aes(treatment, gene_set_name, size = pval2, )) +
    geom_point(stroke = 1.5, shape = 21, alpha = 0.4, fill = "darkblue", colour = "darkblue") +
    scale_color_manual(values = c("darkblue", "goldenrod3")) + 
    theme_bw() +
    labs(
      # title = "Figure 1. Associations between Indicators of Socioeconomic Status 
      #         and mRNA-Based Disease Signatures, Add Health 
      #         (p-values reported, FDR-corrected for whole genome)",
      y = "mRNA Signatures",
      x = "SES Indicators") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          text = element_text(size=10, face = "bold")) +
    scale_size_continuous(name = "P-value",
                          range = c(0, 16),
                          limits = c(0.0000001, 100000), breaks = c(10000, 15000, 25000, 100000),
                          labels = c("p<0.05", "p<0.01", "p<0.001", "p<0.0001"))+
    scale_alpha(guide = 'none') +
    geom_text(aes(label = proportion), vjust = 0.5, colour = 'black', size = 4) +
    guides(#shape = guide_legend(override.aes = list(size = 10)),
      fill = guide_legend(override.aes = list(size = 10)))+
    scale_fill_manual(values = c("darkblue", "goldenrod3"))
  pdf(str_c("sesdenovoPCAmed_", med,".pdf"))
  print(fig)
  dev.off()
  
  
}