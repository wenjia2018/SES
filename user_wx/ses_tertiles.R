#' ###  bubble polot of fig1 panelA for each SES indicator
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
# ses tertiles
example0_with1k <- readRDS("~/ses-1/user_wx/example_tmm_m12_withinflame_ses3.rds")

# ses tertiles
example0_without1k <- readRDS("~/ses-1/user_wx/example_tmm_m12_noinflame_ses3.rds")



ex0_without1k <- example0_without1k %>%
  hoist(out, p = list("result", "m12_fdr", 1, "p")) %>% 
  dplyr::select(treatment, gene_set_name, p) %>% 
  filter(p<0.05) %>% 
  # dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  # dplyr::filter(treatment %in% c("edu_p", "income_pp1_log", "SEI_max_p_w12", "ses_composite_pp1" )) %>% 
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
                              treatment =="ses_composite_ff5"  ~ "SES Composite 3",
                              treatment =="ses_sss_composite_3_1" ~ "SES Composite Low",
                              treatment =="ses_sss_composite_3_3"  ~ "SES Composite High"
                              
                              ))  %>%
  mutate(treatment = factor(treatment, levels = c("Parental SES Composite", "Parental Education","Parental Income",
                                                  "Parental SEI", "Mother's Occupation", "Father's Occupation",
                                                  "SES Composite 3", "SES Composite", "Education", "Income",
                                                  "Occupation", "Subjective Social Status","Occupation Work Collar",
                                                  "SES Composite Low","SES Composite High"
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
  dplyr::select(treatment, gene_set_name, p) %>% 
  filter(p<0.05) %>% 
  # dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  # dplyr::filter(treatment %in% c("edu_p", "income_pp1_log", "SEI_max_p_w12", "ses_composite_pp1" )) %>%
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
                              treatment =="ses_composite_ff5"  ~ "SES Composite 3",
                              treatment =="ses_sss_composite_3_1" ~ "SES Composite Low",
                              treatment =="ses_sss_composite_3_3"  ~ "SES Composite High"
                              
         ))  %>%
  mutate(treatment = factor(treatment, levels = c("Parental SES Composite", "Parental Education","Parental Income",
                                                  "Parental SEI", "Mother's Occupation", "Father's Occupation",
                                                  "SES Composite 3", "SES Composite", "Education", "Income",
                                                  "Occupation", "Subjective Social Status","Occupation Work Collar",
                                                  "SES Composite Low","SES Composite High"
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


fig =
  ggplot(ex_0, aes(treatment, gene_set_name, size = pval2,
                   fill = `1KI Genes`,
                   colour = `1KI Genes`)) +
  geom_point(stroke = 1.5, shape = 21, alpha = 0.4) +
  scale_fill_manual(values = c("darkblue", "goldenrod3")) +
  scale_color_manual(values = c("darkblue", "goldenrod3")) +
  # geom_jitter(height = 0.00000025) +
  # gghighlight(class == "inflam") +
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

fig
if(plot<-FALSE){
  pdf("ses_tertiles.pdf")
  print(fig)
  dev.off()
}
DE_logFC_ploting_ses3= function(example, caption_text){
  a = example %>% 
    # mutate(dim = map_int(.$m, ~dim(.)[1]),
    #        gene_set_name = str_c(gene_set_name, "(", dim, ")")) %>% 
    ungroup %>% 
    unnest(m) %>% 
    dplyr::select(treatment, gene_set_name, logFC) 
  
  axiscolor = c("grey30", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30","grey30", "grey30", "grey30")
  a %>% mutate(
    treatment= case_when(treatment =="ses_sss_composite_3_1" ~ "SES Composite Low",
                         treatment =="ses_sss_composite_3_3"  ~ "SES Composite High" 
  )) %>% 
    ggplot(aes(x = logFC, y = gene_set_name, fill = treatment, color = treatment)) +
    geom_violin()+
    facet_wrap( ~  treatment , scales = "free_x", strip.position = "bottom") +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    scale_color_manual(values=c("goldenrod3", "lightblue")) +
    scale_fill_manual(values=c("goldenrod3", "lightblue")) +
    theme(
      # legend.position = "none",
      panel.margin.x = unit(1, "lines"),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(color=axiscolor),
      axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
      axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
      plot.title = element_text(size = 20, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
      plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
      strip.text = element_text(size = 7),
      text = element_text(family = "Georgia"))+
    labs(y = "mRNA Signatures",
         caption = paste(caption_text))
  
}

data_logfc = example0_without1k %>%  hoist(out, m = list("result",  "m12_fdr", 1, "other", "m"))

fig1 = data_logfc %>% DE_logFC_ploting_ses3(caption_text = "")


fig1
if(plot<-FALSE){
  pdf("ses_logfc_tertiles.pdf")
  print(fig1)
  dev.off()
}
