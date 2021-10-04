#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
library(tidyverse)
threshold = 0.05
nfactors = 10
# rle
# example0_with1k <- readRDS("~/ses-1/user_wx/example_RLE_pca_nomed_withinflame.rds")
# example0_without1k<- readRDS("~/ses-1/user_wx/example_RLE_pca_nomed_noinflame.rds")
# tmm
example0_with1k <-readRDS("~/ses-1/user_wx/example_tmm_m7_withinflame.rds")
example0_without1k <- readRDS("~/ses-1/user_wx/example_tmm_m7_noinflame0209.rds")

plotPC = function(example0_with1k, example0_without1k, rotation, threshold, nfactors) {
  
  exB <- example0_with1k %>%
    hoist(out, p = list("result", rotation, 1, "p")) %>% 
    unnest_longer(p) %>% 
    dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
    mutate(p = p.adjust(p, method = "fdr"),
           pval2=case_when(p<0.0001 ~ 100000,
                           p<0.001 ~ 25000,
                           p<0.01 ~ 15000,
                           p<0.05 ~ 10000,
                           p>0.05 ~ 0.0000001)
    ) %>% 
    dplyr::select(treatment, gene_set_name, p, p_id, pval2) %>% 
    dplyr::filter(p < threshold, p_id %in% c(str_c("d",1:nfactors)) )
  
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
    mutate("1KI Genes" = "With 1KI Genes") %>% 
    mutate(gene_set_name = gene_set_name %>% str_remove("_d.*$") %>% str_trim,
           gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Cvd", "CVD"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Ckd", "CKD"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Copd", "COPD"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI"))
  
  exB_data_with1k %>% kableExtra::kable() %>% kableExtra::kable_styling()
  
  exB <- example0_without1k %>%
    hoist(out, p = list("result", rotation, 1, "p")) %>% 
    unnest_longer(p) %>% 
    dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
    mutate(p = p.adjust(p, method = "fdr"),
           pval2=case_when(p<0.0001 ~ 100000,
                           p<0.001 ~ 25000,
                           p<0.01 ~ 15000,
                           p<0.05 ~ 10000,
                           p>0.05 ~ 0.0000001)
    ) %>% 
    dplyr::select(treatment, gene_set_name, p, p_id, pval2) %>% 
    dplyr::filter(p < threshold, p_id %in% c(str_c("d",1:nfactors)) )
  
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
    mutate("1KI Genes" = "Without 1KI Genes") %>% 
    mutate(gene_set_name = gene_set_name %>% str_remove("_d.*$") %>% str_trim,
           gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Cvd", "CVD"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Ckd", "CKD"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Copd", "COPD"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI"))
  
  exB_data_without1k %>% kableExtra::kable() %>% kableExtra::kable_styling()
  
  exB_data = full_join(exB_data_with1k %>% mutate(gene_set_name = str_c(gene_set_name,"_", p_id)) %>% dplyr::select(1,2,3,6),
                       exB_data_without1k %>% mutate(gene_set_name = str_c(gene_set_name,"_", p_id)) %>% dplyr::select(1,2,3,6),
                       by.x = treatment, by.y = "gene_set_name") %>% 
    # filter(gene_set_name %in% gene_set) %>% 
    mutate(
      # gene_set_name = gene_set_name %>% str_c(" ","PC"),
      pval2 = case_when(p<0.0001 ~ 100000,
                        p<0.001 ~ 25000,
                        p<0.01 ~ 15000,
                        p<0.05 ~ 10000,
                        p>0.05 ~ 0.0000001)
      # gene_set_name = gene_set_name %>% replace(gene_set_name == "COPD PC", "COPD PC1"),
      # gene_set_name = gene_set_name %>% replace(gene_set_name == "1KI PC", "1KI PC1"),
      # gene_set_name = factor(gene_set_name, levels = c("1KI PC1", "1KI PC2", "Alzheimers PC",
      #                                                  "Aortic Aneurysm PC", "Asthma PC", "CKD PC",
      #                                                  "COPD PC1","COPD PC2", "CVD PC",
      #                                                  "Depression PC",
      #                                                  "Diabetes PC",
      #                                                  "Hypertension PC", "Rheumatoid Arthritis PC"
      #                                                  
      # )) %>% fct_rev
    ) 

  
  axiscolor = c("grey30", "grey30", "grey30", "darkseagreen1","darkseagreen4" ,"grey30","red", "darkred")
  panelB = ggplot(exB_data, aes(treatment, gene_set_name, size = pval2, fill =`1KI Genes`, colour = `1KI Genes`)) +
    geom_point(stroke = 1.5, shape = 21, alpha = 0.4) +
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
           fill = guide_legend(override.aes = list(size = 8))) +
    labs(caption = str_c("rotation :", rotation, ",  threshold :", threshold, ", nfactors:", nfactors))
  panelB

} 
if(plot<-FALSE){
  pdf("sesPCA.pdf")
  print(panelB)
  dev.off()
  
}
#' ### 10 PC oblimin
#+ echo=F, eval=T, warning=FALSE, message=FALSE
example0_with1k <-readRDS("~/ses-1/user_wx/example_tmm_m7_withinflame.rds")
example0_without1k <- readRDS("~/ses-1/user_wx/example_tmm_m7_noinflame0209.rds")

plotPC(example0_with1k, example0_without1k, rotation = "m7_ob", threshold = 0.05, nfactors = 10)
plotPC(example0_with1k, example0_without1k, rotation = "m7_ob", threshold = 0.05, nfactors = 5)
plotPC(example0_with1k, example0_without1k, rotation = "m7_ob",threshold = 0.01, nfactors = 5)

#' ### 10 PC no rotation
#+ echo=F, eval=T, warning=FALSE, message=FALSE
example0_with1k <-readRDS("~/ses-1/user_wx/example_tmm_m7_withinflame_nonerotation.rds")
example0_without1k <- readRDS("~/ses-1/user_wx/example_tmm_m7_noinflame_nonerotation.rds")

plotPC(example0_with1k, example0_without1k, rotation = "m7_nn", threshold = 0.05, nfactors = 10)
plotPC(example0_with1k, example0_without1k, rotation = "m7_nn", threshold = 0.05, nfactors = 5)
plotPC(example0_with1k, example0_without1k, rotation = "m7_nn", threshold = 0.01, nfactors = 5)

#' ### 5 PC oblimin
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example0_with1k <-readRDS("~/ses-1/user_wx/example_tmm_m7_withinflame_oblimin5factors.rds")
example0_without1k <- readRDS("~/ses-1/user_wx/example_tmm_m7_noinflame_oblimin5factors.rds")


plotPC(example0_with1k, example0_without1k, rotation = "m7_ob", threshold = 0.05, nfactors = 5)
plotPC(example0_with1k, example0_without1k, rotation = "m7_ob",threshold = 0.01, nfactors = 5)
