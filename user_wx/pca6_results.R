#+ echo=F, eval=T

set.seed(123)
library(here)
library(tidyverse)
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase)
library(enrichplot)
library(dbr) # my package

walk(dir(path = here("R"),full.names = TRUE), source)
load_data(reconciled = FALSE, remove_inflam = TRUE)
define_treatments_and_controls()
recode_variables_in_dat()
print(abbreviations)
funcs = str_subset(abbreviations$shorthand, "^m") 
funcs = funcs %>% str_subset("m[6-8]")

############################################################
# EXAMPLE: SIGNATURES
############################################################
threshold = 0.05/24
    
    example0 = readRDS("/home/share/scratch/xu/example0_w5bmi_removeinflam_6pc_trial2.rds")
  
  ex0 <- example0 %>%
    hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
    filter(controls=="all") %>% 
    dplyr::select(treatment, gene_set_name, p)
  
  ex_0 <- ex0 %>%
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
                                treatment =="work_collar_ff5" ~ "Occupation",
                                treatment =="edu_max" ~ "Education" ,
                                treatment =="income_hh_ff5" ~ "Income"     ,
                                treatment =="SEI_ff5" ~ "SEI"      ,
                                treatment =="ses_sss_composite" ~"SES Composite 4"  ,
                                treatment =="sss_5" ~ "SSS",
                                treatment =="ses_composite_ff5"  ~"SES Composite 3"))  %>%
    mutate(treatment = factor(treatment, levels = c("Parental SES Composite", "Parental Education","Parental Income",
                                                    "Parental SEI", "Mother's Occupation", "Father's Occupation",
                                                    "SES Composite 3", "SES Composite 4", "Education", "Income",
                                                    "SEI", "SSS","Occupation"
    )))
  
  
  panelA = ggplot(ex_0, aes(treatment, gene_set_name, size = pval2, alpha = 0.4)) +
    geom_point(fill = "red", color="navy") +
    theme_bw() +
    labs(
      # title = "Figure 1. Associations between Indicators of Socioeconomic Status 
      #         and mRNA-Based Disease Signatures, Add Health 
      #         (p-values reported, FDR-corrected for whole genome)",
      y = "mRNA Signatures",
      x = "SES Indicators") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_size_continuous(name = "P-value",
                          range = c(0, 11),
                          limits = c(0.0000001, 100000), breaks = c(0.0000001, 10000, 15000, 25000, 100000),
                          labels = c("n.s.", "p<0.05", "p<0.01", "p<0.001", "p<0.0001"))+
    scale_alpha(guide = 'none')
  
  # panel B

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
    dplyr::filter(p < threshold) 
  
  exB_data = exB %>% 
    mutate(gene_set_name = str_c(gene_set_name,"_",p_id),
           treatment= case_when(treatment == "edu_p" ~ "Parental Education",
                                treatment =="income_pp1_log" ~  "Parental Income" ,
                                treatment =="SEI_max_p_w12" ~ "Parental SEI",
                                treatment =="ses_composite_pp1" ~ "Parental SES Composite",
                                treatment =="work_collar_rm_f12" ~ "Mother's Occupation",
                                treatment =="work_collar_rf_f12" ~ "Father's Occupation" ,
                                treatment =="work_collar_ff5" ~ "Occupation",
                                treatment =="edu_max" ~ "Education" ,
                                treatment =="income_hh_ff5" ~ "Income"     ,
                                treatment =="SEI_ff5" ~ "SEI"      ,
                                treatment =="ses_sss_composite" ~"SES Composite 4"  ,
                                treatment =="sss_5" ~ "SSS",
                                treatment =="ses_composite_ff5"  ~"SES Composite 3"))  %>%
    mutate(treatment = factor(treatment, levels = c("Parental SES Composite", "Parental Education","Parental Income",
                                                    "Parental SEI", "Mother's Occupation", "Father's Occupation",
                                                    "SES Composite 3", "SES Composite 4", "Education", "Income",
                                                    "SEI", "SSS","Occupation")))
  panelB = ggplot(exB_data, aes(treatment, gene_set_name, size = pval2, alpha = 0.4)) +
    geom_point(fill = "red", color="navy") +
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
  
  # ggplot(exB, aes(p_id, pval, alpha = 0.4, colour = treatment, shape = gene_set_name)) +
  #   geom_point() 
  
  
  # panel C interesting pc pathways
  
  
  # ESTIMATE VARIOUS PCA "ROTATIONS"
  example0_noerror = remove_errors(example0)
  example0_m7_ob = example0_noerror %>% get_sig_PCs_and_sig_enrichment_on_those_PCs("m7_ob", threshold = threshold)
  # example0_m7_ob = readRDS("/home/xu/ses-1/user_wx/ob_rotation_pca.rds")
  
  m7_ob = example0_m7_ob %>%
    filter(map_lgl(well_loaded_genes_on_significant_PCs, ~ length(.x)!=0)) %>% 
    filter(treatment != "ses_composite_ff5")
  
  all_pathway = m7_ob %>%
    unnest_longer(enrichment_of_well_loaded_genes_on_significant_PCs) %>%
    hoist(enrichment_of_well_loaded_genes_on_significant_PCs, pathway = list("out","enriched_physiology","reactome")) 
  
  all_pathway$pathway %>%
    map(function(x) as_tibble(x) %>% dplyr::select(Description, p.adjust) %>% dplyr::top_n(-3)) %>%
    set_names(all_pathway$gene_set_name) %>% 
    map_df(~as.data.frame(.x), .id="tag") %>%
    as_tibble() %>%
    unique() %>% 
    kableExtra::kable() %>%
    kableExtra::kable_styling()
  
  rownames(temp) <- NULL
  
  # panelC = ggtexttable(temp %>% dplyr::select(2,3))
  
  
  panelT1 = tableGrob(temp[1:5, ] %>% dplyr::select(2,3), rows = NULL)
  panelT2 = tableGrob(temp[6:10, ] %>% dplyr::select(2,3), rows = NULL)
  panelT3 = tableGrob(temp[11:15, ] %>% dplyr::select(2,3), rows = NULL)
  # panelF = tableGrob(temp[15, ] %>% dplyr::select(2,3), rows = NULL,
  #                   theme = ttheme_default(
  #                   core = list(fg_params=list(cex = 1)),
  #                   colhead = list(fg_params=list(cex = 1)),
  #                   rowhead = list(fg_params=list(cex = 1.0))))
  
  # 
  # text <- paste("Figure 1:","\npanel A: Associations between Indicators of Socioeconomic Status and mRNA-Based Disease Signatures, Add Health (p-values reported, FDR-corrected for whole genome)", "panel B: Significant PCs", "\npanel C: pathways")
  # text.p <- ggparagraph(text = text, face = "bold", size = 11, color = "black")
  # lay <- rbind(c(1,1,1,2,2,2),
  #              c(1,1,1,2,2,2),
  #              c(1,1,1,2,2,2),
  #              c(3,3,3,3,3,3),
  #              c(3,3,3,3,3,3),
  #              c(4,4,4,4,4,4))
  
  panelC = grid.arrange(panelT1, panelT2, panelT3, nrow =1)
  # figure = grid.arrange(panelA, panelB, panelC, text.p,
  #                       layout_matrix = lay
  #                       )
  # ggarrange can label each element
  figure = ggarrange(
    ggarrange(panelA, panelB, ncol = 2, labels = c("A", "B")), 
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
  
  
