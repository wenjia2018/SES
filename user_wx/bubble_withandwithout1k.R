
#' #### Panle A
#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(tidyverse)
library(here)
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase)
library(enrichplot)
library(gridExtra)
library(grid)
library(ggpubr)
library(dbr) # my package
walk(dir(path = here("R"),full.names = TRUE), source)
example0_with1k = readRDS("/data/home2/xu/old/scratch/xu/example0_new2signature_withinflam.rds")
example0_without1k = readRDS("/home/share/scratch/xu/example0_new2signature_noinflam.rds")
# some old dir have been removed
# example0_with1k<- readRDS("/home/share/scratch/example0_with_1KI.rds")
# example0_without1k <- readRDS("/home/share/scratch/example0_without_1KI.rds")
threshold_with1k = 0.05/10/5
gene_set = c("1KI",
             "Alzheimers","Aortic Aneurysm", "Asthma", 
             "CKD", "COPD", "CVD",
             "Depression", "Diabetes",
             "Hypertension",
             "Rheumatoid Arthritis"
)

ex0_with1k <- example0_with1k %>%
  hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
  filter(controls=="all") %>% 
  dplyr::select(treatment, gene_set_name, p) %>% 
  filter(p<0.05) %>% 
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
         gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI"),
         "1KI Genes" = "With 1KI Genes" %>% as.factor())

ex0_without1k <- example0_without1k %>%
  hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
  filter(controls=="all") %>% 
  dplyr::select(treatment, gene_set_name, p) %>% 
  filter(p<0.05) %>% 
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
         gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI"),
         "1KI Genes" = "Without 1KI Genes"%>% as.factor())

ex_0 = bind_rows(ex0_with1k, ex0_without1k) %>% filter(gene_set_name %in% gene_set)
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
ex_0$gene_set_name <- factor(ex_0$gene_set_name, levels = c("1KI",
                                                            "Alzheimers","Aortic Aneurysm", "Asthma", 
                                                            "CKD", "COPD", "CVD",
                                                            "Depression", "Diabetes",
                                                            "Hypertension",
                                                            "Rheumatoid Arthritis"
                                                             )) %>% fct_rev


#+ echo=F, eval=F, warning=FALSE, message=FALSE

with1k = ex0_with1k %>% mutate(p_withinflam = p) %>% dplyr::select(1,2,7)
without1k = ex0_without1k %>% mutate(p_withoutinflam = p) %>% dplyr::select(1,2,7)
  with1k %>% left_join(without1k, by.x = treatment, by.y = gene_set_name) %>% 
    kableExtra::kable() %>% kableExtra::kable_styling()

#' #### Panle B, pca dimention is totally different things before and after removing inflam1k singatures, so I remove the dimension information
#' and only keep the signature names
  
#+ echo=F, eval=T, warning=FALSE, message=FALSE
  
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
    mutate("1KI Genes" = "With 1KI Genes") %>% 
    mutate(gene_set_name = gene_set_name %>% str_remove("_d.*$") %>% str_trim,
           gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Cvd", "CVD"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Ckd", "CKD"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Copd", "COPD"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI"))
  
  exB_data_with1k %>% kableExtra::kable() %>% kableExtra::kable_styling()
  
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
    mutate("1KI Genes" = "Without 1KI Genes") %>% 
    mutate(gene_set_name = gene_set_name %>% str_remove("_d.*$") %>% str_trim,
           gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Cvd", "CVD"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Ckd", "CKD"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Copd", "COPD"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI"))
  
  exB_data_without1k %>% kableExtra::kable() %>% kableExtra::kable_styling()
  
  exB_data = full_join(exB_data_with1k %>% dplyr::select(1,2,3,6),
                       exB_data_without1k %>% dplyr::select(1,2,3,6),
                       by.x = treatment, by.y = "gene_set_name") %>% 
    filter(gene_set_name %in% gene_set) %>% 
    mutate(gene_set_name = gene_set_name %>% str_c(" ","PC"),
           pval2 = case_when(p<0.0001 ~ 100000,
                           p<0.001 ~ 25000,
                           p<0.01 ~ 15000,
                           p<0.05 ~ 10000,
                           p>0.05 ~ 0.0000001),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "COPD PC", "COPD PC1"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "1KI PC", "1KI PC1"),
           gene_set_name = factor(gene_set_name, levels = c("1KI PC1", "1KI PC2", "Alzheimers PC",
                                                            "Aortic Aneurysm PC", "Asthma PC", "CKD PC",
                                                            "COPD PC1","COPD PC2", "CVD PC",
                                                            "Depression PC",
                                                            "Diabetes PC",
                                                            "Hypertension PC", "Rheumatoid Arthritis PC"
                                                            
           )) %>% fct_rev
    ) 
  # plot pca dimention seprately, for depression and 1KI there are two dimentions significant
 
  # exB_data[8,2] = "COPD PC1"
  exB_data[9,2] = "COPD PC2"
  # exB_data[11,2] = "1KI PC1"
  exB_data[12,2] = "1KI PC2"
  # exB_data[24,2] = "Depression PC1"
  # exB_data[29,2] = "Depression PC2"
  
  
#' ### well loaded genes for significant pc with inflame

#+ echo=F, eval=T, warning=FALSE, message=FALSE 
if(0){
  example0_noerror = remove_errors(example0_with1k)
  example0_m7_ob = example0_noerror %>% get_sig_PCs_and_sig_enrichment_on_those_PCs("m7_ob", threshold = threshold_with1k)
  m7_ob = example0_m7_ob %>%
    filter(map_lgl(well_loaded_genes_on_significant_PCs, ~ length(.x)!=0))%>% 
    dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5"))
  m7_ob %>% saveRDS("./user_wx/m7_ob.rds")
  m7_ob %>% saveRDS("./user_wx/m7_ob_fullpathways.rds")
  }
  
  m7_ob = readRDS("./user_wx/m7_ob.rds")
  a = m7_ob %>%
    dplyr::select(well_loaded_genes_on_significant_PCs) %>%
    pull(well_loaded_genes_on_significant_PCs) 
    # 
  list(Copd1 = a[[8]][["TC1"]], Copd2 = a[[8]][["TC7"]],
       "1KI1" = a[[10]][["TC8"]], "1KI2" = a[[10]][["TC7"]],
       Depression = a[[2]][["TC6"]],
       Arthritis = a[[6]][["TC9"]],
       Asthma = a[[12]][["TC6"]]) %>% openxlsx::write.xlsx("./user_wx/well_loaded_gene_withinflam.xlsx")
  
  all_pathway = m7_ob %>%
    unnest_longer(enrichment_of_well_loaded_genes_on_significant_PCs) %>%
    hoist(enrichment_of_well_loaded_genes_on_significant_PCs, pathway = list("out","enriched_physiology","reactome")) 
  
  names = all_pathway$gene_set_name
  
  temp_with = all_pathway$pathway %>%
    map(function(x) as_tibble(x) %>%
          dplyr::select(Description, p.adjust) #%>%
          # dplyr::top_n(-4)
        ) #%>%
    # set_names(all_pathway$gene_set_name) %>% 
    # .[c(1,2,5,6,14)] %>% 
    # map_df(~as.data.frame(.x), .id="Signature") %>%
    # as_tibble() %>% 
    # unique() %>% 
    # .[c(1,2,3,4,5,7,8,9,10,11,12,14,15,16),]
  # depression
  # depr=temp_with[c(2,4,7,10,16)] %>% bind_rows() %>% dplyr::distinct(Description, .keep_all = TRUE)
  depr = temp_with[c(2,4,7,10,16)] %>%
    bind_rows() %>% 
    arrange(Description, p.adjust) %>% 
    dplyr::distinct(Description, .keep_all = TRUE) %>% 
    arrange(p.adjust) %>% 
    mutate_at(.vars =c("p.adjust"),
              .funs = list(~ .x %>% format(digits = 3, scientific =T))) 
    
  
  depr %>% kableExtra::kable() %>% kableExtra::kable_styling()
  # copd
  copd = temp_with[c(1,3,8,9,15)] %>%
    bind_rows() %>%
    arrange(Description, p.adjust) %>% 
    dplyr::distinct(Description, .keep_all = TRUE) %>% 
    arrange(p.adjust) %>% 
    mutate_at(.vars =c("p.adjust"),
              .funs = list(~ .x %>% format(digits = 3, scientific =T))) 
  
  copd %>% kableExtra::kable() %>% kableExtra::kable_styling()
  # inflam1k
  i1k=temp_with[c(5,11,12)] %>%
    bind_rows() %>%
    arrange(Description, p.adjust) %>% 
    dplyr::distinct(Description, .keep_all = TRUE) %>% 
    arrange(p.adjust) %>% 
    mutate_at(.vars =c("p.adjust"),
              .funs = list(~ .x %>% format(digits = 3, scientific =T))) 
  
  i1k %>% kableExtra::kable() %>% kableExtra::kable_styling()
  
  
  
  
  # temp_with = rbind(depr,copd,i1k)
  # rownames(temp_with) <- NULL
  # 

#' without 1KI
#+ echo=F, eval=F, warning=FALSE, message=FALSE 

    example0_noerror_without1k = remove_errors(example0_without1k)
  example0_m7_ob = example0_noerror_without1k %>%
    get_sig_PCs_and_sig_enrichment_on_those_PCs("m7_ob",
                                                threshold = threshold_with1k)
  m7_ob = example0_m7_ob %>%
    filter(map_lgl(well_loaded_genes_on_significant_PCs, ~ length(.x)!=0))%>%
    dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5"))
  m7_ob %>%  saveRDS("./user_wx/m7_ob_without1KI.rds")
  m7_ob %>%  saveRDS("./user_wx/m7_ob_without1KI_fullpathways.rds")
  m7_ob = readRDS("./user_wx/m7_ob_without1KI.rds")
  b = m7_ob %>%
    dplyr::select(well_loaded_genes_on_significant_PCs) %>%
    pull(well_loaded_genes_on_significant_PCs) %>%
    setNames(m7_ob$gene_set_name)
  
  # b[c(1,3,4,5,8)] %>% openxlsx::write.xlsx("./user_wx/well_loaded_gene_noinflam.xlsx")
  
  all_pathway = m7_ob %>%
    unnest_longer(enrichment_of_well_loaded_genes_on_significant_PCs) %>%
    hoist(enrichment_of_well_loaded_genes_on_significant_PCs, pathway = list("out","enriched_physiology","reactome")) 
  names = all_pathway$gene_set_name
  temp_without = all_pathway$pathway %>%
    map(function(x) as_tibble(x) %>% dplyr::select(Description, p.adjust))
    #     %>% dplyr::top_n(-2)) %>%
    # set_names(all_pathway$gene_set_name) %>% 
    # .[c(1,3,4,5,8)] %>% 
    # map_df(~as.data.frame(.x), .id="Signature") %>%
    # as_tibble() %>% 
    # unique()
  
  depr = temp_without[c(3,7,12)] %>%
    bind_rows() %>% 
    arrange(Description, p.adjust) %>% 
    dplyr::distinct(Description, .keep_all = TRUE) %>% 
    arrange(p.adjust) %>% 
    mutate_at(.vars =c("p.adjust"),
              .funs = list(~ .x %>% format(digits = 3, scientific =T))) 
  
  
  depr %>% kableExtra::kable() %>% kableExtra::kable_styling()
  # copd
  copd = temp_without[c(1,2,6,11)] %>%
    bind_rows() %>%
    arrange(Description, p.adjust) %>% 
    dplyr::distinct(Description, .keep_all = TRUE) %>% 
    arrange(p.adjust) %>% 
    mutate_at(.vars =c("p.adjust"),
              .funs = list(~ .x %>% format(digits = 3, scientific =T))) 
  
  copd %>% kableExtra::kable() %>% kableExtra::kable_styling()
  # inflam1k
  arthritis = temp_without[c(4,9,13)] %>%
    bind_rows() %>%
    arrange(Description, p.adjust) %>% 
    dplyr::distinct(Description, .keep_all = TRUE) %>% 
    arrange(p.adjust) %>% 
    mutate_at(.vars =c("p.adjust"),
              .funs = list(~ .x %>% format(digits = 3, scientific =T))) 
  
  arthritis %>% kableExtra::kable() %>% kableExtra::kable_styling()
  
  
  hypertension = temp_without[c(8)] %>%
    bind_rows() %>%
    arrange(Description, p.adjust) %>% 
    dplyr::distinct(Description, .keep_all = TRUE) %>% 
    arrange(p.adjust) %>% 
    mutate_at(.vars =c("p.adjust"),
              .funs = list(~ .x %>% format(digits = 3, scientific =T))) 
  
  hypertension %>% kableExtra::kable() %>% kableExtra::kable_styling()
  
  rownames(temp_without) <- NULL
  # check if pathways are equal for with and without inflamation
  all.equal(temp_with, temp_without)

  
#+ echo=F, eval=T, warning=FALSE, message=FALSE   
 #  temp = temp_with %>% mutate_at(.vars =c("p.adjust"), .funs = list(~ .x %>% format(digits = 3, scientific =T)))
 # temp %>% kableExtra::kable() %>% kableExtra::kable_styling()
   mytheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 1),
                bg_params = list(fill="white", col=NA)),
    colhead = list(fg_params=list(cex = 0.7),
                   bg_params = list(fill="white", col=NA)),
    padding = unit(c(0.4, 4), "mm"),
    rowhead = list(fg_params=list(cex = 0.5)))
   
  # panelT1 = tableGrob(temp[1:5, ] %>% dplyr::select(2,3), rows = NULL, theme = mytheme)
  # panelT2 = tableGrob(temp[6:10, ] %>% dplyr::select(2,3), rows = NULL, theme = mytheme)
  # panelT3 = tableGrob(temp[11:14, ] %>% dplyr::select(2,3), rows = NULL, theme = mytheme)
   panelT1 = tableGrob(copd[1:5, ], rows = NULL, theme = mytheme)
   panelT2 = tableGrob(depr[1:5, ], rows = NULL, theme = mytheme)
   panelT3 = tableGrob(i1k[1:5, ], rows = NULL, theme = mytheme)
   
  dd1 <- ggplot() + annotation_custom(panelT1) + labs(title = 'COPD') + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, vjust= -6, face = "bold", color = "darkblue"))
  dd2 <- ggplot() + annotation_custom(panelT2) + labs(title = 'Depression')+ theme_void() +
    theme(plot.title = element_text(hjust = 0.5, vjust= -6, face = "bold", color = "darkblue"))
  dd3 <- ggplot() + annotation_custom(panelT3) + labs(title = '1KI')+ theme_void() +
    theme(plot.title = element_text(hjust = 0.5, vjust= -6, face = "bold", color = "darkblue"))
  
  # grid.arrange(bottom = "cv",panelT1)
  
  panelC = grid.arrange(dd1, dd2, dd3, ncol = 3)
  ggpubr::ggarrange(
    panelC,
    labels = c("C"),
    font.label = list(size = 20)
    # Label of the line plot
  )
  
  ggpubr::ggarrange(panelA, panelB, ncol = 2, labels = c("A", "B"), common.legend = TRUE, legend="right")
  ggpubr::ggarrange(
    panelC,
    font.label = list(size = 20),
    labels = c("C")       # Label of the line plot
  )
  # panelC = grid.arrange(panelT1, panelT2, panelT3, nrow =1)
  panelA = ggplot(ex_0, aes(treatment, gene_set_name, size = pval2,
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
  
  # color order from bottom to top
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
                          limits = c(0.0000001, 100000), breaks = c(15000, 25000),
                          labels = c("p<0.005", "p<0.001"))+
    scale_alpha(guide = 'none')+
    guides(shape = guide_legend(override.aes = list(size = 10)),
           fill = guide_legend(override.aes = list(size = 8)))
  ggpubr::ggarrange(panelA, panelB, ncol = 2, labels = c("A", "B"), common.legend = TRUE, legend="right")
  figure = ggpubr::ggarrange(
    ggpubr::ggarrange(panelA, panelB, ncol = 2, labels = c("A", "B"), common.legend = TRUE, legend="right"), 
    panelC,
    nrow = 2,
    heights = c(3, 1),
    labels = c("","C")       # Label of the line plot
  )
  
  title <- expression(atop("Figure 1:",
                           scriptscriptstyle("Panel A:Associations between Indicators of Socioeconomic Status and mRNA-Based Disease Signatures, Add Health (p-values reported, FDR-corrected for whole genome),Panel B: Significant PCs,Panel C: pathways")))

  
  # title <- expression(atop("Figure 1:",
  #                          scriptscriptstyle("Panel A:Associations between Indicators of Socioeconomic Status and mRNA-Based Disease Signatures, Add Health (p-values reported, FDR-corrected for whole genome),Panel B: Significant PCs,Panel C: pathways")))
  
  # title <- expression(atop(bold("Figure 1:"),
  # atop("panel A:Associations between Indicators of Socioeconomic Status and mRNA-Based Disease Signatures, Add Health (p-values reported, FDR-corrected for whole genome) \n panel B: Significant PCs \n panel C: pathways")))

  annotate_figure(figure, bottom = text_grob(title))
  annotate_figure(figure, fig.lab = "Figure 1:",fig.lab.pos = "bottom")
  
#'`rmarkdown::render("/home/xu/ses-1/user_wx/bubble_withandwithout1k.R")`

  
