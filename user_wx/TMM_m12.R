
#' ---
#' title: Omnibus test results for ses using TMM normalization
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->


#' ###  regression of gene on SES indicator with basic controls
#' * gene ~ ses + controls
#' * about 12,000 genes(differs in each regression)
#' * p values are genowide FDR corrected
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
library(tidyverse)
# TABLE 1
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

source("/home/xu/ses-1/user_wx/utils_logfc.R")
#' ###  intersection of DE genes for each SES indicator
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
example_TMM_DE <- readRDS("~/ses-1/user_wx/example_tmm_genowidebydesign.rds")
data_DE = example_TMM_DE %>% 
  hoist(out, ttT = list("result", "ttT")) %>% 
  filter(map_lgl(.$ttT, ~dim(.)[1]!=0)) %>% 
  mutate(gene_sig = map(ttT, ~ filter(., adj.P.Val<0.05) %>% pull(gene)) %>% set_names(treatment)) %>% 
  # mutate(min_p = map(ttT, ~ dplyr::slice_min(.,adj.P.Val) %>% pull(adj.P.Val)) %>% set_names(treatment)) %>% 
  dplyr::select(-out, -gene_set_name) 

venn::venn(data_DE$gene_sig, zcolor = "style",
           ellipse = FALSE, opacity = 0.15, ilcs = 1,
           box = FALSE,
           sncs = 1)

UpSetR::upset(UpSetR::fromList(data_DE$gene_sig), 
              matrix.color ="darkred", 
              main.bar.color = "gray23",
              sets.bar.color=c("navajowhite4","lightsteelblue1","orange","goldenrod3","cadetblue1"),
              order.by = "degree")

#' ### functional pathway results
#' * do pathway analysis using DE gene and their rank in enrichr website
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
if(0){
  data_DE$ttT %>% 
    map(~ filter(.x, adj.P.Val<0.05)) %>% 
    set_names(data_DE$treatment) %>% 
    openxlsx::write.xlsx("./user_wx/TMM_DE_ses.xlsx")
  
  
  
  for (i in 1:5) {
    temp = 
      data_DE$ttT %>% 
      map(~ filter(.x, adj.P.Val<0.05 & logFC<0) %>% pull(gene)) %>%
      set_names(data_DE$treatment) 
    
    name = names(temp)[i]
    write(temp[[i]], file = str_c(getwd(),"/user_wx/", name, "_down.txt"))
  }
  for (i in 1:5) {
    temp = 
      data_DE$ttT %>% 
      map(~ filter(.x, adj.P.Val<0.05 & logFC>0) %>% pull(gene)) %>%
      set_names(data_DE$treatment) 
    
    name = names(temp)[i]
    write(temp[[i]], file = str_c(getwd(),"/user_wx/", name, "_up.txt"))
  }
}

path = read.delim("~/ses-1/user_wx/KEGG_2021_Human_table.txt", header = TRUE, sep = "\t", dec = ".")

path = path %>% filter(Adjusted.P.value <0.05) %>% dplyr::select(-5,-6)


# direction = path$Genes %>% map(~ str_split(.x, ";")) %>% flatten()

path = path %>% mutate(gene_list = map(.$Genes, ~str_split(.x, ";") %>% unlist),
                       `direction(-, +)`  = gene_list %>% map( ~ filter(data_DE$ttT[[4]], gene %in% .) %>% pull("logFC") %>% sign %>% table))
path %>% dplyr::select(-Genes, -gene_list) %>%  kableExtra::kable() %>% kableExtra::kable_styling()
#' ###  logFC plot for each treatment and each disease signature
#' * test are correlated t test at a correlation level of 0.1
#' * test hypothesis is the mean of the logFCs is 0
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
example0_with1k <- readRDS("~/ses-1/user_wx/example_tmm_m12_withinflame.rds")
# supplement include people with disease
# example0_with1k <- readRDS("~/ses-1/user_wx/example_tmm_m12_withinflame_withdiseasepeople.rds")
data_logfc_with1k = example0_with1k %>%  hoist(out, m = list("result",  "m12_fdr", 1, "other", "m"))

data_logfc_with1k %>% DE_logFC_ploting_notest(caption_text = "withinflamation 1K genes")



example0_without1k  <- readRDS("~/ses-1/user_wx/example_tmm_m12_noinflame.rds")
# supplement include people with disease
# example0_without1k  <- readRDS("~/ses-1/user_wx/example_tmm_m12_noinflame_withdiseasepeople.rds")
data_logfc = example0_without1k %>%  hoist(out, m = list("result",  "m12_fdr", 1, "other", "m"))

data_logfc %>% DE_logFC_ploting(caption_text = "removing inflamation 1K genes")

#' ###  bubble polot of fig1 panelA for each SES indicator
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA

ex0_without1k <- example0_without1k %>%
  hoist(out, p = list("result", "m12_fdr", 1, "p")) %>% 
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


ex0_with1k <- example0_with1k %>%
  hoist(out, p = list("result", "m12_fdr", 1, "p")) %>% 
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
  pdf("myplot.pdf")
  print(fig)
  dev.off()
  
}

outcome_set = readRDS("/home/share/preprocessing/preprocessed_two_batches/allpossiblegene.rds")

outcome_set = outcome_set$outcome_set[table1]



temp = NULL
for(i in 1: (data_DE$treatment %>% length)) {
  temp[[i]] = outcome_set %>% map(~ .x %in% data_DE$gene_sig[[i]] %>% sum) %>% bind_rows
}
temp = temp %>% bind_rows() %>% t %>% `colnames<-`(data_DE$treatment)


outcome_set_noinflame = outcome_set %>% map(~ setdiff(.x, outcome_set$inflam1k_mRNA))

temp_noinflame = NULL
for(i in 1: (data_DE$treatment %>% length)) {
  temp_noinflame[[i]] = outcome_set_noinflame %>% map(~ .x %in% data_DE$gene_sig[[i]] %>% sum) %>% bind_rows
}
temp_noinflame = temp_noinflame %>% bind_rows() %>% t %>% `colnames<-`(data_DE$treatment)

#' ###  Number DE genes in each predefined signature for each SES indicator
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA


temp %>% kableExtra::kable() %>% kableExtra::kable_styling()

#' ### Number DE genes in each predefined signature(excluding inflammation 1k) for each SES indicator
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA



temp_noinflame %>% kableExtra::kable() %>% kableExtra::kable_styling()

#' ### Analog to fig1 panleA in pnas paper 
#' * BUT here the size of bubbles means the number of significant genes in this signature
#' * without cell type controls
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA


temp %>% 
  as.data.frame() %>% 
  rownames_to_column("signature") %>% 
  pivot_longer(cols = (2:6), names_to = "treatment", values_to ="No.SigGene") %>% 
  mutate(No.SigGene = ifelse(No.SigGene==0, NA, No.SigGene),
         `1KIgene` = "with 1KI") %>% 
  bind_rows(temp_noinflame %>% 
              as.data.frame() %>% 
              rownames_to_column("signature") %>% 
              pivot_longer(cols = (2:6), names_to = "treatment", values_to ="No.SigGene") %>% 
              mutate(No.SigGene = ifelse(No.SigGene==0, NA, No.SigGene),
                     `1KIgene` = "without 1KI")) %>% 
  ggplot(mapping = aes(x = treatment, size = No.SigGene, y = signature, fill = `1KIgene`)) + 
  geom_point( shape = 21, alpha = 0.4) +
  scale_fill_manual(values = c("darkblue", "goldenrod3")) +
  scale_color_manual(values = c("darkblue", "goldenrod3")) +
  scale_size_continuous(breaks = c(10, 20, 50, 100)) 

