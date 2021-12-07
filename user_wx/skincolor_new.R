
#' ---
#' title: Omnibus test results for skin color
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->


#' ###  regression of gene on several skin color indicators(5 levels, 3 levels, continuous) with basic controls in ses paper(remove race ethnicity)
#' * gene ~ skincolor + controls
#' * about 12,000 genes(differs in each regression)
#' * p values are FDR corrected within each signature
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
    "CKD_mRNA")
ctra = c(
    "ctra_mRNA",
    "inflame_mRNA",
    "interferon_mRNA",
    "AntBIntF_mRNA",
    # "antibody_mRNA", #only 1 gene
    "inflam1k_mRNA")
aging = c(
    "aging_mRNA",
    "aging_up_mRNA",
    "aging_down_mRNA",
    "aging_cluster_complement_mRNA",
    "aging_down_cl1_mRNA",
    "aging_down_cl1a_mRNA",
    "aging_down_cl1b_mRNA",
    "aging_down_cl1c_mRNA",
    "aging_down_cl2_mRNA",
    "aging_down_cl3_mRNA",
    "aging_up_cl1_mRNA",
    "aging_up_cl2_mRNA",
    "aging_up_cl3_mRNA",
    "aging_up_cl4_mRNA"
  )
source("/home/xu/ses-1/user_wx/utils_logfc.R")
#' ###  intersection of DE genes for skincolor 3 levle in whole blood sample
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA

example_TMM_DE <- readRDS("~/ses-1/user_wx/genowide_with1k_sc_allsample.rds")
data_DE = example_TMM_DE %>% 
  hoist(out, ttT = list("result", "ttT")) %>% 
  filter(map_lgl(.$ttT, ~dim(.)[1]!=0)) %>% 
  mutate(gene_sig = map(ttT, ~ filter(., adj.P.Val<0.05) %>% pull(gene)) %>% set_names(treatment)) %>% 
  # mutate(min_p = map(ttT, ~ dplyr::slice_min(.,adj.P.Val) %>% pull(adj.P.Val)) %>% set_names(treatment)) %>% 
  dplyr::select(-out, -gene_set_name) %>% 
  filter(controls == "basic")

venn::venn(data_DE$gene_sig, zcolor = "style",
           ellipse = FALSE, opacity = 0.15, ilcs = 1,
           box = FALSE,
           sncs = 1)

UpSetR::upset(UpSetR::fromList(data_DE$gene_sig), 
              matrix.color ="darkred", 
              main.bar.color = "gray23",
              # sets.bar.color=c("navajowhite4","lightsteelblue1","orange","goldenrod3","cadetblue1"),
              sets.bar.color=c("navajowhite4","lightsteelblue1"),
              order.by = "degree")


#' ###  bubble polot of fig1 panelA for each SES indicator
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA



res_fun = function(example0_with1k, example0_without1k) {
  ex0_without1k <- example0_without1k %>%
    filter(gene_set_name!="inflam1k_mRNA" ) %>% 
    hoist(out, m = list("result", "m12_fdr", 1, "other", "m")) %>% 
    filter(!map_lgl(m, ~is.null(.x))) %>% 
    mutate(p = m %>% map(~ dplyr::slice(.x, which.min(adj.p.within))) %>% map_dbl(~ .x %>% pluck("adj.p.within"))) %>% 
    dplyr::select(treatment, gene_set_name, controls, p) %>% 
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
                           p>0.05 ~ 0.0000001)
    )%>%
    mutate(
      # gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
      #      gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
      #      gene_set_name = gene_set_name %>% replace(gene_set_name == "Cvd", "CVD"),
      #      gene_set_name = gene_set_name %>% replace(gene_set_name == "Copd", "COPD"),
      #      gene_set_name = gene_set_name %>% replace(gene_set_name == "Ckd", "CKD"),
      #      gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI"),
      "1KI Genes" = "Without 1KI Genes"%>% as.factor())
  
  
  ex0_with1k <- example0_with1k %>%
    hoist(out, m = list("result", "m12_fdr", 1, "other", "m")) %>% 
    filter(!map_lgl(m, ~is.null(.x))) %>% 
    mutate(p = m %>% map(~ dplyr::slice(.x, which.min(adj.p.within))) %>% map_dbl(~ .x %>% pluck("adj.p.within"))) %>% 
    dplyr::select(treatment, gene_set_name, controls, p) %>% 
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
                           p>0.05 ~ 0.0000001)
    )%>%
    mutate(
      # gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
      #      gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
      #      gene_set_name = gene_set_name %>% replace(gene_set_name == "Cvd", "CVD"),
      #      gene_set_name = gene_set_name %>% replace(gene_set_name == "Copd", "COPD"),
      #      gene_set_name = gene_set_name %>% replace(gene_set_name == "Ckd", "CKD"),
      #      gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI"),
      "1KI Genes" = "With 1KI Genes"%>% as.factor())
  
  bind_rows(ex0_with1k, ex0_without1k)
}

# all sample

example0_without1k <- readRDS("~/ses-1/user_wx/m12_without1k_sc3levels_allsample.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_sc3levels_allsample.rds")

example0_without1k <- readRDS("~/ses-1/user_wx/m12_without1k_sccont_allsample.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_sccont_allsample.rds")

example0_without1k <- readRDS("~/ses-1/user_wx/m12_without1k_sc5levles_allsample.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_sc5levels_allsample.rds")

# nonhispanic black sample
example0_without1k <- readRDS("~/ses-1/user_wx/m12_without1k_scbinary_NonBstrata.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_scbinary_NonBstrata.rds")


res1 = res_fun(example0_with1k, example0_without1k)

example0_without1k <- readRDS("~/ses-1/user_wx/m12_without1k_sc4levels_NonBstrata.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_sc4levels_NonBstrata.rds")

example0_without1k <- readRDS("~/ses-1/user_wx/m12_without1k_aging_sc4levels_NonBstrata_correctioninunion.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_aging_sc4levels_NonBstrata_correctioninunion.rds")

example0_without1k <- readRDS("~/ses-1/user_wx/m12_without1k_aging_sc4levels_NonBstrata_correctioninunion.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_aging_sc4levels_NonBstrata_correctioninunion.rds")


res2 = res_fun(example0_with1k, example0_without1k)

example0_without1k <- readRDS("~/ses-1/user_wx/m12_without1k_sccont_NonBstrata.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_sccont_NonBstrata.rds")

res3 = res_fun(example0_with1k, example0_without1k)

res = bind_rows(res1,res2, res3) %>% 
  # filter(gene_set_name != "inflame_mRNA") %>% 
  dplyr::select(-pval, - pval2) %>%
  filter(controls=="all", gene_set_name %in% c(aging,table1)) %>%
  dplyr::select("treatment", "gene_set_name", "p", "1KI Genes")


# nonhispanic black sample with PC1- PC4
example0_without1k <- readRDS("~/ses-1/user_wx/m12_without1k_scbinary_NonBstrata_PC4.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_scbinary_NonBstrata_PC4.rds")

res1 = res_fun(example0_with1k, example0_without1k)

example0_without1k <- readRDS("~/ses-1/user_wx/m12_without1k_sc4levels_NonBstrata_PC4.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_sc4levels_NonBstrata_PC4.rds")

res2 = res_fun(example0_with1k, example0_without1k)

example0_without1k <- readRDS("~/ses-1/user_wx/m12_without1k_sccont_NonBstrata_PC4.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_sccont_NonBstrata_PC4.rds")

res3 = res_fun(example0_with1k, example0_without1k)

res = bind_rows(res1,res2, res3) %>% 
  # filter(gene_set_name != "inflame_mRNA") %>% 
  dplyr::select(-pval, - pval2) %>%
  filter(controls=="all", gene_set_name %in% c(aging,table1)) %>%
  dplyr::select("treatment", "gene_set_name", "p", "1KI Genes")
res %>% kableExtra::kable() %>% kableExtra::kable_styling()

res$pval2 = -log10(res$p)
res$pval2[which(res$pval2>-log10(0.0001))] = -log10(0.0001)
res$treatment = as.character(res$treatment)
res$gene_set_name = as.character(res$gene_set_name)


colnames(res) = c(colnames(res)[1:3], "Cluster", colnames(res)[5])
res = res[order(res$Cluster),]

res$Instance = 0
for (i in 1:length(res$treatment)) {
  res$Instance[i] = length(which(res$treatment==res$treatment[i] & res$gene_set_name==res$gene_set_name[i]))
}



ggplot(res[res$Instance==1,], aes(x = treatment,y  = gene_set_name, color = Cluster, fill = Cluster,  size = pval2)) +
  geom_point(alpha  = 0.5,shape = 21, stroke = 1) +
  
  geom_point(data = res[res$Instance==2 & !duplicated(res[,c(1,2)]),], inherit.aes = T, position = position_nudge(y = +0.15), alpha = 0.5, stroke = 1,show.legend = F) + 
  geom_point(data = res[res$Instance==2 & duplicated(res[,c(1,2)]),], inherit.aes = T, position = position_nudge(y = -0.15), alpha = 0.5, stroke = 1,show.legend = F) +
  
  theme_bw(base_size = 12) +
  scale_color_manual(values = c("#01665e","#ff7f0e"), name = "1KI Genes") +
  scale_fill_manual(values = c("#01665e","#ff7f0e"), name = "1KI Genes") +
  scale_size_continuous(range = c(1,12),
                        name = "Adjusted\n p-value", 
                        limits = c(-log10(0.05), -log10(0.0001)),
                        breaks = c(-log10(0.05),-log10(0.01),-log10(0.001), -log10(0.0001)), 
                        labels = c("p<0.05", "p<0.01","p<0.001","p<0.0001")) + 
  xlab("SES Indicators") + ylab("mRNA Signatures") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title = element_text(size = 11, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=10, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=10,face = "bold" ,family = "Calibri", angle = 30,hjust = 1)) +
  theme(legend.text=element_text(size=10, family = "Calibri")) +
  theme(legend.title = element_text(size = 10, face = "bold", family = "Calibri")) +
  theme(legend.position="right") +
  scale_y_discrete(limits = rev)  +
  guides(color = guide_legend(override.aes = list(size = 7))) +
  guides(fill = guide_legend(override.aes = list(size = 7))) +
  guides(shape = guide_legend(override.aes = list(size = 7))) 

