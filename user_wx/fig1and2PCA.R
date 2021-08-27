
#' ---
#' title: figure 1 and 2 PCA results
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

#' #### Figure1 PCA 
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10
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
walk(dir(path = here("R"), full.names = TRUE), source)
source("/home/xu/ses-1/user_wx/extract_v2.R")
pca_corrected = function(data, p_eqtl, adjmethod) {
  skincolor = outm7pca(p = p_eqtl, control = "ancestryPC_ses", data)
  ex0_m7pca = skincolor %>% 
    filter(gene_set_name %>% str_detect("aging")) %>% 
    filter(gene_set_name!="aging_down_cl1_mRNA") %>% 
    dplyr::select(treatment, gene_set_name, p_pca, p_id) %>% 
    filter(p_pca!="NA")
  if(adjmethod=="fdr") {
    ex0_m7pca = ex0_m7pca %>% 
      group_by(treatment) %>% 
      # using fdr for correction
      mutate(p_adj= p_pca %>% p.adjust("fdr")) %>% 
      ungroup
  } else if(adjmethod =="bonferroni") {
    ex0_m7pca = ex0_m7pca %>% 
      mutate(p_adj = case_when(gene_set_name == "aging_up_cl2_mRNA" ~ p_pca*9*13,
                           gene_set_name == "aging_up_cl3_mRNA" ~ p_pca*7*13,
                           gene_set_name == "aging_up_cl4_mRNA" ~ p_pca*5*13,
                           TRUE ~ p_pca*10*13))
  }
  ex0_m7pca = ex0_m7pca %>%
    # filter(p<0.05) %>% 
    mutate(
      treatment= case_when(treatment == "raceethnicity_Hispanic" ~ "Hispanic",
                           treatment =="raceethnicity_NonHblack" ~  "Non Hispanic Black",
                           treatment =="color_byinterviewer3_LightMed"  ~ "Light Medium",  
                           treatment =="color_byinterviewer3_DarkBlack"  ~ "Dark Black")) %>%
    mutate(treatment = factor(treatment, levels = c("Hispanic", "Non Hispanic Black","Light Medium","Dark Black"))) %>%
    mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
           gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
           # gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI"),
           # gene_set_name = gene_set_name %>% replace(gene_set_name == "Ctra", "CTRA"),
           # gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflame", "Inflamation in CTRA"),
           gene_set_name= case_when(gene_set_name =="Aging Up Cl1" ~ "Innate/Adaptive Immunity",
                                    gene_set_name =="Aging Up Cl2" ~ "Actin Regulation",#/Focal Adhesion/Tight Junctions",
                                    gene_set_name =="Aging Up Cl3" ~ "Fatty Acid Metabolism",#/Peroxisome Activity",
                                    gene_set_name =="Aging Up Cl4" ~ "Lysosome Metabolism",#/Glycosaminoglycan Degradation",
                                    # gene_set_name =="Aging Down Cl1" ~ "Aging Down Cl1",
                                    gene_set_name =="Aging Down Cl1a" ~ "RNA Metabolism",#/Ribosome Biogenesis/Purine Metabolism",
                                    gene_set_name =="Aging Down Cl1b" ~ "Mitochondrial",
                                    gene_set_name =="Aging Down Cl1c" ~ "DNA Replication",#/Elongation/Mismatch Repair",
                                    gene_set_name =="Aging Down Cl2" ~ "Ribosome",
                                    gene_set_name =="Aging Down Cl3" ~ "Immune Genes",
                                    gene_set_name =="Aging" ~ "Aging Composite",
                                    gene_set_name == "Inflam1k" ~ "1KI",
                                    gene_set_name == "Ctra" ~ "CTRA",
                                    gene_set_name == "Inflame" ~ "Inflamation in CTRA",
                                    TRUE ~ gene_set_name),
           gene_set_name = factor(gene_set_name, levels = c("1KI", "CTRA", "Antbintf",  "Interferon", "Inflamation in CTRA",
                                                            
                                                            "Aging Down Cl1",
                                                            "Lysosome Metabolism",#/Glycosaminoglycan Degradation"
                                                            "Fatty Acid Metabolism",#/Peroxisome Activity", 
                                                            "Actin Regulation",#/Focal Adhesion/Tight Junctions",
                                                            "Ribosome",
                                                            "DNA Replication",#/Elongation/Mismatch Repair",
                                                            "Mitochondrial",
                                                            "RNA Metabolism",#/Ribosome Biogenesis/Purine Metabolism",
                                                            "Immune Genes",
                                                            "Innate/Adaptive Immunity",
                                                            "Aging Up",
                                                            "Aging Down",
                                                            "Aging Cluster Complement",
                                                            "Aging Composite"
           )) %>% fct_rev
           # gene_set_name = str_c(gene_set_name, "_", p_id)
    )
  
  ex0_m7pca =
    ex0_m7pca %>%
    mutate(gene_set_name = fct_recode(gene_set_name, "Lysosome Metabolism (n=5)"="Lysosome Metabolism",
                                      "Innate/Adaptive Immunity (n=68)" = "Innate/Adaptive Immunity",
                                      "Fatty Acid Metabolism (n=7)"="Fatty Acid Metabolism",
                                      "Actin Regulation (n=9)"="Actin Regulation",
                                      "Immune Genes (n=44)"="Immune Genes",
                                      "Ribosome (n=12)"="Ribosome",
                                      "Mitochondrial (n=29)"="Mitochondrial",
                                      "RNA Metabolism (n=40)"="RNA Metabolism",
                                      "DNA Replication (n=16)"="DNA Replication",
                                      "Aging Up (n=438)"="Aging Up",
                                      "Aging Down (n=610)"="Aging Down",
                                      "Aging Cluster Complement (n=818)"="Aging Cluster Complement",
                                      "Aging Composite (n=1048)"="Aging Composite"))
    
}
p_eqtl <- c(0.05)
example_skincolor3 = readRDS("~/ses-1/user_wx/skincolor_eqtl005_aging_composite_ancestry_11.05.2021.rds")
pca_wholeblood = pca_corrected(example_skincolor3, p_eqtl, adjmethod = "fdr")
pca_wholeblood %>% 
  filter(p_adj<0.05) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()

#' #### Figure2 PCA 
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

example_skincolor3_black = readRDS("~/ses-1/user_wx/skincolor3_NonHblack_strata_eqtl005_aging_composite_ancestry_12.05.2021.rds")


pca_black = pca_corrected(example_skincolor3_black, p_eqtl, adjmethod = "fdr")
pca_black %>% 
  mutate(treatment = "Dark Black") %>% 
  filter(p_adj<0.05) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()


