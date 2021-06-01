
#' ---
#' title: numbers of PCA well loaded genes and intersection between dimensions
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE


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

control = "ancestryPC_ses"
example_skincolor3 = readRDS("~/ses-1/user_wx/skincolor_eqtl005_aging_composite_ancestry_11.05.2021.rds")
p_eqtl = 0.05
# outm71 = p_eqtl %>% map(outm7pca, control, example_skincolor3)

outm71 = f0(p_eqtl, example_skincolor3)
temp = 
  outm71 %>%
  hoist(out, loadings = list("result", "m7_ob", 1, "other", "loadings")) %>%
  filter(loadings!="NULL") %>% 
  mutate(loadings = loadings %>% map(~.[]) %>% map(~ as.data.frame(.)),
         loadings = loadings %>% map(~ set_names(., str_c("d", 1:length(.)))),
         loadings = loadings %>% map(~ split.default(., seq_along(.)))) %>% 
  hoist(out, well_loaded = list("result", "m7_ob", 1, "other", "well_loaded")) %>% 
  mutate(well_loaded = well_loaded %>% map(~ set_names(.x, str_c("d", 1:length(.x))))) %>% 
  mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
         gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
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
         )) %>% fct_rev) %>% 
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

varying_threshold = function(data, threshold){
  temp2 = 
    data %>%
    hoist(out, loadings = list("result", "m7_ob", 1, "other", "loadings")) %>%
    filter(loadings!="NULL") %>% 
    mutate(loadings = loadings %>% map(~.[]) %>% map(~ as.data.frame(.)),
           loadings = loadings %>% map(~ set_names(., str_c("d", 1:length(.)))),
           loadings = loadings %>% map(~ split.default(., seq_along(.)))) %>% 
    select(treatment, gene_set_name, loadings) %>% 
    mutate(loadings = loadings %>% map_depth(2, ~.x %>% rownames_to_column("gene"))) %>% 
    mutate(loadings = loadings %>% map_depth(2, ~ .x %>% filter(across(where(is.numeric), ~ abs(.)> threshold)))) %>% 
    mutate(well_loaded = loadings %>% map_depth(2, ~ .x %>% pull(gene))) %>% 
    mutate(well_loaded = well_loaded %>% map(~ set_names(.x, str_c("d", 1:length(.x))))) %>% 
    mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
           gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
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
           )) %>% fct_rev) %>% 
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
                                      "Aging Composite (n=1048)"="Aging Composite")) %>% 
    filter(treatment %>% str_detect("DarkBlack")) 
  temp3 = 
    temp2 %>% 
    pull(well_loaded) %>%
    set_names(temp2$gene_set_name)
  
  temp3 = temp3[c("Aging Composite (n=1048)",
                  "Aging Cluster Complement (n=818)",
                  "Aging Down (n=610)",
                  "Aging Up (n=438)",
                  "Innate/Adaptive Immunity (n=68)",
                  "Immune Genes (n=44)",
                  "RNA Metabolism (n=40)",
                  "Mitochondrial (n=29)",
                  "DNA Replication (n=16)",
                  "Ribosome (n=12)",
                  "Actin Regulation (n=9)",
                  "Fatty Acid Metabolism (n=7)",
                  "Lysosome Metabolism (n=5)")]
  
  temp3 %>% 
    map(~ crossprod(table(stack(.))))
}

# https://stackoverflow.com/questions/24614391/intersect-all-possible-combinations-of-list-elements

#' ### threshold 0.1
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10
varying_threshold(outm71, 0.1)


#' ### threshold 0.2
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10
varying_threshold(outm71, 0.2)


#' ### threshold 0.3
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10
varying_threshold(outm71, 0.3)


#' ### threshold 0.4
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10
varying_threshold(outm71, 0.4)


#' ### threshold 0.5
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10
varying_threshold(outm71, 0.5)
