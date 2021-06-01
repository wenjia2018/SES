#' ---
#' title: aging up PCA loading per cluster plotting
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

#' #### Figure aging up PCA loadings
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

walk(dir(path = here("R"),full.names = TRUE), source)
source("/home/xu/ses-1/user_wx/extract_v2.R")
load_data(reconciled = FALSE, remove_inflam = FALSE)
example_race = readRDS("~/ses-1/user_wx/race_bespoke_12.02.2021.rds")
example_skincolor3 = readRDS("~/ses-1/user_wx/color3_bespoke_18.02.2021.rds")
threshold = 0.005
p_eqtl = c(0.05, 0.01)
race = outm7pca(p = 0.05, control = "ancestryPC_ses", example_race) 
skincolor = outm7pca(p = 0.05, control = "ancestryPC_ses", example_skincolor3)

bespoke = example_skincolor3 %>% rbind(example_race)
outm71 = p_eqtl %>% map(outm7pca, control, bespoke)


outm71 %>%
  bind_rows() %>%
  dplyr::filter(p_pca < threshold) %>%
  select(-coef) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_pca) %>%
  filter(gene_set_name=="aging_up_mRNA") %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_pca=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  select(-well_loaded, -loadings) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()


fig1A = function(dim){
  ex0_m7pca = race %>% 
    rbind(skincolor) %>%
    filter(gene_set_name %>% str_detect("aging")) %>% 
    filter(gene_set_name!="aging_down_cl1_mRNA") %>% 
    dplyr::select(treatment, gene_set_name, p_pca, p_id, var_explained, loadings) %>% 
    mutate(p = case_when(gene_set_name == "aging_up_cl2_mRNA" ~ p_pca*9,
                         gene_set_name == "aging_up_cl3_mRNA" ~ p_pca*7,
                         gene_set_name == "aging_up_cl4_mRNA" ~ p_pca*5,
                         TRUE ~ p_pca*10)) %>% 
    # filter(p<0.05) %>% 
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
                                    gene_set_name =="Aging Up Cl2" ~ "Actin Cytoskeleton",#/Focal Adhesion/Tight Junctions",
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
                                                            "RNA Metabolism",#/Ribosome Biogenesis/Purine Metabolism",
                                                            "Mitochondrial",
                                                            "DNA Replication",#/Elongation/Mismatch Repair",
                                                            "Ribosome",
                                                            "Immune Genes",
                                                            "Innate/Adaptive Immunity", 
                                                            "Actin Cytoskeleton",#/Focal Adhesion/Tight Junctions",
                                                            "Fatty Acid Metabolism",#/Peroxisome Activity", 
                                                            "Lysosome Metabolism",#/Glycosaminoglycan Degradation"
                                                            "Aging Up",
                                                            "Aging Down",
                                                            "Aging Composite"
           )) %>% fct_rev
           # gene_set_name = str_c(gene_set_name, "_", p_id)
    )
  a <-
    ex0_m7pca %>%
    filter(gene_set_name == "Aging Up" & treatment == "Light Medium") %>%
    filter(p_id == dim) %>%
    pull(loadings) %>%
    as.data.frame() %>%
    rownames_to_column("gene")
  
  cluster_union = reduce(names(signatures$outcome_set) %>% 
                           str_detect('cl') %>%
                           keep(signatures$outcome_set, .), union)
  
  signatures$outcome_set$aging_cluster_complement =  setdiff(signatures$outcome_set$aging_mRNA, cluster_union)
  b = names(signatures$outcome_set) %>% 
    str_detect('aging') %>%
    keep(signatures$outcome_set, .)
  
  col_names = b %>% names
  for(i in 1: length(col_names)){
    var = col_names[i]
    a = 
      a %>% mutate(!!var := ifelse(gene %in% b[[!!var]], !!sym(dim), NA))
  }
  
  a %>% 
    select(-1, -2) %>%
    pivot_longer(everything(),
                 names_to = "gene_set_name",
                 values_to = dim
    ) %>%
    filter(gene_set_name != "aging_down_cl1_mRNA") %>%
    mutate(
      gene_set_name = gene_set_name %>% str_replace_all("_mRNA", ""),
      gene_set_name = gene_set_name %>% str_replace_all("_", " ") %>% str_to_title(),
      # gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI"),
      # gene_set_name = gene_set_name %>% replace(gene_set_name == "Ctra", "CTRA"),
      # gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflame", "Inflamation in CTRA"),
      gene_set_name = case_when(
        gene_set_name == "Aging Up Cl1" ~ "Innate/Adaptive Immunity",
        gene_set_name == "Aging Up Cl2" ~ "Actin Cytoskeleton", # /Focal Adhesion/Tight Junctions",
        gene_set_name == "Aging Up Cl3" ~ "Fatty Acid Metabolism", # /Peroxisome Activity",
        gene_set_name == "Aging Up Cl4" ~ "Lysosome Metabolism", # /Glycosaminoglycan Degradation",
        # gene_set_name =="Aging Down Cl1" ~ "Aging Down Cl1",
        gene_set_name == "Aging Down Cl1a" ~ "RNA Metabolism", # /Ribosome Biogenesis/Purine Metabolism",
        gene_set_name == "Aging Down Cl1b" ~ "Mitochondrial",
        gene_set_name == "Aging Down Cl1c" ~ "DNA Replication", # /Elongation/Mismatch Repair",
        gene_set_name == "Aging Down Cl2" ~ "Ribosome",
        gene_set_name == "Aging Down Cl3" ~ "Immune Genes",
        gene_set_name == "Aging" ~ "Aging Composite",
        gene_set_name == "Inflam1k" ~ "1KI",
        gene_set_name == "Ctra" ~ "CTRA",
        gene_set_name == "Inflame" ~ "Inflamation in CTRA",
        TRUE ~ gene_set_name
      ),
      gene_set_name = factor(gene_set_name, levels = c(
        "1KI", "CTRA", "Antbintf", "Interferon", "Inflamation in CTRA",
        
        "Aging Down Cl1",
        "RNA Metabolism", # /Ribosome Biogenesis/Purine Metabolism",
        "Mitochondrial",
        "DNA Replication", # /Elongation/Mismatch Repair",
        "Ribosome",
        "Immune Genes",
        "Innate/Adaptive Immunity",
        "Actin Cytoskeleton", # /Focal Adhesion/Tight Junctions",
        "Fatty Acid Metabolism", # /Peroxisome Activity",
        "Lysosome Metabolism", # /Glycosaminoglycan Degradation"
        "Aging Up",
        "Aging Down",
        "Aging Cluster Complement",
        "Aging Composite"
      )) 
      # gene_set_name = str_c(gene_set_name, "_", p_id)
    ) %>% 
    ggplot(aes(gene_set_name, !!sym(dim))) +
    geom_violin()+
    ggplot2::stat_summary(
      fun.data = mean_sdl,
      fun.args = list(mult = 1),
      geom = "pointrange") +
    ggplot2::theme(legend.position = "none") +
    theme(
      axis.text.x = element_text(angle = 10, hjust = 1),
      # plot.margin=unit(c(1, 1, 0.1, 1), "cm"),
      text = element_text(size = 10, face = "bold")) +
    labs(
      # caption = paste("p_eqtl = ", p_eqtl),
      title = paste("Figure . PCA", dim, "loading from Aging Up for aging subsets and clusters"),
      y = paste("PCA", dim, "loadings"),
      x = "Aging Gene Sets"
    )  
}
no.dim = str_c("d", 1:10)

map(no.dim, fig1A)