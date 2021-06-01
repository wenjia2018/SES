
#' ---
#' title: figure 1
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

#' #### Figure1 panel B PCA
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
source("/home/xu/ses-1/user_wx/extract_v2.R")
# example_race = readRDS("~/ses-1/user_wx/race_bespoke_12.02.2021.rds")
# example_skincolor3 = readRDS("~/ses-1/user_wx/color3_bespoke_18.02.2021.rds")
example_race = readRDS("~/ses-1/user_wx/race_bespoke_15.03.2021.rds")
example_skincolor3 = readRDS("~/ses-1/user_wx/color3_bespoke_15.03.2021.rds")

p_eqtl = c(0.01, 0.05)

fig1A = function(p_eqtl){
  race = outm7pca(p = p_eqtl, control = "ancestryPC_ses", example_race) 
  skincolor = outm7pca(p = p_eqtl, control = "ancestryPC_ses", example_skincolor3)
  ex0_m7pca = race %>% 
    rbind(skincolor) %>%
    filter(gene_set_name %>% str_detect("aging")) %>% 
    filter(gene_set_name!="aging_down_cl1_mRNA") %>% 
    dplyr::select(treatment, gene_set_name, p_pca, p_id) %>% 
    mutate(p = case_when(gene_set_name == "aging_up_cl2_mRNA" ~ p_pca*9,
                         gene_set_name == "aging_up_cl3_mRNA" ~ p_pca*7,
                         gene_set_name == "aging_up_cl4_mRNA" ~ p_pca*5,
                         TRUE ~ p_pca*10)) %>% 
    filter(p<0.05) %>% 
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
                                                            "Aging Cluster Complement",
                                                            "Aging Up",
                                                            "Aging Down",
                                                            "Aging Composite"
           )) %>% fct_rev
           # gene_set_name = str_c(gene_set_name, "_", p_id)
    )
  
  ex0_m7pca =
    ex0_m7pca %>%
    group_by(treatment, gene_set_name) %>% 
    mutate(noPC = n()) %>% 
    mutate(noPC = ifelse(noPC>1, noPC, NA)) %>% 
    mutate(group = ifelse(treatment %>% str_detect("Hispanic"), "Race", "Skin Color"))
  # ex0_m8fdr$gene_set_name <- factor(ex0_m8fdr$gene_set_name,
  #                                   levels = c("1KI", "CTRA", "Antbintf",  "Interferon", "Inflamation in CTRA",                                                                    "Aging", "Aging Down", "Aging Up",
  #                                              "Aging Down Cl1",  "Aging Down Cl1a", "Aging Down Cl1b", "Aging Down Cl1c",
  #                                              "Aging Down Cl2",  "Aging Down Cl3",
  #                                              "Aging Up Cl1",    "Aging Up Cl2",    "Aging Up Cl3",    "Aging Up Cl4"  )) %>% fct_rev
  axiscolor = c("darkblue","darkblue","darkblue", "darkblue", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30","grey30", "grey30", "grey30")
  ggplot(ex0_m7pca, aes(treatment, gene_set_name, size = pval2,
                        colour = treatment,
                        # fill = treatment
  )) +
    # geom_point(
    #   # stroke = 1.5, 
    #            shape = 21,
    #            alpha = 0.4, 
    #            # colour = "darkblue",
    #            fill = "navy"
    #            # position = position_jitterdodge()
    #            ) +
    geom_jitter(position = position_jitter(width = 0.5, height = 0.15, seed = 1)) +
    geom_text(aes(label = p_id), 
              show.legend=FALSE, 
              size=2,
              color ="black",
              position = position_jitter(width = 0.5, height = 0.15, seed = 1)) +
      facet_wrap(~group, scales = "free_x") +
    # geom_text(aes(label=noPC), show.legend=FALSE, size=3, color ="white")+
    # scale_fill_manual(values = c("darkblue", "goldenrod3")) +
    # scale_color_manual(values = c("darkblue", "goldenrod3")) +
    # geom_jitter(height = 0.00000025) +
    # gghighlight(class == "inflam") +
    theme_bw() +
    labs(
      caption = paste("p_eqtl = ", p_eqtl),
      title = "Figure . PCA regression",
      colour = "Race Skin Color",
      y = "mRNA Signatures PC",
      x = "Race Skincolor") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          # plot.margin=unit(c(1, 1, 0.1, 1), "cm"),
          axis.text.y = element_text(color=axiscolor),
          text = element_text(size=10, face = "bold")) +
    scale_color_manual(values = c("goldenrod1", "goldenrod3", "lightblue","cornflowerblue")) +
    # scale_color_brewer(palette = "Dark2")
    # scale_color_viridis(option = "D", discrete = TRUE)+
    scale_size_continuous(name = "Adjusted P-value", range = c(0, 21), 
                          limits = c(0.0000001, 100000), breaks = c(
                            # 0.0000001,
                            10000, 15000, 25000, 100000),
                          labels = c(
                            # "n.s.", 
                            "p<0.05", "p<0.01", "p<0.001", "p<0.0001")) +
    scale_alpha(guide = 'none') +
    guides(shape = guide_legend(override.aes = list(size = 10)),
           fill = guide_legend(override.aes = list(size = 8)),
           # colour = FALSE,
           size = guide_legend(override.aes = list(colour = "lightgrey")))    
}

p_eqtl %>% map(fig1A)