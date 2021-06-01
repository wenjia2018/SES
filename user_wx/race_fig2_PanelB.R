
#' ---
#' title: figure 2
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

#' #### Figure2 panel B PCA
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
walk(dir(path = here("R"), full.names = TRUE), source)
source("/home/xu/ses-1/user_wx/extract_v2.R")
black <- readRDS("/home/xu/ses-1/user_wx/color3_bespoke_NonHblack_strata_16.03.2021.rds")
white <- readRDS("/home/xu/ses-1/user_wx/color3_bespoke_NonHwhite_strata_15.03.2021.rds")
hispanic <- readRDS("/home/xu/ses-1/user_wx/color3_bespoke_Hispanic_strata_15.03.2021.rds")
p_eqtl <- c(0.05, 0.01)

fig1A <- function(p_eqtl) {
  data_black <- outm7pca(p = p_eqtl, control = "ancestryPC_ses", black) %>% mutate(group = "Non Hispanic Black")
  data_white <- outm7pca(p = p_eqtl, control = "ancestryPC_ses", white) %>% mutate(group = "Non Hispanic White")
  data_hispanic <- outm7pca(p = p_eqtl, control = "ancestryPC_ses", hispanic) %>% mutate(group = "Hispanic")

  ex0_m7pca <- data_black %>%
    rbind(data_white) %>%
    rbind(data_hispanic) %>%
    filter(gene_set_name %>% str_detect("aging")) %>%
    filter(gene_set_name != "aging_down_cl1_mRNA") %>%
    dplyr::select(treatment, gene_set_name, p_pca, p_id, group) %>%
    mutate(p = case_when(
      gene_set_name == "aging_up_cl2_mRNA" ~ p_pca * 9,
      gene_set_name == "aging_up_cl3_mRNA" ~ p_pca * 7,
      gene_set_name == "aging_up_cl4_mRNA" ~ p_pca * 5,
      TRUE ~ p_pca * 10
    )) %>%
    filter(p < 0.05) %>%
    mutate(
      pval = case_when(
        p < 0.0001 ~ 0.0001,
        p < 0.001 ~ 0.001,
        p < 0.01 ~ 0.01,
        p < 0.05 ~ 0.05,
        p > 0.05 ~ 100
      ),
      pval2 = case_when(
        p < 0.0001 ~ 100000,
        p < 0.001 ~ 25000,
        p < 0.01 ~ 15000,
        p < 0.05 ~ 10000,
        p > 0.05 ~ 0.0000001
      ),
      treatment = case_when(
        treatment == "raceethnicity_Hispanic" ~ "Hispanic",
        treatment == "raceethnicity_NonHblack" ~ "Non Hispanic Black",
        treatment == "color_byinterviewer3_LightMed" ~ "Light Medium",
        treatment == "color_byinterviewer3_DarkBlack" ~ "Dark Black"
      )
    ) %>%
    mutate(treatment = factor(treatment, levels = c("Hispanic", "Non Hispanic Black", "Light Medium", "Dark Black"))) %>%
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
        "Aging Cluster Complement",
        "Aging Up",
        "Aging Down",
        "Aging Composite"
      )) %>% fct_rev()
      # gene_set_name = str_c(gene_set_name, "_", p_id)
    )

  example <-
    ex0_m7pca %>%
    group_by(treatment, gene_set_name, group) %>%
    mutate(noPC = n()) %>%
    mutate(
      noPC = ifelse(noPC > 1, noPC, NA),
      treatment = case_when(
        group == "Non Hispanic Black" ~ "Light Medium v.s Dark Black",
        TRUE ~ "Light Medium v.s White"
      )
    )
  # ex0_m8fdr$gene_set_name <- factor(ex0_m8fdr$gene_set_name,
  #                                   levels = c("1KI", "CTRA", "Antbintf",  "Interferon", "Inflamation in CTRA",                                                                    "Aging", "Aging Down", "Aging Up",
  #                                              "Aging Down Cl1",  "Aging Down Cl1a", "Aging Down Cl1b", "Aging Down Cl1c",
  #                                              "Aging Down Cl2",  "Aging Down Cl3",
  #                                              "Aging Up Cl1",    "Aging Up Cl2",    "Aging Up Cl3",    "Aging Up Cl4"  )) %>% fct_rev
  axiscolor <- c("darkblue", "darkblue", "darkblue", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30")
  ggplot(example, aes(treatment, gene_set_name,
    size = pval2
    # fill = `1KI Genes`,
    # colour = `1KI Genes`
  )) +
    geom_point(stroke = 1.5, shape = 21, alpha = 0.4, colour = "darkblue", fill = "navy") +
    geom_text(aes(label = noPC), show.legend = FALSE, size = 3, color = "white") +
    facet_wrap(~group, scales = "free_x") +
    # scale_fill_manual(values = c("darkblue", "goldenrod3")) +
    # scale_color_manual(values = c("darkblue", "goldenrod3")) +
    # geom_jitter(height = 0.00000025) +
    # gghighlight(class == "inflam") +
    theme_bw() +
    labs(
      caption = paste("p_eqtl = ", p_eqtl),
      title = "Figure . PCA regression",
      y = "mRNA Signatures PC",
      x = "Skincolor") +
    theme(
      axis.text.x = element_text(angle = 0),
      # plot.margin=unit(c(1, 1, 0.1, 1), "cm"),
      axis.text.y = element_text(color = axiscolor),
      text = element_text(size = 10, face = "bold")
    ) +
    scale_size_continuous(
      name = "Adjusted P-value", range = c(0, 21),
      limits = c(0.0000001, 100000), breaks = c(0.0000001, 10000, 15000, 25000, 100000),
      labels = c("n.s.", "p<0.05", "p<0.01", "p<0.001", "p<0.0001")
    ) +
    scale_alpha(guide = "none") +
    guides(
      shape = guide_legend(override.aes = list(size = 10)),
      fill = guide_legend(override.aes = list(size = 8))
    )
}
p_eqtl %>% map(fig1A)
