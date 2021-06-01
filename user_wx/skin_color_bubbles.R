
#' #### skin color different measures bubble plots
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
example = readRDS("/home/xu/ses-1/user_wx/color_conti_binary_dummy3and5_bespoke.rds")


p_eqtl <- c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)

fig1A = function(p_eqtl){
  data = m8_present(p = p_eqtl, control = "ancestryPC_ses", example) %>% 
    filter(!is.na(p_omnibus))

  ex0_m8fdr = data %>% 
    dplyr::select(treatment, gene_set_name, p_omnibus) %>% 
    rename(p = p_omnibus) %>% 
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
           treatment= case_when(treatment == "color_byinterviewer_continuous" ~ "Continuous",
                                treatment == "color_byinterviewer_binary" ~  "Binary",
                                treatment == "color_byinterviewer3_LightMed"  ~ "Light Medium",  
                                treatment == "color_byinterviewer3_DarkBlack"  ~ "Dark Black",
                                treatment == "color_byinterviewer5_Light"  ~ "Light",  
                                treatment == "color_byinterviewer5_Medium"  ~ "Medium",
                                treatment == "color_byinterviewer5_Dark"  ~ "Dark",  
                                treatment == "color_byinterviewer5_Black"  ~ "Black")
    )  %>%
    mutate(treatment = factor(treatment, levels = c("Continuous", "Binary","Light Medium","Dark Black", "Light", "Medium", "Dark", "Black"
    ))) %>%
    mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
           gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Ctra", "CTRA"),
           gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflame", "Inflamation in CTRA"),)
  
  # need to specify shape, because some shape does not have borders. 
  # reorder y axis
  ex0_m8fdr$gene_set_name <- factor(ex0_m8fdr$gene_set_name, levels = c("1KI", "CTRA", "Antbintf",  "Interferon", "Inflamation in CTRA",
                                                                        "Aging", "Aging Down", "Aging Up",
                                                                        "Aging Down Cl1",  "Aging Down Cl1a", "Aging Down Cl1b", "Aging Down Cl1c",
                                                                        "Aging Down Cl2",  "Aging Down Cl3",
                                                                        "Aging Up Cl1",    "Aging Up Cl2",    "Aging Up Cl3",    "Aging Up Cl4"
                                                                        
  )) %>% fct_rev
  
  
  ggplot(ex0_m8fdr, aes(treatment, gene_set_name, size = pval2
                        # fill = `1KI Genes`,
                        # colour = `1KI Genes`
  )) +
    geom_point(stroke = 1.5, shape = 21, alpha = 0.4, colour = "darkblue", fill = "navy") +
    # scale_fill_manual(values = c("darkblue", "goldenrod3")) +
    # scale_color_manual(values = c("darkblue", "goldenrod3")) +
    # geom_jitter(height = 0.00000025) +
    # gghighlight(class == "inflam") +
    theme_bw() +
    labs(
      caption = paste("p_eqtl = ", p_eqtl),
      title = "Associations between different measures of Skin Color
            and mRNA-Based Disease Signatures, Add Health
            (p-values reported, FDR-corrected for whole genome)",
      y = "mRNA Signatures",
      x = "Skincolor") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          # plot.margin=unit(c(1, 1, 0.1, 1), "cm"),
          text = element_text(size=10, face = "bold")) +
    scale_size_continuous(name = "Adjusted P-value", range = c(0, 21), 
                          limits = c(0.0000001, 100000), breaks = c(0.0000001, 10000, 15000, 25000, 100000),
                          labels = c("n.s.", "p<0.05", "p<0.01", "p<0.001", "p<0.0001")) +
    scale_alpha(guide = 'none') +
    guides(shape = guide_legend(override.aes = list(size = 10)),
           fill = guide_legend(override.aes = list(size = 8)))  
}
plots = p_eqtl %>% map(fig1A)
plots