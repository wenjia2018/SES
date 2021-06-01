
#' ---
#' title: skin color3 
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

#' ### treatments: 
#' * "color_byinterviewer3": "color_byinterviewer3_DarkBlack" and "color_byinterviewer3_LightMed" with "color_byinterviewer3_white" as reference
#' *  Multiple comparison is adjusted using Bonferroni
control = "ancestryPC_ses"

#' ### controls: 
#' 
#' * ancestryPC_ses : basic + ses + bespoke ancestryPC
#' 
#' ### outcome: aging and clusters


#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(here)
library(tidyverse)
source("/home/xu/ses-1/user_wx/extract_v2.R")

data <- readRDS("/home/xu/ses-1/user_wx/color3_bespoke_m1m2m3_12.04.2021.rds")
bespoke = data %>% filter(table1 %>% str_detect("aging|whole_genome"))

n = data %>% filter(table1 %>% str_detect("aging")) %>% pull(table1) %>% unique %>% length


#' ### m1 multivariate outcome = set of individual genes regression (aging, aging_up, aging_down, aging_cluster_complement are not calculated)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
p_eqtl = c(0.01)

out1 = p_eqtl %>%
  map_df(outm1, control, bespoke) %>% 
  mutate(adj.p = ifelse(p*(n-4) < 1, p*(n-4), 1)) %>% 
  filter(adj.p <0.05) %>% 
  mutate(model = "m1(multivariate)")


#' ### m2 univariate outcome (median + Yeo-Johnson)  regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
out2 = p_eqtl %>%
  map_df(outm2, control, bespoke) %>% 
  mutate(adj.p = ifelse(p*n < 1, p*n, 1)) %>% 
  filter(adj.p <0.05) %>% 
  mutate(model = "m2(median)")


#' ### m3 univariate outcome (mean)   regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE

out3 = p_eqtl %>%
  map_df(outm3, control, bespoke) %>% 
  mutate(adj.p = ifelse(p*n < 1, p*n, 1)) %>% 
  filter(adj.p <0.05) %>% 
  mutate(model = "m3(mean)")

data = out1 %>% 
  rbind(out2) %>% 
  rbind(out3)


ex0_m123 = data %>% 
  # rbind(race) %>%
  filter(gene_set_name %>% str_detect("aging")) %>% 
  filter(gene_set_name!="aging_down_cl1_mRNA") %>% 
  dplyr::select(treatment, gene_set_name, adj.p, model) %>% 
  rename(p = adj.p) %>% 
  mutate(pval=case_when(
    # p<0.0001 ~ 0.0001,
    p<0.001 ~ 0.001,
    p<0.01 ~ 0.01,
    p<0.05 ~ 0.05,
    p>0.05 ~ 100),
    pval2=case_when(
      # p<0.0001 ~ 100000,
      p<0.001 ~ 55000,
      p<0.01 ~ 25000,
      p<0.05 ~ 10000,
      p>0.05 ~ 0.0000001),
    treatment= case_when(treatment == "raceethnicity_Hispanic" ~ "Hispanic",
                         treatment =="raceethnicity_NonHblack" ~  "Non Hispanic Black",
                         treatment =="color_byinterviewer3_LightMed"  ~ "Light Medium",  
                         treatment =="color_byinterviewer3_DarkBlack"  ~ "Dark Black")
  )  %>%
  mutate(treatment = factor(treatment, levels = c("Hispanic", "Non Hispanic Black","Light Medium","Dark Black"
  ))) %>%
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
                                  TRUE ~ gene_set_name
         ))
# need to specify shape, because some shape does not have borders. 
# reorder y axis
ex0_m123$gene_set_name <- factor(ex0_m123$gene_set_name, levels = c("1KI", "CTRA", "Antbintf",  "Interferon", "Inflamation in CTRA",
                                                                      
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

# 
# ex0_m123$gene_set_name <- factor(ex0_m123$gene_set_name,
#                                   levels = c("1KI", "CTRA", "Antbintf",  "Interferon", "Inflamation in CTRA",                                                                    "Aging", "Aging Down", "Aging Up",
#                                              "Aging Down Cl1",  "Aging Down Cl1a", "Aging Down Cl1b", "Aging Down Cl1c",
#                                              "Aging Down Cl2",  "Aging Down Cl3",
#                                              "Aging Up Cl1",    "Aging Up Cl2",    "Aging Up Cl3",    "Aging Up Cl4"  )) %>% fct_rev
axiscolor = c("darkblue","darkblue","darkblue","darkblue","grey30", "grey30", "grey30", "grey30", "grey30", "grey30","grey30", "grey30", "grey30")
ex0_m123 %>% 
  # mutate(group = ifelse(treatment %>% str_detect("Hispanic"), "Race", "Skin Color")) %>% 
  ggplot(aes(treatment, gene_set_name, size = pval2
             # fill = `1KI Genes`,
             # colour = model
  )) +
  geom_point(stroke = 1.5, shape = 21, alpha = 0.4, colour = "darkblue", fill = "navy") +
  facet_wrap( ~ model, scales = "free") + 
  # facet_wrap(~group, scales = "free_x") +
  # scale_fill_manual(values = c("darkblue", "goldenrod3")) +
  # scale_color_manual(values = c("darkblue", "goldenrod3")) +
  # geom_jitter(height = 0.00000025) +
  # gghighlight(class == "inflam") +
  theme_bw() +
  labs(
    caption = paste("p_eqtl = ", p_eqtl),
    # title = "Figure 1. Associations between Race Ethnicity/Skin Color
    #         and mRNA-Based Disease Signatures, Add Health
    #         (p-values reported, FDR-corrected for whole genome)",
    y = "mRNA Signatures",
    x = "Skincolor") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        # plot.margin=unit(c(1, 1, 0.1, 1), "cm"),
        axis.text.y = element_text(color=axiscolor),
        text = element_text(size=10, face = "bold")) +
  scale_size_continuous(name = "Adjusted P-value", range = c(0, 21), 
                        limits = c(0.0000001, 55000), breaks = c(0.0000001, 10000, 25000, 55000),
                        labels = c("n.s.", "p<0.05", "p<0.01", "p<0.001")) +
  scale_alpha(guide = 'none') +
  guides(shape = guide_legend(override.aes = list(size = 10)),
         fill = guide_legend(override.aes = list(size = 8)),
         size = FALSE)  

