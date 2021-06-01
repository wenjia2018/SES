#' ---
#' title: skincolor3 Aging signature logFC plotting
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

#' ## skincolor3 gene by gene regression within each signature logFC plots
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
load_data(reconciled = FALSE, remove_inflam = FALSE)
source("/home/xu/ses-1/user_wx/extract_v2.R")

#' ### skincolor3 (p_eqtl = 0.01)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
example = readRDS("/home/xu/ses-1/user_wx/color3_bespoke_28.03.2021.rds")
example_skincolor3 = example %>% filter(table1 %>% str_detect("aging|whole_genome"))
p_eqtl = 0.01
control = "ancestryPC_ses"

bespoke = example_skincolor3
outm8 = p_eqtl %>%
  f0(data = example_skincolor3) %>% 
  filter(control_set==control) %>% 
  hoist(out, m = list("result", "m8_fdr", 1, "other", "m"))

a = outm8 %>%
  unnest(m) %>% 
    filter(gene_set_name %>% str_detect("aging")) %>% 
    filter(gene_set_name!="aging_down_cl1_mRNA") %>% 
    dplyr::select(treatment, gene_set_name, logFC) %>% 
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
                                                            "Aging Cluster Complement",
                                                            "Aging Composite"
           )) %>% fct_rev
           # gene_set_name = str_c(gene_set_name, "_", p_id)
    )
axiscolor = c("darkblue","darkblue","darkblue", "darkblue", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30","grey30", "grey30", "grey30")
a %>% 
  ggplot(aes( x = logFC, y = gene_set_name, color = treatment, fill = p_val)) +
  geom_violin()+
  facet_grid( ~  treatment , scales = "free")  +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")+
  scale_color_manual(values=c("goldenrod3", "lightblue"))+
  theme(legend.position = "none",
        panel.margin.x = unit(1, "lines"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color=axiscolor),
        axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
        axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
        plot.title = element_text(size = 20, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
        plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
        strip.text = element_text(size = 7),
        text = element_text(family = "Georgia"))

#' ### logFC double side t test p value 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
two_side_t = outm8 %>%
  select(p_eqtl, treatment, gene_set_name, m) %>%
  filter(m!="NULL") %>% 
  mutate(test = m %>% map(~ t.test(.$logFC, mu = 0, alternative = "two.sided")),
         logFC.t.test.p = test %>% map_dbl(~ .$p.value) 
         # %>% format(digits = 3, scientific =T)
         ) %>% 
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

two_side_t %>% 
  select(-m,-test) %>% 
  arrange(logFC.t.test.p) %>%
  kableExtra::kable() %>% 
  kableExtra::kable_styling()

# dt = 
#   two_side_t %>% 
#   filter(gene_set_name %>% str_detect("aging")) %>% 
#   filter(gene_set_name!="aging_down_cl1_mRNA") %>% 
#   dplyr::select(treatment, gene_set_name, p) %>% 
#   mutate(pval=case_when(
#     # p<0.0001 ~ 0.0001,
#     p<0.001 ~ 0.001,
#     p<0.01 ~ 0.01,
#     p<0.05 ~ 0.05,
#     p>0.05 ~ 100),
#     pval2=case_when(
#       # p<0.0001 ~ 100000,
#       p<0.001 ~ 55000,
#       p<0.01 ~ 25000,
#       p<0.05 ~ 10000,
#       p>0.05 ~ 0.0000001),
#     treatment= case_when(treatment == "raceethnicity_Hispanic" ~ "Hispanic",
#                          treatment =="raceethnicity_NonHblack" ~  "Non Hispanic Black",
#                          treatment =="color_byinterviewer3_LightMed"  ~ "Light Medium",  
#                          treatment =="color_byinterviewer3_DarkBlack"  ~ "Dark Black")
#   )  %>%
#   mutate(treatment = factor(treatment, levels = c("Hispanic", "Non Hispanic Black","Light Medium","Dark Black"
#   ))) %>%
#   mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
#          gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
#          # gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI"),
#          # gene_set_name = gene_set_name %>% replace(gene_set_name == "Ctra", "CTRA"),
#          # gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflame", "Inflamation in CTRA"),
#          gene_set_name= case_when(gene_set_name =="Aging Up Cl1" ~ "Innate/Adaptive Immunity",
#                                   gene_set_name =="Aging Up Cl2" ~ "Actin Cytoskeleton",#/Focal Adhesion/Tight Junctions",
#                                   gene_set_name =="Aging Up Cl3" ~ "Fatty Acid Metabolism",#/Peroxisome Activity",
#                                   gene_set_name =="Aging Up Cl4" ~ "Lysosome Metabolism",#/Glycosaminoglycan Degradation",
#                                   # gene_set_name =="Aging Down Cl1" ~ "Aging Down Cl1",
#                                   gene_set_name =="Aging Down Cl1a" ~ "RNA Metabolism",#/Ribosome Biogenesis/Purine Metabolism",
#                                   gene_set_name =="Aging Down Cl1b" ~ "Mitochondrial",
#                                   gene_set_name =="Aging Down Cl1c" ~ "DNA Replication",#/Elongation/Mismatch Repair",
#                                   gene_set_name =="Aging Down Cl2" ~ "Ribosome",
#                                   gene_set_name =="Aging Down Cl3" ~ "Immune Genes",
#                                   gene_set_name =="Aging" ~ "Aging Composite",
#                                   gene_set_name == "Inflam1k" ~ "1KI",
#                                   gene_set_name == "Ctra" ~ "CTRA",
#                                   gene_set_name == "Inflame" ~ "Inflamation in CTRA",
#                                   TRUE ~ gene_set_name
#          ))
# # need to specify shape, because some shape does not have borders. 
# # reorder y axis
# dt$gene_set_name <- factor(dt$gene_set_name, levels = c("1KI", "CTRA", "Antbintf",  "Interferon", "Inflamation in CTRA",
#                                                         
#                                                         "Aging Down Cl1",
#                                                         "RNA Metabolism",#/Ribosome Biogenesis/Purine Metabolism",
#                                                         "Mitochondrial",
#                                                         "DNA Replication",#/Elongation/Mismatch Repair",
#                                                         "Ribosome",
#                                                         "Immune Genes",
#                                                         "Innate/Adaptive Immunity", 
#                                                         "Actin Cytoskeleton",#/Focal Adhesion/Tight Junctions",
#                                                         "Fatty Acid Metabolism",#/Peroxisome Activity", 
#                                                         "Lysosome Metabolism",#/Glycosaminoglycan Degradation"
#                                                         "Aging Cluster Complement",
#                                                         "Aging Up",
#                                                         "Aging Down",
#                                                         "Aging Composite"
#                                                         
# )) %>% fct_rev
# 
# # 
# # dt$gene_set_name <- factor(dt$gene_set_name,
# #                                   levels = c("1KI", "CTRA", "Antbintf",  "Interferon", "Inflamation in CTRA",                                                                    "Aging", "Aging Down", "Aging Up",
# #                                              "Aging Down Cl1",  "Aging Down Cl1a", "Aging Down Cl1b", "Aging Down Cl1c",
# #                                              "Aging Down Cl2",  "Aging Down Cl3",
# #                                              "Aging Up Cl1",    "Aging Up Cl2",    "Aging Up Cl3",    "Aging Up Cl4"  )) %>% fct_rev
# axiscolor = c("darkblue","darkblue","darkblue","darkblue","grey30", "grey30", "grey30", "grey30", "grey30", "grey30","grey30", "grey30", "grey30")
# dt %>% 
#   # mutate(group = ifelse(treatment %>% str_detect("Hispanic"), "Race", "Skin Color")) %>% 
#   ggplot(aes(treatment, gene_set_name, size = pval2
#              # fill = `1KI Genes`,
#              # colour = `1KI Genes`
#   )) +
#   geom_point(stroke = 1.5, shape = 21, alpha = 0.4, colour = "darkblue", fill = "navy") +
#   # facet_wrap(~group, scales = "free_x") +
#   # scale_fill_manual(values = c("darkblue", "goldenrod3")) +
#   # scale_color_manual(values = c("darkblue", "goldenrod3")) +
#   # geom_jitter(height = 0.00000025) +
#   # gghighlight(class == "inflam") +
#   theme_bw() +
#   labs(
#     caption = paste("p_eqtl = ", p_eqtl),
#     # title = "Figure 1. Associations between Race Ethnicity/Skin Color
#     #         and mRNA-Based Disease Signatures, Add Health
#     #         (p-values reported, FDR-corrected for whole genome)",
#     y = "mRNA Signatures",
#     x = "Skincolor") +
#   theme(axis.text.x = element_text(angle = 30, hjust = 1),
#         # plot.margin=unit(c(1, 1, 0.1, 1), "cm"),
#         axis.text.y = element_text(color=axiscolor),
#         text = element_text(size=10, face = "bold")) +
#   scale_size_continuous(name = "Adjusted P-value", range = c(0, 21), 
#                         limits = c(0.0000001, 55000), breaks = c(0.0000001, 10000, 25000, 55000),
#                         labels = c("n.s.", "p<0.05", "p<0.01", "p<0.001")) +
#   scale_alpha(guide = 'none') +
#   guides(shape = guide_legend(override.aes = list(size = 10)),
#          fill = guide_legend(override.aes = list(size = 8)),
#          size = FALSE)  
#' ### skincolor3 (p_eqtl = 0.05)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

p_eqtl = 0.05
control = "ancestryPC_ses"

bespoke = example_skincolor3
outm8 = p_eqtl %>%
  f0(data = example_skincolor3) %>% 
  filter(control_set==control) %>% 
  hoist(out, m = list("result", "m8_fdr", 1, "other", "m"))

a = outm8 %>%
  unnest(m) %>%
  filter(gene_set_name %>% str_detect("aging")) %>% 
  filter(gene_set_name!="aging_down_cl1_mRNA") %>% 
  dplyr::select(treatment, gene_set_name, logFC) %>% 
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
                                                          "Aging Cluster Complement",
                                                          "Aging Composite"
         )) %>% fct_rev
         # gene_set_name = str_c(gene_set_name, "_", p_id)
  )
axiscolor = c("darkblue","darkblue","darkblue", "darkblue", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30","grey30", "grey30", "grey30")
a %>% 
  ggplot(aes(x = gene_set_name, y = logFC, fill = treatment)) +
  geom_violin()+
  facet_wrap( ~ treatment, scales = "free")  +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")+
  scale_fill_manual(values=c("goldenrod3", "lightblue"))+
  coord_flip()+
  theme(legend.position = "none",
        panel.margin.x = unit(1, "lines"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color=axiscolor),
        axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
        axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
        plot.title = element_text(size = 20, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
        plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
        strip.text = element_text(size = 7),
        text = element_text(family = "Georgia"))

#' ### logFC double side t test p value 
#+ echo=F, eval=T, warning=FALSE, message=FALSE

two_side_t = outm8 %>%
  select(p_eqtl, treatment, gene_set_name, m) %>%
  filter(m!="NULL") %>% 
  mutate(test = m %>% map(~ t.test(.$logFC, mu = 0, alternative = "two.sided")),
         logFC.t.test.p = test %>% map_dbl(~ .$p.value) 
         # %>% format(digits = 3, scientific =T)
  ) %>% 
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


two_side_t %>% 
  select(-m,-test) %>% 
  arrange(logFC.t.test.p) %>%
  kableExtra::kable() %>% 
  kableExtra::kable_styling()


#' ## skincolor3 in black race strata gene by gene regression within each signature logFC plots
#' ### skincolor3 in black race strata (p_eqtl = 0.01)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
example = readRDS("/home/xu/ses-1/user_wx/color3_bespoke_NonHblack_strata_28.03.2021.rds")
example_skincolor3 = example %>% filter(table1 %>% str_detect("aging|whole_genome"))
p_eqtl = 0.01
control = "ancestryPC_ses"

bespoke = example_skincolor3
outm8 = p_eqtl %>%
  f0(data = example_skincolor3) %>% 
  filter(control_set==control) %>% 
  hoist(out, m = list("result", "m8_fdr", 1, "other", "m"))

a = outm8 %>%
  unnest(m) %>% 
  filter(gene_set_name %>% str_detect("aging")) %>% 
  filter(gene_set_name!="aging_down_cl1_mRNA") %>% 
  dplyr::select(treatment, gene_set_name, logFC) %>% 
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
                                                          "Aging Cluster Complement",
                                                          "Aging Composite"
         )) %>% fct_rev
         # gene_set_name = str_c(gene_set_name, "_", p_id)
  )
axiscolor = c("darkblue","darkblue","darkblue", "darkblue", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30","grey30", "grey30", "grey30")
a %>% 
  ggplot(aes(x = gene_set_name, y = logFC, fill = treatment)) +
  geom_violin()+
  facet_wrap( ~ treatment, scales = "free")  +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")+
  scale_fill_manual(values=c("goldenrod3", "lightblue"))+
  coord_flip()+
  theme(legend.position = "none",
        panel.margin.x = unit(1, "lines"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color=axiscolor),
        axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
        axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
        plot.title = element_text(size = 20, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
        plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
        strip.text = element_text(size = 7),
        text = element_text(family = "Georgia"))


#' ### logFC double side t test p value 
#+ echo=F, eval=T, warning=FALSE, message=FALSE

two_side_t = outm8 %>%
  select(p_eqtl, treatment, gene_set_name, m) %>%
  filter(m!="NULL") %>% 
  mutate(test = m %>% map(~ t.test(.$logFC, mu = 0, alternative = "two.sided")),
         logFC.t.test.p = test %>% map_dbl(~ .$p.value) 
         # %>% format(digits = 3, scientific =T)
  ) %>% 
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


two_side_t %>% 
  select(-m,-test) %>% 
  arrange(logFC.t.test.p) %>%
  kableExtra::kable() %>% 
  kableExtra::kable_styling()

#' ### skincolor3 in black race strata (p_eqtl = 0.05)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

p_eqtl = 0.05
control = "ancestryPC_ses"

bespoke = example_skincolor3
outm8 = p_eqtl %>%
  f0(data = example_skincolor3) %>% 
  filter(control_set==control) %>% 
  hoist(out, m = list("result", "m8_fdr", 1, "other", "m"))

a = outm8 %>%
  unnest(m) %>% 
  filter(gene_set_name %>% str_detect("aging")) %>% 
  filter(gene_set_name!="aging_down_cl1_mRNA") %>% 
  dplyr::select(treatment, gene_set_name, logFC) %>% 
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
                                                          "Aging Cluster Complement",
                                                          "Aging Composite"
         )) %>% fct_rev
         # gene_set_name = str_c(gene_set_name, "_", p_id)
  )
axiscolor = c("darkblue","darkblue","darkblue", "darkblue", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30","grey30", "grey30", "grey30")
a %>% 
  ggplot(aes(x = gene_set_name, y = logFC, fill = treatment)) +
  geom_violin()+
  facet_wrap( ~ treatment, scales = "free")  +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")+
  scale_fill_manual(values=c("goldenrod3", "lightblue"))+
  coord_flip()+
  theme(legend.position = "none",
        panel.margin.x = unit(1, "lines"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color=axiscolor),
        axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
        axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
        plot.title = element_text(size = 20, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
        plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
        strip.text = element_text(size = 7),
        text = element_text(family = "Georgia"))

#' ### logFC double side t test p value 
#+ echo=F, eval=T, warning=FALSE, message=FALSE

two_side_t = outm8 %>%
  select(p_eqtl, treatment, gene_set_name, m) %>%
  filter(m!="NULL") %>% 
  mutate(test = m %>% map(~ t.test(.$logFC, mu = 0, alternative = "two.sided")),
         logFC.t.test.p = test %>% map_dbl(~ .$p.value) 
         # %>% format(digits = 3, scientific =T)
  ) %>% 
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


two_side_t %>% 
  select(-m,-test) %>% 
  arrange(logFC.t.test.p) %>%
  kableExtra::kable() %>% 
  kableExtra::kable_styling()
