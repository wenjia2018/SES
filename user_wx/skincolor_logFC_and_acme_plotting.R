
#' ---
#' title: logFC and aceme plotting aging composite ancestry control with fdr correction
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE


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
library(stringi)
walk(dir(path = here("R"),full.names = TRUE), source)
source("/home/xu/ses-1/user_wx/extract_v2.R")
load_data(reconciled = FALSE, remove_inflam = FALSE)
source("/home/xu/ses-1/user_wx/plotting_utils.R")
control = "ancestryPC_ses"
# skincolor_eqtl005_aging_composite_ancestry_11.05.2021.rds m7 m8 without genowide mediation
example_skincolor3 = readRDS("/home/xu/ses-1/user_wx/skincolor_eqtl005_aging_composite_ancestry_11.05.2021.rds")

p_eqtl = 0.05


#' ### 1 omnibus logFC with aging composite ancestry controls
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=10, fig.height=8

logfc = logFC_ploting(example = example_skincolor3, p_eqtl = p_eqtl)

logfc

logfc_sig = logfc$data %>% 
  group_nest(treatment, gene_set_name) %>% 
  hoist(data, psig = list("p_sig")) %>% 
  mutate(psig = psig %>% map_chr(unique))

#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=10, fig.height=8
# skincolor3_wholegenome_aging_composite_ancestry_14.05.2021.rds genowide mediation
example_genowide_med_skincolor3 = readRDS("/home/xu/ses-1/user_wx/skincolor3_wholegenome_aging_composite_ancestry_14.05.2021.rds")
gene_sets= c("aging_mRNA",
             "aging_up_mRNA",
             "aging_down_mRNA",
             "aging_cluster_complement_mRNA",
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

mediators = 
  c(
    "stress_perceived_lm",
    "bills_binary",
    "currentsmoke_binary",
    "w5bmi_lm",
    "insurance_lack_binary",
    "totdiscrim2_binary",
    "discrim2_binary",  
    "totdiscrim1_category")
temp =
  example_genowide_med_skincolor3$out[[1]]$result$example0 %>% 
  hoist(out, med= list("result", "m8_fdr", 1, "mediation_single")) %>%
  unnest_longer(med) %>% 
  unnest_longer(med) %>% 
  hoist(med, acme_p = list("result", "p")) %>%
  hoist(med, prop = list("result", "other", "med_prop")) %>%
  hoist(med, acme = list("result", "other", "med_ACME")) %>%
  hoist(med, ade = list("result", "other", "med_ADE")) %>% 
  hoist(med, ade_p = list("result", "other", "med_ADE_p")) %>% 
  dplyr::select(-controls, -med, -out) 
# old way to split mediator and gene
# temp1 =
#   temp %>% 
#   dplyr::mutate(mediator=stri_extract(med_id, regex='[^_]*'),
#                 gene =stri_extract_last(med_id, regex='([^_]+$)')) %>% 
#   mutate(treatment= case_when(treatment =="color_byinterviewer3_LightMed" ~ "Light Medium",
#                               treatment =="color_byinterviewer3_DarkBlack" ~ "Dark Black") %>% 
#            factor(levels = c("Light Medium","Dark Black")))

# another way to split mediator and gene
# temp1 =
#   temp %>% 
#   dplyr::mutate(mediator=stri_extract(med_id, regex='^([^_]*_[^_]*)'),
#                 gene =str_remove(med_id, str_c(mediator,"_"))) %>% 
#   mutate(treatment= case_when(treatment =="color_byinterviewer3_LightMed" ~ "Light Medium",
#                               treatment =="color_byinterviewer3_DarkBlack" ~ "Dark Black") %>% 
#            factor(levels = c("Light Medium","Dark Black")))


temp1 =
  temp %>% 
  dplyr::mutate(gene_start = (med_id  %>% stri_locate_first_regex( "[A-Z]"))[,1],
                gene = med_id %>% str_sub(gene_start),
                mediator = med_id %>% str_sub(1, gene_start-2)) %>% 
  mutate(treatment= case_when(treatment =="color_byinterviewer3_LightMed" ~ "Light Medium",
                              treatment =="color_byinterviewer3_DarkBlack" ~ "Dark Black") %>% 
           factor(levels = c("Light Medium","Dark Black")))



# after genowide correction there is no significant mediation results
  # temp1 %>% 
  # group_by(treatment) %>% 
  # mutate(across(.cols = c(acme_p, ade_p), .fns = ~ p.adjust(.x, method = "fdr"), .names = "{.col}_adj")) %>% 
  # ungroup %>% 
  # filter(acme_p_adj<0.05)

temp_data = 
  tibble(gene_set_name = gene_sets,
         data = gene_sets %>% 
           map(~ filter(temp1, gene %in% signatures$outcome_set[[.]]) %>% 
                 select(-gene_set_name))) %>% 
  unnest(data) %>% 
  mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
         gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
         prop= as.numeric(prop),
         acme= as.numeric( acme),
         ade = as.numeric(ade),
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
temp_data1 =
  temp_data %>% 
  group_by(gene_set_name, treatment, mediator) %>% 
  mutate(across(c(prop, acme, ade),  ~ out_lier_NA(.x))) %>% 
  ungroup %>% 
  group_nest(gene_set_name, treatment) %>% 
  left_join(logfc_sig %>% select(1:3), by = c("gene_set_name", "treatment")) %>% 
  unnest(data) %>% 
  group_nest(gene_set_name, treatment, mediator) %>%
  mutate(ade_test = data %>% map_dbl(~ wilcox.test(.$ade, mu = 0, alternative = "two.sided")$p.value),
         acme_test = data %>% map_dbl(~ wilcox.test(.$acme, mu = 0, alternative = "two.sided")$p.value),
         prop_test= data %>% map_dbl(~ wilcox.test(.$prop, mu = 0, alternative = "two.sided")$p.value)) %>% 
  mutate(across(c(prop_test, acme_test, ade_test),  ~ p.adjust(.x, method = "fdr"))) %>% 
  unnest(data) %>% 
  mutate(ACME_test = ifelse(psig=="Sig", ifelse(acme_test< 0.05, "ACME sig", "ACME not sig"), "Total effect not sig"))


#' ### 2 omnibus mediation with aging composite ancestry controls acme plotting
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=10, fig.height=8
acme_plotting = function(data) {
  axiscolor = c("darkblue","darkblue","darkblue","darkblue","grey30", "grey30", "grey30", "grey30", "grey30", "grey30","grey30", "grey30", "grey30")
  
data %>% 
    ggplot(aes( x = acme, y = gene_set_name, fill = ACME_test, color = ACME_test)) +
    geom_violin() +
    facet_wrap( ~  treatment , scales = "free_x", strip.position = "bottom")  +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
    scale_color_manual(values=c("goldenrod3", "lightblue", "white"))+
    scale_fill_manual(values=c("goldenrod3", "lightblue", "white"))+
    labs(
      title = paste("mediator = ", data$mediator %>% unique),
      # title = "Figure 1. Associations between Race Ethnicity/Skin Color
      #         and mRNA-Based Disease Signatures, Add Health
      #         (p-values reported, FDR-corrected for whole genome)",
      y = "mRNA Signatures",
      x = "Skincolor") +
    theme(
      legend.position = "bottom",
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
}


acme_plotting(temp_data1 %>% filter(mediator == mediators[1]))
acme_plotting(temp_data1 %>% filter(mediator == mediators[2]))
acme_plotting(temp_data1 %>% filter(mediator == mediators[3]))
acme_plotting(temp_data1 %>% filter(mediator == mediators[4]))
acme_plotting(temp_data1 %>% filter(mediator == mediators[5]))
acme_plotting(temp_data1 %>% filter(mediator == mediators[6]))
acme_plotting(temp_data1 %>% filter(mediator == mediators[7]))
acme_plotting(temp_data1 %>% filter(mediator == mediators[8]))
