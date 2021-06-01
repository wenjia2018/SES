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
control = "ancestryPC_ses"

logFC_ploting = function(example, p_eqtl){
  isnt_out_z <- function(x, thres = 3, na.rm = TRUE) {
    abs(x - mean(x, na.rm = na.rm)) <= thres * sd(x, na.rm = na.rm)
  }


  outm8 = p_eqtl %>%
    f0(data = example) %>% 
    filter(control_set==control) %>% 
    hoist(out, m = list("result", "m8_fdr", 1, "other", "m"))
  
  a = outm8 %>%
    select(p_eqtl, treatment, gene_set_name, m) %>%
    filter(m!="NULL") %>% 
    filter(gene_set_name %>% str_detect("aging")) %>% 
    filter(gene_set_name!="aging_down_cl1_mRNA") %>% 
    mutate_at(.vars = vars("m"),
              .funs = list(~ map(., ~ mutate(., logFC_outlier = isnt_out_z(logFC)) %>% filter(logFC_outlier==TRUE)))) %>% 
    mutate(test = m %>% map(~ wilcox.test(.$logFC, mu = 0, alternative = "two.sided")),
           logFC.t.test.p = test %>% map_dbl(~ .$p.value) %>% p.adjust(method = "fdr") 
           # %>% format(digits = 3, scientific =T)
    ) %>% 
    unnest(m) %>% 

    dplyr::select(treatment, gene_set_name, logFC, logFC.t.test.p) %>% 
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
                                                            "Lysosome Metabolism",#/Glycosaminoglycan Degradation"
                                                            "Fatty Acid Metabolism",#/Peroxisome Activity", 
                                                            "Actin Cytoskeleton",#/Focal Adhesion/Tight Junctions",
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

  axiscolor = c("darkblue","darkblue","darkblue", "darkblue", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30","grey30", "grey30", "grey30")
  a %>% 
    mutate(p_sig = ifelse(logFC.t.test.p>0.05, "Non Sig", "Sig")) %>% 
    ggplot(aes( x = logFC, y = gene_set_name,
                color = p_sig,
                fill = p_sig)) +
    geom_violin()+
    facet_wrap( ~  treatment , scales = "free_x", strip.position = "bottom")  +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
    scale_color_manual(values=c("goldenrod3", "lightblue"))+
    scale_fill_manual(values=c("goldenrod3", "lightblue"))+
    theme(
      # legend.position = "none",
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
#' ### skincolor3 (p_eqtl = 0.01)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example = readRDS("/home/xu/ses-1/user_wx/color3_bespoke_28.03.2021.rds")
example = example %>% filter(table1 %>% str_detect("aging|whole_genome"))

p_eqtl = 0.01
logFC_ploting(example = example, p_eqtl = p_eqtl)

#' ### skincolor3 (p_eqtl = 0.05)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

p_eqtl = 0.05

logFC_ploting(example = example, p_eqtl = p_eqtl)

#' ## skincolor3 in black race strata gene by gene regression within each signature logFC plots
#' ### skincolor3 in black race strata (p_eqtl = 0.01)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
example = readRDS("/home/xu/ses-1/user_wx/color3_bespoke_NonHblack_strata_28.03.2021.rds")
example= example %>% filter(table1 %>% str_detect("aging|whole_genome"))
p_eqtl = 0.01
logFC_ploting(example = example, p_eqtl = p_eqtl)

#' ### skincolor3 in black race strata (p_eqtl = 0.05)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

p_eqtl = 0.05
logFC_ploting(example = example, p_eqtl = p_eqtl)
