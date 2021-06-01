#' ---
#' title: skincolor3 Aging signature logFC correlated t test
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

#' ## skincolor3 gene by gene regression within each signature logFC t-test significant at maximum correlation
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
isnt_out_z <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - mean(x, na.rm = na.rm)) <= thres * sd(x, na.rm = na.rm)
}

enf <-function(trimmedLogFCs){
  r = seq(0.1, 0.9, by = 0.1)
  n = length(trimmedLogFCs)
  mu = mean(trimmedLogFCs)
  stdv = sd(trimmedLogFCs)
  t = map_dbl(r, function(rho) mu/sqrt(stdv^2 * (1/n + rho/(1- rho))))
  thr_left = qt(0.025, n-1)
  thr_right = qt(0.975, n-1)
  
  p = enframe(t, value = "t.statistic") %>% 
    transmute(cor = r,
              # sig = ifelse(between(t.statistic, thr_left, thr_right), "Non Sig", "Sig"),
              p = 2*pt(-abs(t.statistic), df = n-1))
}
logFC_test = function(example, p_eqtl){

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
    mutate(test_statistics = m %>% map(~ enf(.$logFC))) %>% 
    unnest(test_statistics) %>% 
    select(-m) %>% 
    group_by(cor) %>% 
    mutate(p.adj =p.adjust(p,method= "fdr")) %>% 
    ungroup %>% 
    filter(p.adj<0.05) %>% 
    group_by(treatment, gene_set_name) %>% 
    slice_max(cor) %>% 
    ungroup %>% 
    dplyr::select(treatment, gene_set_name, cor, p.adj) %>% 
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
}
#' ### skincolor3 (p_eqtl = 0.01)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example = readRDS("/home/xu/ses-1/user_wx/color3_bespoke_28.03.2021.rds")
example = example %>% filter(table1 %>% str_detect("aging|whole_genome"))

p_eqtl = 0.01
logFC_test(example = example, p_eqtl = p_eqtl) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()

#' ### skincolor3 (p_eqtl = 0.05)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

p_eqtl = 0.05

logFC_test(example = example, p_eqtl = p_eqtl) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()

#' ## skincolor3 in black race strata gene by gene regression within each signature logFC t-test significant at maximum correlation
#' ### skincolor3 in black race strata (p_eqtl = 0.01)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
example = readRDS("/home/xu/ses-1/user_wx/color3_bespoke_NonHblack_strata_28.03.2021.rds")
example= example %>% filter(table1 %>% str_detect("aging|whole_genome"))
p_eqtl = 0.01
logFC_test(example = example, p_eqtl = p_eqtl) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()

#' ### skincolor3 in black race strata (p_eqtl = 0.05)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

p_eqtl = 0.05
logFC_test(example = example, p_eqtl = p_eqtl) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()
