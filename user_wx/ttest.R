#' ---
#' title: Supplementary Material
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE, include=FALSE
#+ echo=FALSE

set.seed(123)
library(here)
library(tidyverse)
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase)
library(enrichplot)
library(stringi)
library(dbr) # my package
library(SASxport)
#'## Descriptives
example0 = readRDS("/home/share/projects/aging/example0_mediation_all_5_tot_4.rds")

ex2 <- example0 %>%
  hoist(out, p = list("result", "m8_fdr", 1, "other", "m")) %>% 
  filter(controls=="all") %>% 
  dplyr::select(treatment, gene_set_name, p) %>% 
  dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  mutate(treatment= case_when(treatment =="work_collar_ff5" ~ "Occupation Work Collar",
                              treatment =="edu_max" ~ "Education" ,
                              treatment =="income_hh_ff5" ~ "Income"     ,
                              treatment =="SEI_ff5" ~ "Occupation"      ,
                              treatment =="ses_sss_composite" ~ "SES Composite"  ,
                              treatment =="sss_5" ~ "Subjective Social Status",
                              treatment =="ses_composite_ff5"  ~ "SES Composite 3")) %>% 
  mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
         gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
         gene_set_name= case_when(gene_set_name =="Aging Up Cl4" ~ "Lysosome Metabolism" ,
                                  gene_set_name =="Aging Up Cl1" ~ "Innate/Adaptive Immunity",
                                  gene_set_name =="Aging Up Cl3" ~ "Fatty Acid Metabolism" ,
                                  gene_set_name =="Aging Up Cl2" ~ "Actine Regulation",
                                  gene_set_name =="Aging Down Cl3" ~ "Immune Genes",
                                  gene_set_name =="Aging Down Cl2" ~ "Ribosome"  ,
                                  gene_set_name =="Aging Down Cl1b" ~ "Mitochondrial",
                                  gene_set_name =="Aging Down Cl1a" ~ "RNA Metabolism",
                                  gene_set_name =="Aging Down Cl1c" ~ "DNA Replication",
                                  gene_set_name =="Aging Cluster Complement" ~ "Aging Complement",
                                  gene_set_name =="Aging" ~ "Aging All Genes"))
#outlier detection https://www.r-bloggers.com/2017/12/combined-outlier-detection-with-dplyr-and-ruler/
isnt_out_z <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - mean(x, na.rm = na.rm)) <= thres * sd(x, na.rm = na.rm)
}

exx <- ex2 %>% #group_by(treatment, gene_set_name) %>%
  mutate(t.test=map(.x=ex2$p, .f=~ isnt_out_z(.x$logFC)))

exx <- ex2 %>% 
  unnest(p) %>% 
  group_by(treatment, gene_set_name) %>%
  mutate(t.test=isnt_out_z(logFC)) %>% 
  ungroup %>% 
  dplyr::select(treatment, gene_set_name, logFC, t.test)

exxx <-ex2 %>%
  unnest(p) %>%
  left_join(exx, by=c("treatment", "gene_set_name", "logFC")) %>% 
  dplyr::filter(t.test==T)

# exx<-ex2 %>% #group_by(treatment, gene_set_name) %>%
#   mutate(outlier=map(.x=ex2$, .f=~ (boxplot(.x$logFC)$out))
# boxplot(mtcars$disp)$out
# 
# exx <- exx %>%
#   mutate(tidy = map(wilcoxtest, broom::tidy))
# 
# exx <-exx %>%
#   unnest(tidy, .drop = T) %>% 
#   dplyr::select(treatment, gene_set_name, p, p.value) %>% 
#   mutate(padj.tval=p.adjust(p.value, method="BH"))

library(tidyverse)

#######################################################
#Objective: a cheap and cheerful solution to the following question: does the
#*distribution* of effects across the signature have zero mean.   the validity
#of inference in the setting rests on three independent corrections: correction
#for potential non-normality/skewness/outliers, correction for the possible
#statistical dependence between the gene specific effects logFC estimates, and
#multiple comparisons across the multiple treatment signature combinations. the
#below simply offers a cheap scheme for the second correction. The first and
#third should be relatively straightforward, e.g. Bonferonni or FDR, and remove
#outliers.
#######################################################

enf <-function(trimmedLogFCs){
  r = seq(0.1, 0.9, by = 0.1)
  n = length(trimmedLogFCs)
  mu = mean(trimmedLogFCs)
  stdv = sd(trimmedLogFCs)
  t = map_dbl(r, function(rho) mu/sqrt(stdv^2 * (1/n + rho/(1- rho))))
  thr_left = qt(0.025, n-1)
  thr_right = qt(0.975, n-1)

  enframe(t, value = "t.statistic") %>% 
    mutate(r = r,
           sig = ifelse(between(t.statistic, thr_left, thr_right), "Non Sig", "Sig"),
           p = 2*pt(-abs(t.statistic), df = n-1)) %>% 
    filter(sig=="Sig") %>% 
    slice_max(r) %>% 
    pull(r)
}

  ######################################################## Plot
  #dependence-corrected t-value against correlation against. Horizontal lines
  #indicate the threshold for significance. naturally, because there is less
  #information increasingly correlated data, the t-value and significance
  #decreases on the x-axis (dependence).
  ########################################################


enf_2 <- function(trimmedLogFCs){
  r = seq(0.1, 0.9, by = 0.1)
  n = length(trimmedLogFCs)
  mu = mean(trimmedLogFCs)
  stdv = sd(trimmedLogFCs)
  t = map_dbl(r, function(rho) mu/(stdv^2 * (1/n + rho/(1- rho))))
  thr_pos =  qt(0.025, n-1)
  thr_neg = -qt(0.025, n-1)
  
  ######################################################## Plot
  #dependence-corrected t-value against correlation against. Horizontal lines
  #indicate the threshold for significance. naturally, because there is less
  #information increasingly correlated data, the t-value and significance
  #decreases on the x-axis (dependence).
  ########################################################
  library(ggformula)
  enf<- enframe(t) 
  enf <- enf %>% mutate(p.value = 2*pt(-abs(value), df=n-1),
                        padj.tval=p.adjust(p.value, method="BH"))
  #%>% dplyr::filter(value>thr_pos) %>% dplyr::slice(1)
  return(enf)
}


rob<-exxx %>% 
  group_by(treatment, gene_set_name) %>%
  group_map(~ enf(.x$logFC)) 

#calculation p-value too
robust<-exxx %>% 
  group_by(treatment, gene_set_name) %>%
  group_map(~ enf_2(.x$logFC))

#adjusting the p-value
library(purrr)
map(robust, ~filter(.x, padj.tval>0.05))

