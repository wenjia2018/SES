#' ---
#' title: check sig genes set in fig1 and fig2, same, overlap?
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

#' #### overlap
#' significant genes from Figure1(FDR corrected regression within each signature)
#'  and Figure2(within in each race strata skincolor responsive significant genes)
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

p_eqtl <- c(0.05, 0.01)
# example_race = readRDS("~/ses-1/user_wx/race_bespoke_12.02.2021.rds")
# example_skincolor3 = readRDS("~/ses-1/user_wx/color3_bespoke_18.02.2021.rds")
# example_black = readRDS("/home/xu/ses-1/user_wx/color3_bespoke_NonHblack_strata_25.02.2021.rds")

example_race = readRDS("~/ses-1/user_wx/race_bespoke_15.03.2021.rds")
example_skincolor3 = readRDS("~/ses-1/user_wx/color3_bespoke_15.03.2021.rds")
example_black = readRDS("/home/xu/ses-1/user_wx/color3_bespoke_NonHblack_strata_16.03.2021.rds")
gene_set_check = function(p_eqtl){
  race = outm8_allsig (p = p_eqtl, control = "ancestryPC_ses", example_race) 
  skincolor = outm8_allsig (p = p_eqtl, control = "ancestryPC_ses", example_skincolor3)
  fig1_sig = race %>% 
    rbind(skincolor) %>%
    filter(gene_set_name %>% str_detect("aging")) %>% 
    filter(gene_set_name!="aging_down_cl1_mRNA") %>% 
    select(p_eqtl, treatment, gene_set_name, gene_sig) %>% 
    mutate(no_sig = gene_sig %>% map_int(~ length(.x))) %>% 
    filter(no_sig>0)
  
  
  
  
  black = outm8_allsig (p = p_eqtl, control = "ancestryPC_ses", example_black)
  fig2_sig = black %>% 
    filter(gene_set_name %>% str_detect("aging")) %>% 
    filter(gene_set_name!="aging_down_cl1_mRNA") %>% 
    select(p_eqtl, treatment, gene_set_name, gene_sig) %>% 
    mutate(no_sig = gene_sig %>% map_int(~ length(.x))) %>% 
    filter(no_sig>0)
 
   par(mar= c(4,4,4,2))
  venn::venn(list(`Aging Down Fig1` = fig1_sig %>%
                    filter(gene_set_name=="aging_down_mRNA") %>% 
                    pull(gene_sig) %>% `[[`(1), 
                  `Aging Down Fig2` = fig2_sig %>%
                    filter(gene_set_name=="aging_down_mRNA") %>% 
                    pull(gene_sig) %>% `[[`(1)),
             ilabels = TRUE, box = F, par = F,
             zcolor = "style", ellipse = FALSE,
             opacity = 0.15)
  title(main = paste("p_eqtl = ", p_eqtl))
  
  par(mar= c(4,4,4,2))
  venn::venn(list(`Aging Fig1` = fig1_sig %>%
                    filter(gene_set_name=="aging_mRNA") %>% 
                    pull(gene_sig) %>% `[[`(1), 
                  `Aging Fig2` = fig2_sig %>%
                    filter(gene_set_name=="aging_mRNA") %>% 
                    pull(gene_sig) %>% `[[`(1)),
             ilabels = TRUE, box = F, par = F,
             zcolor = "style", ellipse = FALSE,
             opacity = 0.15)

  title(main = paste("p_eqtl = ", p_eqtl))
  
  par(mar= c(4,4,4,2))
  venn::venn(list(`Aging Cluster Complement Fig1` = fig1_sig %>%
                    filter(gene_set_name=="aging_cluster_complement") %>% 
                    pull(gene_sig) %>% `[[`(1), 
                  `Aging Cluster Complement Fig2` = fig2_sig %>%
                    filter(gene_set_name=="aging_cluster_complement") %>% 
                    pull(gene_sig) %>% `[[`(1)),
             ilabels = TRUE, box = F, par = F,
             zcolor = "style", ellipse = FALSE,
             opacity = 0.15)
  
  title(main = paste("p_eqtl = ", p_eqtl))
}

p_eqtl %>% map(gene_set_check)
