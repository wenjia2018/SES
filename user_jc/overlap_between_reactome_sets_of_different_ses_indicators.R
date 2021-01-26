#' ---
#' title: Examples
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#' Set global options 
#+ setup, warning=FALSE, message=FALSE
# knitr::opts_chunk$set(echo = FALSE)

library(tidyverse)
library(here)

gsea_genesetnames <- readRDS(here("gsea_genesetnames.rds"))
complete_tables <- readRDS(here("gsea_removefig1A.rds"))
complete_tables <- complete_tables %>% reduce(rbind) %>% select(geneSet, description) %>% unique()
gsea_genesetnames_complement <- map(gsea_genesetnames, ~setdiff(complete_tables$geneSet, .x)) # the complement of these reactome terms

# the universe of reactome terms considered here is "all reactome sets deemed significantly related to at least one ses predictor"
# there are then 2^5 intersections in the venn diagram (i.e. intersections over the 5 ses predictors, with each being the complement or not)
# 

universe= 
  list(gsea_genesetnames_complement,
       gsea_genesetnames) %>% 
  transpose()

crossing(edu =0:1, 
         income =0:1, 
         SEI =0:1,
         ses4 =0:1,
         sss =0:1) %>%
  mutate(x = 
           pmap(., 
                function(edu, income, SEI, ses4, sss) 
                  list(universe$edu[edu + 1], 
                       universe$income[income + 1], 
                       universe$SEI[SEI + 1],
                       universe$ses4[ses4 + 1], 
                       universe$sss[sss + 1]))  %>%
           map(flatten)  %>% 
           map(reduce, intersect)) %>% 
  unnest(x) %>% 
  left_join(complete_tables, by = c("x" = "geneSet"))  %>% 
  knitr::kable()

