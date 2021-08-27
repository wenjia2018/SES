#' ---
#' title: DE and tfbm for skin color
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
library(readxl)
library(parsnip)
library(workflows)
library(Biobase)
library(enrichplot)
library(gridExtra)
library(grid)
library(ggpubr)
library("clusterProfiler")
walk(dir(path = here("R"),full.names = TRUE), source)
load_data(reconciled = FALSE, remove_inflam = FALSE)
source("/home/xu/ses-1/user_wx/extract_v2.R")
# example_race = readRDS("~/ses-1/user_wx/race_bespoke_15.03.2021.rds")
example = readRDS("/home/xu/ses-1/user_wx/color3_bespoke_28.03.2021.rds")

example_skincolor3 = example %>% filter(table1 %>% str_detect("aging|whole_genome"))
control = "ancestryPC_ses"
p_eqtl = 0.05


#' ### gsea
#+ echo=F, eval=T, warning=FALSE, message=FALSE

all = read.gmt("./msigdb.v7.4.symbols.gmt")
path = split(all$gene,all$term)

DE_ttT = p_eqtl %>%
  outm10_whole_genome(., control, data = example_skincolor3) %>% 
  select(-tfbm_all, -tfbm_immue)

make_list = function(ttT) {
  delist = ttT$logFC
  names(delist) = ttT$gene
  delist = sort(delist,decreasing = T)
}



DE_ttT = DE_ttT %>% 
 mutate(delist = ttT %>% map(~ make_list(.x)))


set.seed(123456)
skincolor_fgsea = DE_ttT %>% 
  mutate(gsea = delist %>% map(~ fgsea::fgseaMultilevel(pathways = path, stats = .x, nPermSimple = 10000)))

skincolor_fgsea %>% saveRDS("./user_wx/fgsea_skincolor.rds")


