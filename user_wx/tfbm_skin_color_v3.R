#' ---
#' title: DE and tfbm for skin color for eqtl =0.05
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
library(dbr) # my package
walk(dir(path = here("R"),full.names = TRUE), source)
load_data(reconciled = FALSE, remove_inflam = FALSE)
source("/home/xu/ses-1/user_wx/extract_v2.R")
# example_race = readRDS("~/ses-1/user_wx/race_bespoke_15.03.2021.rds")
example = readRDS("/home/xu/ses-1/user_wx/color3_bespoke_28.03.2021.rds")
example_skincolor3 = example %>% filter(table1 %>% str_detect("aging|whole_genome"))
control = "ancestryPC_ses"
p_eqtl = 0.05
tempfolder = "temp_webgestalt"

temp = p_eqtl %>%
  outm10_whole_genome(., control, data = example_skincolor3) %>% 
  filter(ttT!="NULL") %>% 
  mutate(DE_gene = map(ttT, ~ filter(., adj.P.Val< 0.05) %>% pull(gene))) %>%
  filter(map_lgl(DE_gene, ~ length(.)>0))
DE_gene = list2(!!str_c(temp$treatment, temp$p_eqtl) := temp %>% pull(DE_gene)) %>% 
  unlist(recursive = F)

# temp$ttT[[1]] %>% select(1) %>% write_csv("./user_wx/degene_darkblack005.csv")
# pass the genelist to gorilla got GO_ldarkblack005.xlsx


GO_db = read_excel("/home/xu/ses-1/user_wx/GO_darkblack005.xlsx")

temp2 = p_eqtl %>%
  outm10_whole_genome(., control, data = example_skincolor3) %>% 
  filter(ttT!="NULL") %>% 
  filter(treatment %>% str_detect("LightMed"))
# temp2$ttT[[1]] %>% select(1) %>% write_csv("./user_wx/degene_lightmed005.csv")
# pass the genelist to gorilla got GO_lightmed005.xlsx
GO_lm = read_excel("/home/xu/ses-1/user_wx/GO_lightmed005.xlsx")



#' ### intersection of GOResults biological pathways(FDR <0.05) for darkblack and lightmed
#+ echo=F, eval=T, warning=FALSE, message=FALSE
db = GO_db %>% 
  filter(`FDR q-value` < 0.05) %>% 
  select(1,2)
  
lm = GO_lm %>% 
  filter(`FDR q-value` < 0.05) %>% 
  select(1,2)
db %>% 
  inner_join(lm) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()
