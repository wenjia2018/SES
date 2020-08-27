#' ---
#' title: DE results for ses
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
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase)
library(ggformula)
library(here)
library(dbr)
library(venn)
walk(dir(path = here("R"),full.names = TRUE), source)
############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################

load_data(reconciled = FALSE, remove_inflam = FALSE)
define_treatments_and_controls()
recode_variables_in_dat()
# print(abbreviations)
funcs = str_subset(abbreviations$shorthand, "^m")
funcs = funcs %>% str_subset("m[6-8]")
# explicitly assign ncomp as the smallest number of table signatures gene numbers

ncomp = signatures$outcome_set[table1]%>% map_dfc(length) %>% unlist() %>%  min
if(0){
  example1 =
  args %>%
    filter(treatment %in%c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5", "sss_5",
                           # "ses_composite_pp1", "edu_p", "SEI_max_p_w12", "income_pp1_log"
                           ),
           gene_set_name == "whole_genome_and_tfbm",
           names(controls) == "all") %>% slice_(2) %>% 
    mutate(out = pmap(., safely(model_fit), funcs),
           controls = names(controls))
}
example1 = example1 %>% filter(treatment %in%c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5"))

#' ## DE results for ses4, edu, income, SEI
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example1 = readRDS("/home/xu/ses-1/user_wx/DE_with1k.rds")

# load our whole gene name list
entrezgeneid =  readRDS("/home/share/preprocessed_two_batches/entrezgeneid.rds")

# ERRORS?
error = example1 %>%
  hoist(out, "error") %>%
  mutate(error = map(error, as.character)) %>%
  unnest(error)

# RESULTS
temp = example1 %>%
  hoist(out, "result")


result = temp %>%
  pluck("result") %>%
  set_names(temp$treatment)

#' #### ses DE genes
#' 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
DE_updown = result %>% 
  map(pluck("ttT")) %>% 
  map(~ filter(.x,adj.P.Val<=0.05) %>% mutate(updown = ifelse(logFC>0,"up","down")) %>% select(1,8))
      
DE = result %>% 
  map(pluck("ttT")) %>% 
  map(~ filter(.x,adj.P.Val<=0.05) %>% pull(gene))

DE %>% print(row.names = FALSE)

#' #### venn diagram for ses DE genes
#' 
#+ echo=F, eval=T, warning=FALSE, message=FALSE

venn::venn(DE, ilabels = TRUE, zcolor = "style", ellipse = FALSE, opacity = 0.15)

#' #### pathway for each ses DE genes
#' 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
pathway = result %>%
  map(pluck("fig1")) %>%
  map(pluck("reactome_pathway")) %>%
  map(pluck("result")) %>%
  map( ~ dplyr::select(.x, 2 ,6) %>% filter(p.adjust<=0.05))
pathway %>% print(row.names = FALSE)

#' #### pathway for each ses unique DE genes
#' income and ses4 has unique genes, sss only have 1 unique gene "ZCCHC10"
#+ echo=F, eval=T, warning=FALSE, message=FALSE

ses4_unique = Reduce(setdiff, list(DE$ses_sss_composite, DE$income_hh_ff5, DE$sss_5, DE$edu_max, DE$SEI_ff5))
income_unique = Reduce(setdiff, list(DE$income_hh_ff5, DE$ses_sss_composite, DE$sss_5, DE$edu_max, DE$SEI_ff5))

#' ses4 unique DE genes

print(ses4_unique)

#' income unique DE genes
#' 
print(income_unique)

#+ echo=F, eval=T, warning=FALSE, message=FALSE

ses4_unique = data.frame("hgnc_symbol" = ses4_unique)
income_unique = data.frame("hgnc_symbol" = income_unique)

pathway = list(ses4 = ses4_unique, income = income_unique) %>% 
  map(~ left_join(.x, entrezgeneid, by = "hgnc_symbol"))%>% 
  map(~ dplyr::pull(.x,entrezgene_id)) %>% 
  map(~ ReactomePA::enrichPathway(.x)) %>% 
  map(pluck("result")) %>% 
  map( ~ dplyr::filter(.x, p.adjust<=0.05)) %>% 
  map( ~ dplyr::select(.x, 2, 6))

pathway %>% print(row.names = FALSE)


#' ## DE results for ses4, edu, income, SEI, after removing all signatures in table 1
#+ echo=F, eval=T, warning=FALSE, message=FALSE


#' #### ses DE genes

#+ echo=F, eval=T, warning=FALSE, message=FALSE
sig = Reduce(union, signatures$outcome_set[table1])

DE_remove1k = map(DE, ~ setdiff(.x, sig))

DE_remove1k %>% print(row.names = FALSE)

#' #### venn diagram for ses DE genes

#+ echo=F, eval=T, warning=FALSE, message=FALSE

venn::venn(DE_remove1k, ilabels = TRUE, zcolor = "style", ellipse = FALSE, opacity = 0.15)

# x = VennDiagram::venn.diagram(DE_remove1k, ilabels = TRUE, zcolor = "style", ellipse = FALSE, opacity = 0.15, filename = NULL)

#' #### pathway for each ses DE genes
#'
#+ echo=F, eval=T, warning=FALSE, message=FALSE


pathway = DE_remove1k %>%
  map(~ data.frame("hgnc_symbol" =.x)) %>% 
  map(~ left_join(.x, entrezgeneid, by = "hgnc_symbol"))%>% 
  map(~ dplyr::pull(.x, entrezgene_id)) %>% 
  map(~ ReactomePA::enrichPathway(.x)) %>% 
  map(pluck("result")) %>% 
  map( ~ dplyr::filter(.x, p.adjust<=0.05)) %>% 
  map( ~ dplyr::select(.x, 2, 6))

pathway %>% print(row.names = FALSE)

#' #### pathway for each ses unique DE genes
#' income and ses4 has unique genes, sss only have 1 unique gene "ZCCHC10"
#+ echo=F, eval=T, warning=FALSE, message=FALSE

# ses 4 and income unique genes after removing all genes in our table 1 disease signatures
ses4_unique = Reduce(setdiff, list(DE_remove1k$ses_sss_composite, DE_remove1k$income_hh_ff5, DE_remove1k$sss_5, DE_remove1k$edu_max, DE_remove1k$SEI_ff5))
income_unique = Reduce(setdiff, list(DE_remove1k$income_hh_ff5, DE_remove1k$ses_sss_composite, DE_remove1k$sss_5, DE_remove1k$edu_max, DE_remove1k$SEI_ff5))

# ses 4 up and down
ses4_up = DE_updown$ses_sss_composite %>% filter(updown =="up") %>% pull(gene) %>% intersect(DE_remove1k$ses_sss_composite)
ses4_down = DE_updown$ses_sss_composite %>% filter(updown =="down") %>% pull(gene) %>% intersect(DE_remove1k$ses_sss_composite)

# income unique up and down
income_unique_up = DE_updown$income_hh_ff5 %>% filter(updown =="up") %>% pull(gene) %>% intersect(income_unique)
income_unique_down = DE_updown$income_hh_ff5 %>% filter(updown =="down") %>% pull(gene) %>% intersect(income_unique)

# ses 4 and income intersection up and down

ses_income_intersection = intersect(DE_remove1k$ses_sss_composite, DE_remove1k$income_hh_ff5)
ses4_income_intersection_up = DE_updown$income_hh_ff5 %>% filter(updown =="up") %>% pull(gene) %>% intersect(ses_income_intersection)
ses4_income_intersection_down = DE_updown$income_hh_ff5 %>% filter(updown =="down") %>% pull(gene) %>% intersect(ses_income_intersection)

list(ses4_up = ses4_up,
     ses4_down = ses4_down,
     income_unique_up = income_unique_up,
     income_unique_down = income_unique_down,
     ses4_income_intersection_up = ses4_income_intersection_up,
     ses4_income_intersection_down = ses4_income_intersection_down
       ) %>% openxlsx::write.xlsx("./user_wx/DE_removetable1signatures_ses4income.xlsx")

# sss

DE_remove1k$sss_5
Reduce(intersect, list(DE_remove1k$ses_sss_composite, DE_remove1k$income_hh_ff5, DE_remove1k$sss_5, DE_remove1k$edu_max, DE_remove1k$SEI_ff5))
Reduce(intersect, list(DE_remove1k$ses_sss_composite, DE_remove1k$income_hh_ff5, DE_remove1k$sss_5)) %>%
  setdiff(DE_remove1k$edu_max) %>%
  setdiff(DE_remove1k$SEI_ff5)
Reduce(intersect, list(DE_remove1k$ses_sss_composite, DE_remove1k$income_hh_ff5, DE_remove1k$sss_5,DE_remove1k$edu_max))%>%
  setdiff(DE_remove1k$SEI_ff5)
#' ses4 unique DE genes

print(ses4_unique)

#' income unique DE genes
#' 
print(income_unique)

#+ echo=F, eval=T, warning=FALSE, message=FALSE

ses4_unique = data.frame("hgnc_symbol" = ses4_unique)
income_unique = data.frame("hgnc_symbol" = income_unique)

pathway = list(ses4 = ses4_unique, income = income_unique) %>% 
  map(~ left_join(.x, entrezgeneid, by = "hgnc_symbol"))%>% 
  map(~ dplyr::pull(.x,entrezgene_id)) %>% 
  map(~ ReactomePA::enrichPathway(.x)) %>% 
  map(pluck("result")) %>% 
  map( ~ dplyr::filter(.x, p.adjust<=0.05)) %>% 
  map( ~ dplyr::select(.x, 2, 6))

pathway %>% print(row.names = FALSE)


#'`rmarkdown::render("/home/xu/ses-1/user_wx/DE_withandwithout1K.R")`
