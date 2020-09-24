#' ---
#' title: DE up and down pathway analysis
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


#' ## income and ses 4 up down pathways and unique up down pathways
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


DE_updown = result %>% 
  map(pluck("ttT")) %>% 
  map(~ filter(.x,adj.P.Val<=0.05) %>% mutate(updown = ifelse(logFC>0,"up","down")) %>% select(1,8))

DE = result %>% 
  map(pluck("ttT")) %>% 
  map(~ filter(.x,adj.P.Val<=0.05) %>% pull(gene))


#+ echo=F, eval=T, warning=FALSE, message=FALSE
sig = Reduce(union, signatures$outcome_set[table1])

DE_remove1k = map(DE, ~ setdiff(.x, sig))


# ses 4 and income unique genes after removing all genes in our table 1 disease signatures
ses4_unique = Reduce(setdiff, list(DE_remove1k$ses_sss_composite, DE_remove1k$income_hh_ff5, DE_remove1k$sss_5, DE_remove1k$edu_max, DE_remove1k$SEI_ff5))
income_unique = Reduce(setdiff, list(DE_remove1k$income_hh_ff5, DE_remove1k$ses_sss_composite, DE_remove1k$sss_5, DE_remove1k$edu_max, DE_remove1k$SEI_ff5))

# ses 4 up and down
ses4_up = DE_updown$ses_sss_composite %>% filter(updown =="up") %>% pull(gene) %>% intersect(DE_remove1k$ses_sss_composite)
ses4_down = DE_updown$ses_sss_composite %>% filter(updown =="down") %>% pull(gene) %>% intersect(DE_remove1k$ses_sss_composite)

# income up and down
income_up = DE_updown$income_hh_ff5 %>% filter(updown =="up") %>% pull(gene) %>% intersect(DE_remove1k$income_hh_ff5)
income_down = DE_updown$income_hh_ff5 %>% filter(updown =="down") %>% pull(gene) %>% intersect(DE_remove1k$income_hh_ff5)

# ses4 unique up and down
ses4_unique_up = DE_updown$ses_sss_composite %>% filter(updown =="up") %>% pull(gene) %>% intersect(ses4_unique)
ses4_unique_down = DE_updown$ses_sss_composite %>% filter(updown =="down") %>% pull(gene) %>% intersect(ses4_unique)


# income unique up and down
income_unique_up = DE_updown$income_hh_ff5 %>% filter(updown =="up") %>% pull(gene) %>% intersect(income_unique)
income_unique_down = DE_updown$income_hh_ff5 %>% filter(updown =="down") %>% pull(gene) %>% intersect(income_unique)

#  income ses 4 intersection up and down
income_ses4_up = intersect(income_up,ses4_up)
income_ses4_down = intersect(income_down,ses4_down)
# ses 4 and income intersection up and down

ses_income_intersection = intersect(DE_remove1k$ses_sss_composite, DE_remove1k$income_hh_ff5)
ses4_income_intersection_up = DE_updown$income_hh_ff5 %>% filter(updown =="up") %>% pull(gene) %>% intersect(ses_income_intersection)
ses4_income_intersection_down = DE_updown$income_hh_ff5 %>% filter(updown =="down") %>% pull(gene) %>% intersect(ses_income_intersection)



#+ echo=F, eval=T, warning=FALSE, message=FALSE

ses4_up = data.frame("hgnc_symbol" = ses4_up)
ses4_down = data.frame("hgnc_symbol" = ses4_down)
income_up = data.frame("hgnc_symbol" = income_up)
income_down = data.frame("hgnc_symbol" = income_down)

ses4_unique_up = data.frame("hgnc_symbol" = ses4_unique_up)
ses4_unique_down = data.frame("hgnc_symbol" = ses4_unique_down)
income_unique_up = data.frame("hgnc_symbol" = income_unique_up)
income_unique_down = data.frame("hgnc_symbol" = income_unique_down)


pathway = list(ses4_up = ses4_up,
               ses4_down = ses4_down,
               income_up = income_up,
               income_down =income_down,
               ses4_unique_up = ses4_unique_up,
               ses4_unique_down = ses4_unique_down,
               income_unique_up = income_unique_up,
               income_unique_down =income_unique_down) %>% 
  map(~ left_join(.x, entrezgeneid, by = "hgnc_symbol"))%>% 
  map(~ dplyr::pull(.x,entrezgene_id)) %>% 
  map(~ ReactomePA::enrichPathway(.x)) %>% 
  map(pluck("result")) %>% 
  map( ~ dplyr::filter(.x, p.adjust<=0.05)) %>% 
  map( ~ dplyr::select(.x, 2, 6))
pathway %>% print(row.names = FALSE)
# for(i in seq_along(pathway)) {
#   print(
#     kableExtra::kable(pathway[[i]], format = "html", caption = names(pathway)[i], longtable = TRUE, row.names = FALSE) %>%
#       kableExtra::kable_styling()
#   )
# }



#'`rmarkdown::render("/home/xu/ses-1/user_wx/DE_withoutsig.R")`
#'