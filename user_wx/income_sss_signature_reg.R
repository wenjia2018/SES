
#' ---
#' title: regression with income and SSS
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE, echo = FALSE



#' We do a regression with income and SSS in one model controlling for full covariates
#' to see which is significant, controlling the other


#+ echo=F, eval=T
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
library(dbr) # my package

# #+ new_chunk
# walk(dir(path = here("R"),full.names = TRUE), source)
# fit_m4 = partial(fit_m4, n_perm = 1000) # specify n_perm
# 
# 
# ############################################################
# # LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
# ############################################################
# 
# load_data(reconciled = FALSE)
# define_treatments_and_controls()
# recode_variables_in_dat()
# print(abbreviations)
# funcs = str_subset(abbreviations$shorthand, "^m")
# funcs = funcs %>% str_subset("m8")

############################################################
# EXAMPLE: SIGNATURES
############################################################
#+ echo=F, eval=T
if(0) {
  example0 =
    args %>%
    filter(is.element(gene_set_name, table1),
           # treatment =="sss_5",
           names(controls) == "all") %>% 
    mutate(out = pmap(., safely(model_fit), funcs),
           controls = names(controls))
  
  saveRDS(example0, "/home/share/scratch/xu/example0_m8fdr_withdetail.rds")
}
   
#' ## regression on income after controlling for sss   
#+ echo=F, eval=T
    income = readRDS("/home/share/scratch/example0_income.rds")
    
    
    income %>%  hoist(out, detail = list("result", "m8_fdr", 1, "detail")) %>% 
      filter(treatment =="income_hh_ff5") %>% 
      unnest(detail) %>% 
      dplyr::select(treatment, gene_set_name, gene, logFC, t, adj.P.Val) %>%
      kableExtra::kable() %>%
      kableExtra::kable_styling()
#' \newpage
#'     
#' ## regression on sss after controlling for income
#+ echo=F, eval=T
    sss = readRDS("/home/share/scratch/example0_sss.rds")
  
    sss %>% hoist(out, detail = list("result", "m8_fdr", 1, "detail")) %>% 
      filter(treatment =="sss_5") %>% 
      unnest(detail) %>% 
      dplyr::select(treatment, gene_set_name, gene, logFC, t, adj.P.Val) %>%
      kableExtra::kable() %>%
      kableExtra::kable_styling()
#' \newpage
#'     
#' ## regression on income alone
#+ echo=F, eval=T
    example0 = readRDS("/home/share/scratch/xu/example0_m8fdr_withdetail.rds")

    example0 %>% hoist(out, detail = list("result", "m8_fdr", 1, "detail")) %>%
      filter(treatment =="income_hh_ff5") %>%
      unnest(detail) %>%
      dplyr::select(treatment, gene_set_name, gene, logFC, t, adj.P.Val) %>%
      kableExtra::kable() %>%
      kableExtra::kable_styling()

#' \newpage
#' 
#' ## regression on sss alone
#+ echo=F, eval=T
    example0 %>% hoist(out, detail = list("result", "m8_fdr", 1, "detail")) %>%
      filter(treatment =="sss_5") %>%
      unnest(detail) %>%
      dplyr::select(treatment, gene_set_name, gene, logFC, t, adj.P.Val) %>%
      kableExtra::kable() %>%
      kableExtra::kable_styling()

    
#'  `rmarkdown::render("/home/xu/ses-1/user_wx/income_sss_signature_reg.R")`