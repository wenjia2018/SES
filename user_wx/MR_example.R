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
library(MendelianRandomization)
walk(dir(path = here("R"),full.names = TRUE), source)
fit_m4 = partial(fit_m4, n_perm = 1000) # specify n_perm


############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################

load_data(reconciled = FALSE, remove_inflam = TRUE)
define_treatments_and_controls()
recode_variables_in_dat()
print(abbreviations)
funcs = str_subset(abbreviations$shorthand, "^m") 
funcs = "m95"
# explicitly assign ncomp as the smallest number of table signatures gene numbers

ncomp = 10
fit_pca_util = partial(fit_pca_util, ncomp = ncomp) # specify n_perm
# debugonce(model_fit)
# debugonce(model_MR)
example0 =
  args %>%
  filter(is.element(gene_set_name, table1),
         treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5", "sss_5"), 
         names(controls) == "all") %>% 
  mutate(out = pmap(., safely(model_fit), funcs),
         controls = names(controls))

example0 %>% saveRDS("./user_wx/example_MR_table1_completecontrols.rds")

example_MR_table1_completecontrols %>%
  hoist(out, "result") %>%
  unnest_longer(result) %>%
  filter(gene_set_name!="inflam1k_mRNA") %>%
  mutate(result = result %>% unlist(use.names=FALSE)) %>%
  hoist(result, p_unadj = "Pvalue") %>%
  hoist(result, estimate = "Estimate") %>%
  mutate(controls = str_c(controls, "+", treatment),
         treatment = "w5bmi",
         IV = "PGSBMI") %>%
  dplyr::select(1:5, 9) %>%
  select(IV, everything()) %>%
  filter(p_unadj < 0.05)
example_MR_table1_nocontrols %>%
  mutate(treatment = "w5bmi") %>%
  hoist(out, "result") %>%
  unnest_longer(result) %>%
  filter(gene_set_name!="inflam1k_mRNA") %>%
  mutate(result = result %>% unlist(use.names=FALSE)) %>%
  hoist(result, p_unadj = "Pvalue") %>%
  hoist(result, estimate = "Estimate") %>%
  mutate(controls = "NULL",
         IV = "PGSBMI") %>%
  filter(p_unadj <0.05) %>%
  dplyr::select(1:5, 9) %>%
  select(IV, everything()) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()
