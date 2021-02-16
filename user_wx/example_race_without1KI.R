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
define_treatments_and_controls_race()
recode_variables_in_dat()
print(abbreviations)
funcs = str_subset(abbreviations$shorthand, "^m") 
funcs = funcs %>% str_subset("m[7-8]")
# explicitly assign ncomp as the smallest number of table signatures gene numbers

ncomp = 10
fit_pca_util = partial(fit_pca_util, ncomp = ncomp) # specify n_perm
# debugonce(model_fit)
# debugonce(model_MR)
example0 =
  args %>%
  filter(is.element(gene_set_name, table1),
         treatment %in% c("raceethnicity", "raceethnicity2", "raceethnicity3"),
         names(controls) == "all") %>% 
  mutate(out = pmap(., safely(model_fit), funcs),
         controls = names(controls))
# With controls used in SES paper: predicts the signatures from SES paper + Peters aging signature. 
example0 %>% saveRDS("./user_wx/example_race_without1KI_1027.rds") 


# example1 =
#   args %>%
#   filter(treatment =="raceethnicity",
#          gene_set_name == "whole_genome_and_tfbm",
#          names(controls) == "all") %>% 
#   mutate(out = pmap(., safely(model_fit), funcs),
#          controls = names(controls))
