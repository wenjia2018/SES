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
library(dbr) # my package
library(MendelianRandomization)
walk(dir(path = here("R"),full.names = TRUE), source)
fit_m4 = partial(fit_m4, n_perm = 1000) # specify n_perm


############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################

load_data(reconciled = FALSE, remove_inflam = FALSE)
define_treatments_and_controls_skincolor_dummy()
recode_variables_in_dat_racedummy()
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
  filter(is.element(gene_set_name, table1)) %>% 
  mutate(out = pmap(., safely(model_fit), funcs),
         controls = names(controls))
# With controls used in SES paper: predicts the signatures from SES paper + Peters aging signature. 
example0 %>% saveRDS("./user_wx/example_skincolor_dummy_13.11.rds") 
# example_skincolor_dummy_12.11.rds basic control case, does not exclude asian and other nonhispanic, as
# previously use race ethnicity to exclude those people automatically

# # 
# example1 =
#   args %>%
#   filter(gene_set_name == "whole_genome_and_tfbm") %>%
#   mutate(out = pmap(., safely(model_fit), funcs),
#          controls = names(controls))
# example1 %>% saveRDS("./user_wx/skincolor_detfbm_dummy_13.11.rds")
# # 
