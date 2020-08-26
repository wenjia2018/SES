# mediation for significant PCs

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

walk(dir(path = here("R"),full.names = TRUE), source)
fit_m4 = partial(fit_m4, n_perm = 1000) # specify n_perm


############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################

load_data(reconciled = FALSE, remove_inflam = FALSE)
define_treatments_and_controls()
recode_variables_in_dat()
print(abbreviations)
funcs = str_subset(abbreviations$shorthand, "^m") 
funcs = funcs %>% str_subset("m7")

# explicitly assign ncomp as the smallest number of table signatures gene numbers

ncomp = signatures$outcome_set[table1] %>% map_dfc(length) %>% unlist() %>%  min
ncomp = 10
fit_pca_util = partial(fit_pca_util, ncomp = ncomp) # specify n_perm

# debugonce(model_fit)
example0 =
  args %>%
  filter(gene_set_name =="Rheumatoid_Arthritis_mRNA",
         treatment == "income_hh_ff5",
         names(controls) == "all") %>% 
  mutate(out = pmap(., safely(model_fit), funcs),
         controls = names(controls))

# saveRDS(example0, "/home/share/scratch/xu/example0_4mediator.rds")


