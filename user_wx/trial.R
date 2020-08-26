set.seed(123)
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
walk(dir(path = "R",full.names = TRUE), source)
fit_m4 = partial(fit_m4, n_perm = 1000) # specify n_perm

############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################

load_data(reconciled = FALSE)
define_treatments_and_controls()
recode_variables_in_dat()
funcs = "m7_ob"
# debugonce(model_fit)
example3 =
args %>%
filter(gene_set_name == "Rheumatoid_Arthritis_mRNA") %>%
sample_n(2) %>%
mutate(out = pmap(., safely(model_fit), funcs),
controls = names(controls))

