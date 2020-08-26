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

load_data(reconciled = FALSE)
define_treatments_and_controls()
recode_variables_in_dat()
# print(abbreviations)
# funcs = str_subset(abbreviations$shorthand, "^m") 
# funcs = c("m96")


example3 =
  args %>% 
  filter(gene_set_name == "Rheumatoid_Arthritis_mRNA",
         names(controls) == "basic") %>%
  mutate(out = pmap(., safely(model_fit), "m96"),
         controls = names(controls))

saveRDS(example3, "/home/share/scratch/xu/example3_celltype.rds")

