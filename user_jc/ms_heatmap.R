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
print(abbreviations)
funcs = str_subset(abbreviations$shorthand, "^m") %>% setdiff(c("m4", "m97"))  

############################################################
# LOAD RESULTS
############################################################

example0 = 
  readRDS("/home/share/scratch/example0.rds") %>% 
  hoist(out, "result") %>% 
  drop_na()

############################################################
# TABLE 1
############################################################

table1 = get_table1(example0) 
table2 = table1 %>% 
  filter(controls == "all") %>%
  select(treatment,  gene_set_name, m8_fdr_p) %>%
  pivot_wider(names_from = treatment, values_from = m8_fdr_p)
table2
