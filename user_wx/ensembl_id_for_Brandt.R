
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
define_treatments_and_controls_skincolor_fullsibFE()
recode_variables_in_dat_racedummy()
print(abbreviations)
funcs = str_subset(abbreviations$shorthand, "^m") 
funcs = funcs %>% str_subset("m[7-8]")

G_list = readRDS("/home/share/data_input/genename_15062020.rds") 

a = G_list %>% filter(hgnc_symbol %in% signatures$outcome_set$aging_mRNA) %>% dplyr::pull(ensembl_gene_id)

a %>% openxlsx::write.xlsx("./user_wx/aging_list.xlsx")
