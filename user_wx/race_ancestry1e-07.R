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
# load custom ancestry PC
pval = "1e-07"
custom_PCA <- readRDS(str_c("/home/share/dna_ancestry/dna/custom_PCA_", pval,".rds"))
custom_PCA = custom_PCA %>% mutate(AID = AID %>% as.character())

load_data(reconciled = FALSE, remove_inflam = FALSE)
define_treatments_and_controls_race_dummy()
recode_variables_in_dat_racedummy()
print(abbreviations)
funcs = str_subset(abbreviations$shorthand, "^m") 
funcs = funcs %>% str_subset("m[7-8]")
# explicitly assign ncomp as the smallest number of table signatures gene numbers

ncomp = 10
fit_pca_util = partial(fit_pca_util, ncomp = ncomp) # specify n_perm
# debugonce(model_fit)


example0 =
  args %>%
  filter(is.element(gene_set_name, table1),
         names(controls) %in% c("basic", "ses", "ancestryPC")) %>% 
  mutate(out = pmap(., safely(model_fit), funcs),
         controls = names(controls))
# With controls used in SES paper: predicts the signatures from SES paper + Peters aging signature. 

# example0 %>% saveRDS("./user_wx/example_race_dummy_09.11.rds") with ancestry pca from addhealth

# with ancestry pca by us
example0 %>% saveRDS(str_c("./user_wx/example_race_dummy", pval,".rds"))

example1 =
  args %>%
  filter(names(controls) %in% c("basic", "ses", "ancestryPC"),
         gene_set_name == "whole_genome_and_tfbm") %>%
  mutate(out = pmap(., safely(model_fit), funcs),
         controls = names(controls))
# example1 %>% saveRDS("./user_wx/race_detfbm_dummy_06.11.rds") with ancestry pca from addhealth

# with ancestry pca by us
example1 %>% saveRDS(str_c("./user_wx/race_detfbm_dummy", pval,".rds"))
