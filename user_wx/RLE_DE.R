
library(here)

library(tidyverse)
library(EValue)
library(rlang)
library(skimr)
library(furrr)
library(limma)
# library(recipes)
# recipes and Evalue has conflict, after loading recipes, evalue package doesnot work with the error Error: $ operator is invalid for atomic vectors
library(parsnip)
library(workflows)
library(Biobase)
library(dbr) # my package

walk(dir(path = here("R"),full.names = TRUE), source)
fit_m4 = partial(fit_m4, n_perm = 1000) # specify n_perm
source("./dda_v0.1/dda_resdist.r")
source("./dda_v0.1/dda_vardist.r")
source("./dda_v0.1/dda_indep.r")
source("./dda_v0.1/boot_hsic_test.R")
source("./dda_v0.1/nlcor_test.r")
n_boot = 5000


############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################
log2cpm = FALSE
normalization_bydesign = TRUE
load_data(reconciled = TRUE, remove_inflam = FALSE)
define_treatments_and_controls()
recode_variables_in_dat()
print(abbreviations)
mediation_mean = FALSE
mediation_each_gene = FALSE
funcs = str_subset(abbreviations$shorthand, "^m") 
funcs <- c("m1", "m2","m3", "m7", "m8", "m10")
funcs = "m12"
# explicitly assign ncomp as the smallest number of table signatures gene numbers

ncomp = 10
fit_pca_util = partial(fit_pca_util, ncomp = ncomp) # specify n_perm

example0 =
  args %>%
  filter(
    is.element(gene_set_name, table1) &
      # gene_set_name== "inflam1k_mRNA" &
      #   treatment == "ses_sss_composite" &
      names(controls) == "basic") %>%
  mutate(out = pmap(., safely(model_fit), funcs),
         controls = names(controls))

# With controls used in SES paper: predicts the signatures from SES paper + Peters aging signature. 
example0 %>% saveRDS("./user_wx/example_RLE_DE_m12_nonparalle.rds") 
