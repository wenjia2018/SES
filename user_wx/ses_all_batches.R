# main function to start running analysis
# load packages
library(here)
# library(data.table)
# setDTthreads(threads = 20)
library(tidyverse)
library(EValue)
library(rlang)
library(skimr)
library(furrr)
library(limma)
# library(recipes)
# recipes and Evalue has conflict, after loading recipes,
# evalue package doesnot work with the error Error: $ operator is invalid for atomic vectors
library(parsnip)
library(workflows)
library(Biobase)
library(dbr) # my package

# source all the useful files
walk(dir(path = here("R"),full.names = TRUE), source)

############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################
# set some parmeter values
# which PCA to perform
oblimin = FALSE
nn = TRUE
# explicitly assign ncomp for PCA analysiss
ncomp = 10
# for doing genowide DE analysis only
normalization_bydesign = FALSE
# specify if subjects with disease shall be removed
remove_diseased_subjects = TRUE
load_data(reconciled = FALSE, remove_inflam = TRUE)
# for skin color
# define_treatments_and_controls_sc()
define_treatments_and_controls()
recode_variables_in_dat()

mediation_mean = FALSE
mediation_each_gene = FALSE
# choose which function you want to run
funcs <- c("m1", "m2","m3", "m7", "m8", "m10")
funcs = "m12"
if(funcs == "m13"){
  boot = TRUE
  N = 1000
  
}else{
  boot = FALSE
}
fit_pca_util = partial(fit_pca_util, ncomp = ncomp) # specify n_perm

plan(multicore, workers = 40)
# debugonce(model_fit)
  example0 =
  args %>%
  filter(
    is.element(gene_set_name, table1) &
    # gene_set_name %in% c( "Hypertension_mRNA", "Rheumatoid_Arthritis_mRNA"),
    # treatment %in% c("edu_max"),
         names(controls) %in% c("basic","basic_less","basic_less_ancestry") ) %>%
  mutate(out = furrr::future_pmap(., safely(model_fit), funcs),
         controls = names(controls))

# With controls used in SES paper: predicts the signatures from SES paper + Peters aging signature. 
example0 %>% saveRDS("./user_wx/m12_withoutinflam_ses_3setscontrols.rds")

