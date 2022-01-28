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
# set some parameters values
# which PCA to perform
oblimin = FALSE
nn = TRUE
# explicitly assign ncomp for PCA analysis
ncomp = 10
fit_pca_util = partial(fit_pca_util, ncomp = ncomp) # specify n_perm
# type of mediation
mediation_mean = FALSE
mediation_each_gene = FALSE
# for doing genowide DE analysis only
normalization_bydesign = TRUE
# specify if subjects with disease shall be removed
remove_diseased_subjects = TRUE

# load data, specify remove_inflam (whether to remove inflammatory genes)
load_data(reconciled = FALSE, remove_inflam = TRUE)

# define your own treatment, controls, and re code variables to suit your needs

# for ses
define_treatments_and_controls()
recode_variables_in_dat()


# choose which function you want to run
funcs <- c("m1", "m2","m3", "m7", "m8", "m10")
funcs = "m12"
if(funcs == "m13"){
  boot = TRUE
  N = 1000
  
}else{
  boot = FALSE
}

# for parallel computing
# plan(multicore, workers = 40)
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

# save the results
example0 %>% saveRDS("./user_wx/m12_withoutinflam_ses_3setscontrols.rds")

