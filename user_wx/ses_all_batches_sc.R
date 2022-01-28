# main function to start running analysis
# for skin color with addhealth ancestry controls
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
normalization_bydesign = FALSE
# specify if subjects with disease shall be removed
remove_diseased_subjects = TRUE

# load data, specify remove_inflam (whether to remove inflammatory genes)
load_data(reconciled = FALSE, remove_inflam = TRUE)

# define your own treatment, controls, and re code variables to suit your needs
define_treatments_and_controls_sc()
recode_variables_in_dat_racedummy()

# choose which function you want to run
funcs <- c("m1", "m2","m3", "m7", "m8", "m10")
funcs = "m12"
if(funcs == "m13"){
  boot = TRUE
  N = 1000
  
}else{
  boot = FALSE
}


plan(multicore, workers = 36)
# debugonce(model_fit)
  example0 =
  args %>%
  filter(
    is.element(gene_set_name, table1) &
      # gene_set_name=="aging_mRNA"&
    # gene_set_name %in% c("whole_genome_and_tfbm"),
    # treatment %in% c("edu_max"),
         (
           # names(controls) == "basic" | 
             names(controls) == "all")) %>%
  mutate(out = furrr::future_pmap(., safely(model_fit), funcs),
         controls = names(controls))

# With controls used in SES paper: predicts the signatures from SES paper + Peters aging signature. 
example0 %>% saveRDS("./user_wx/m12_with1k_aging_sccont_NonBstrata_correctioninunionxxx.rds")

