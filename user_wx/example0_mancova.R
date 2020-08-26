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

# WHICH EXAMPLES TO RUN?
example1 <- example0 <- TRUE

############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################

load_data(reconciled = FALSE)
define_treatments_and_controls()
recode_variables_in_dat()
print(abbreviations)
funcs = "m1"

############################################################
# EXAMPLE: SIGNATURES
############################################################

  
  if(from_disk <- FALSE){
    
    example0 = readRDS("/home/xu/ses-1/user_wx/example0.rds")
    
  } else {
    
    example0 =
      args %>%
      filter(is.element(gene_set_name, table1)) %>%
      mutate(out = pmap(., safely(model_fit), funcs),
             controls = names(controls))
    
    saveRDS(example0, "/home/xu/ses-1/user_wx/example0.rds")
    
  }