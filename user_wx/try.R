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
funcs = str_subset(abbreviations$shorthand, "^m") %>% setdiff(c("m4", "m97"))  

############################################################
# EXAMPLE: SIGNATURES
############################################################


    example0 =
      args %>%
      filter(is.element(gene_set_name, table1)) %>%
      mutate(out = pmap(., safely(model_fit), funcs),
             controls = names(controls))
    
    saveRDS(example0, "/home/share/scratch/xu/example0.rds")
    
    
    
    example1 =
      args %>%
      filter(treatment %in% c("ses_composite_pp1","ses_sss_composite","income_hh_ff5","edu_max","SEI_ff5","sss_5"),
             gene_set_name == "whole_genome_and_tfbm",
             names(controls) == "all") %>%
      mutate(out = pmap(., safely(model_fit), funcs),
             controls = names(controls))

    
    saveRDS(example1, "/home/share/scratch/xu/example1.rds")
    