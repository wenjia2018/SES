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

if(example0){
  
  if(from_disk <- TRUE){
    
    example0 = readRDS("/home/share/scratch/example0.rds")
    
  } else {
    
    example0 =
      args %>%
      filter(is.element(gene_set_name, table1)) %>%
      mutate(out = pmap(., safely(model_fit), funcs),
             controls = names(controls))
    
    saveRDS(example0, "/home/share/scratch/example0.rds")
    
  }
  
  ############################################################
  # ANY ESTIMATION ERRORS?
  ############################################################
  
  # ERRORS:
  example0 %>%
    hoist(out, "error") %>%
    mutate(error = map(error, as.character)) %>%
    unnest(error) %>%
    group_by(error) %>%
    slice(1)
  
  # WHAT CAUSES ERROR? RELATE NA TO ARGS OF model_fit()
  example0 %>%
    hoist(out, p = list("result", "m1", 1, "p")) %>%
    with(table(gene_set_name, is.na(p)))
  
  # REMOVE MODELS THAT ERR
  example0 = example0 %>% hoist(out, "result") %>% drop_na()
  
  ############################################################
  # GET TABLE 1
  ############################################################
  
  get_table1(example0) %>% print(n = Inf)
  
  ############################################################
  # HARVEST SIGNIFICANT PCs and THEIR POSSIBLY SIGNIFICANT ENRICHMENT
  ############################################################
  
  # FIRST PICK A PCA "ROTATION"
  m7_model = "m7_nn" # of "m7_nn", "m7_vx", "m7_ob"
  example0 = example0 %>% get_sig_PCs_and_sig_enrichment_on_those_PCs(m7_model)
  
  # INSPECT MODELS WHICH HAVE SIGNIFICANT PCs
  interesting_PCS =
    example0 %>%
    select(treatment, controls, gene_set_name, matches("well_loaded")) %>%
    drop_na()
  print(interesting_PCS, n = Inf)
  
  # PICK YOUR FAVORITE ROW OF PRECEDING TABLE TO VISUALIZE
  particularly_interesting_row = 1
  
  # WHAT ARE THE WELL-LOADED GENES FOR EACH SIGNIFICANT PC IN THIS ROW?
  interesting_PCS %>% slice(particularly_interesting_row) %>% pluck(4)
  
  # WHAT IS THE ENRICHMENT STATUS OF THESE PCs
  interesting_PCS %>% slice(particularly_interesting_row) %>% pluck(5)
}

############################################################
# EXAMPLE: WHOLE GENOME
############################################################

if(example1){
  
  
  # FILTER OUT YOUR FAVORITE LHS-RHS AND LOOK AT TFBM, DE AND ENRICHMENT PLOTS
  example1 =
    args %>%
    filter(treatment == "ses_sss_composite",
           gene_set_name == "whole_genome_and_tfbm",
           names(controls) == "basic") %>%
    mutate(out = pmap(., safely(model_fit), funcs),
           controls = names(controls))
  
  # ERRORS?
  example1 %>%
    hoist(out, "error") %>%
    mutate(error = map(error, as.character)) %>%
    unnest(error)
  
  # RESULTS
  example1 %>%
    hoist(out, "result") %>%
    pluck("result")
  
}
