
#' ---
#' title: Nice code
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE


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

walk(dir(path = here("R"),full.names = TRUE), source)
fit_m4 = partial(fit_m4, n_perm = 1000) # specify n_perm

# WHICH EXAMPLES TO RUN? 
example4 <- example3 <- example2 <- example1 <- example0 <- FALSE
example4 <- TRUE 

############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################

load_data(reconciled = FALSE)
define_treatments_and_controls()
recode_variables_in_dat()
print(abbreviations)
funcs = str_subset(abbreviations$shorthand, "^m") 
funcs = funcs %>% setdiff(c("m4","m5b","m99"))

############################################################
# EXAMPLE: SIGNATURES
############################################################

if(example0){
  
  if(from_disk <- TRUE){
    
    example0 = readRDS("/home/share/scratch/example0.rds")
    
  } else {
    
    example0 =
      args %>%
      filter(is.element(gene_set_name, table1),
             names(controls) == "all") %>% 
      mutate(out = pmap(., safely(model_fit), funcs),
             controls = names(controls))
    
    saveRDS(example0, "/home/share/scratch/example0.rds")
    
  }
  
  example0 = remove_errors(example0) 
  get_table1(example0)
  
  # ESTIMATE VARIOUS PCA "ROTATIONS"
  example0_m7_nn = example0 %>% get_sig_PCs_and_sig_enrichment_on_those_PCs("m7_nn")  
  example0_m7_vx = example0 %>% get_sig_PCs_and_sig_enrichment_on_those_PCs("m7_vx")
  example0_m7_ob = example0 %>% get_sig_PCs_and_sig_enrichment_on_those_PCs("m7_ob")
  
  # PICK A ROTATION
  which_rotation = example0_m7_nn
  # remove unused packages (including dependency) in this session which has conflicts with tidyverse
  list("clusterProfiler", "DO.db", "ReactomePA", "reactome.db", "DOSE", "graphite", "enrichplot",  "GO.db", "GOSemSim" ,
       "org.Hs.eg.db", "AnnotationDbi", "IRanges", "S4Vectors") %>% map(detach_package, TRUE)
  
  # INSPECT MODELS WHICH HAVE SIGNIFICANT PCs
  interesting_PCS =
    which_rotation %>%
    dplyr::select(treatment, controls, gene_set_name, matches("well_loaded")) %>%
    drop_na()
  print(interesting_PCS, n = Inf)
  
  # PICK YOUR FAVORITE ROW OF PRECEDING TABLE TO VISUALIZE
  particularly_interesting_row = 1
  
  # WHAT ARE THE WELL-LOADED GENES FOR EACH SIGNIFICANT PC IN THIS ROW?
  interesting_PCS %>% dplyr::slice(particularly_interesting_row) %>% pluck(4)
  
  # WHAT IS THE ENRICHMENT STATUS OF THESE PCs
  interesting_PCS %>% dplyr::slice(particularly_interesting_row) %>% pluck(5)
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
           names(controls) == "all") %>%
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

############################################################
# EXAMPLE: cibersort compositional 
############################################################

if(example3){ 
  # note gene_set_name is not used in this example, but in order to be consistent
  # with the model_fit arguments, any gene_set_name needs to be filled, also controls
  # has to be fixed as basic for cell type compositional analysis
  # (no need to control cell type)
  # debugonce(model_fit)
  # example3 =
  #   args %>% 
  #   filter(treatment == "ses_sss_composite",
  #          names(controls) == "basic",
  #          gene_set_name == "Rheumatoid_Arthritis_mRNA") %>%
  #   mutate(out = pmap(., safely(model_fit), "m96"),
  #          controls = names(controls))
  
  example3 =
    args %>% 
    filter(gene_set_name == "Rheumatoid_Arthritis_mRNA",
           names(controls) == "basic") %>%
    mutate(out = pmap(., safely(model_fit), "m96"),
           controls = names(controls))
}

############################################################
# EXAMPLE: PCA component mediational analysis
############################################################
if(example4){
  
  if(from_disk <- TRUE){
    
    example4 = readRDS("/home/share/scratch/example4.rds")
    
  } else { 
    
   # rename please: lets keep try to keep most variables immutable to reduce unintended side effects
    
    # debugonce(model_fit)
    example4 =
      args %>% 
      filter(treatment == "ses_sss_composite",
             names(controls) == "all",
             gene_set_name == "Rheumatoid_Arthritis_mRNA") %>%
      mutate(out = pmap(., safely(model_fit), "m7"),
             controls = names(controls))
    
    example4 %>% saveRDS("/home/share/scratch/example4.rds")
    
  }
  
  example4 %>%
    hoist(out, "error") %>%
    mutate(error = map(error, as.character)) %>%
    unnest(error)
  
  example4 %>%
    hoist(out, mediation = list("result", "m7_nn", 1, "mediation")) %>% 
    unnest_wider("mediation") %>% 
    hoist(phys_activ_ff5, "result") %>% 
    unnest_wider("result") %>% 
    unnest(matches("^d")) %>%
    unnest(matches("^d"))
  
}


rmarkdown::render("/home/xu/ses-1/user_wx/output.R")
