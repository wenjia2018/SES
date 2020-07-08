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

example1 <- example0 <- TRUE

############################################################
# LOAD DATA
############################################################

load_data(reconciled = FALSE) %>% 
  list2env(.GlobalEnv)

if(discritize_exposures <- TRUE) recode_variables_in_dat()
source("user_ms/define_treatments_controls_outcomes.R")
print(abbreviations)
print("Select which models to estimate from the above table.")
funcs = str_subset(abbreviations$shorthand, "^m") %>% setdiff(c("m4", "m98")) # m98 breaks for some reason 
 
############################################################
# EXAMPLES
############################################################

if(example0){ 
  
  ############################################################
  # SIGNATURES
  ############################################################
  
  # fit_pca_util %>% debugonce()
  # model_fit %>% debugonce()
  example0 = 
    args %>% 
    filter(is.element(gene_set_name, table1)) %>% 
    mutate(out = pmap(., safely(model_fit), funcs), 
           controls = names(controls))
  
  saveRDS(example0, "rds/example0.rds")
  # example0 = readRDS("rds/example0.rds")
  
  ############################################################
  # ERRORS?
  ############################################################

  # group errors
  example0 %>% 
    hoist(out, "error") %>% 
    mutate(error = map_chr(error, as.character)) %>%
    unnest(error) %>%
    group_by(error)
  
  # where is the problem? relate NA to each of the arguments to model_fit(),
  # controls, treatment, gene_set_name:
  example0 %>% 
    hoist(out, p = list("result", "m1", 1, "p")) %>% 
    with(table(gene_set_name, is.na(p))) 
  
  ############################################################
  # unpack PCA
  ############################################################
  (
    tabPCA =
      example0 %>% 
      hoist(out, p = list("result" )) %>%
      unnest(p) %>% 
      unnest(matches("m6|m7")) %>%
      filter(names(m6_vx) == "p") %>%
      unnest(matches("m6")) %>%
      unnest_wider(m7_nn, names_sep = "_") %>% 
      unnest_wider(m7_vx, names_sep = "_") %>% 
      unnest_wider(m7_ob, names_sep = "_")
  )
  
  ############################################################
  # Unpack other
  ############################################################
  
  tab1a = 
    example0 %>% 
    hoist(out, "result") %>% 
    hoist(result, !!!funcs)  %>% 
    unnest(!!funcs) %>% 
    hoist(m1, pm1 = "p") %>% 
    hoist(m2, pm2 = "p") %>% 
    hoist(m3, pm3 = "p") %>% 
    hoist(m5, pm5 = "p")  %>% 
    discard(is.list)
  
  # mediation
  tab1b = 
    example0 %>% 
    hoist(out, "result") %>% 
    hoist(result, !!!funcs)  %>% 
    unnest(m99) %>%
    unnest_wider(m99) %>% 
    hoist(w5bmi, w5bmi_p = c("result", "p"))  %>% 
    hoist(bingedrink, bingedrink_p = c("result", "p"))  %>% 
    hoist(currentsmoke, currentsmoke_p = c("result", "p"))  %>% 
    hoist(phys_activ_ff5, phys_activ_ff5_p = c("result", "p"))  %>%   
    discard(is.list)
  
  (
    tab1a %>% left_join(tab1b)
  )
  
}

if(example1){
  
  ############################################################
  # whole genome
  ############################################################
  
  example1 = 
    args %>% 
    filter(gene_set_name %>% str_detect("whole")) %>%
    slice(1:3) %>%   
    mutate(out = pmap(., safely(model_fit), funcs), 
           controls = names(controls))
  
  # no errors
  example1 %>%
    hoist(out, "error") %>% 
    mutate(error = map(error, as.character)) %>%
    unnest(error) 
  
  # results
  example1 %>% hoist(out, c("result")) %>% pluck("result")
  
}
