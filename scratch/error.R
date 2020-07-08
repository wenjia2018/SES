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
library(dbr) # my package
walk(dir(path = "R",full.names = TRUE), source)
fit_m4 = partial(fit_m4, n_perm = 1000) # specify n_perm

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
funcs = c("m6", "m7")
funcs = c("m99")
 
example3 = 
  args %>% 
  filter(is.element(gene_set_name, table1)) %>% 
  slice(1) %>% 
  mutate(out = pmap(., safely(model_fit), funcs), 
         controls = names(controls))

############################################################
# Extract: behold the NAs!
############################################################

example3 %>% 
  hoist(out, "result") %>% 
  hoist(result, !!!funcs)  %>% 
  unnest(m99) %>%
  unnest_wider(m99) %>% 
  hoist(w5bmi, w5bmi_p = c("result", "p"))  %>% 
  hoist(bingedrink, bingedrink_p = c("result", "p"))  %>% 
  hoist(currentsmoke, currentsmoke_p = c("result", "p"))  %>% 
  hoist(phys_activ_ff5, phys_activ_ff5_p = c("result", "p"))  %>%   
  discard(is.list)
