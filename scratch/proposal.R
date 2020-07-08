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
example_wx = 
  args %>% 
  filter(gene_set_name %>% str_detect("whole")) %>%
  slice(1:3) %>%   
  mutate(out = pmap(., safely(model_fit), funcs), 
         controls = names(controls))


# no errors
example_wx %>%
  hoist(out, "error") %>% 
  mutate(error = map(error, as.character)) %>%
  unnest(error) 

# results
example_wx %>% hoist(out, c("result")) %>% pluck("result")

