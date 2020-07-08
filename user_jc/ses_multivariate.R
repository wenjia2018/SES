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

# fit_pca_util %>% debugonce()
# model_fit %>% debugonce()
example3 = 
  args %>% 
  filter(is.element(gene_set_name, table1)) %>% 
  mutate(out = pmap(., safely(model_fit), funcs), 
         controls = names(controls))

saveRDS(example3, "rds/example3.rds")
# example3 = readRDS("rds/example3.rds")

############################################################
# Peek
############################################################

if(0) data.tree::FromListSimple(example3)
example3 %>% hoist(out, "error") %>% pluck("error") %>% unique
# three groups of errors
example3 %>%
  hoist(out, "error") %>% 
  mutate(error = map(error, as.character)) %>%
  unnest(error) %>%
  group_by(error) %>% 
  slice(1)

example3 %>% 
  hoist(out, "error") %>% 
  mutate(error = map(error, as.character)) %>%
  unnest(error) %>%
  group_by(error) 

# where is the problem? relate NA to each of the arguments to model_fit(),
# controls, treatment, gene_set_name:
example3 %>% 
  hoist(out, p = list("result", "m1", 1, "p")) %>% 
  with(table(gene_set_name, is.na(p)))

example3 %>% hoist(out, p = list("result", "m1", 1, "p"))  

example3 %>% hoist(out, p = list("result", "m6", 1, "p")) 


############################################################
# toy
############################################################
( 
  tabPCA =
    example3 %>% 
    hoist(out, p = list("result" )) %>%
    unnest(p) %>% 
    # select(-c("m1","m2","m3","m5", "m99")) %>% 
    unnest(contains("m")) %>%
    filter(names(m6_vx) == "p") %>%
    unnest(contains("m6")) %>%
    unnest_wider(m7_nn, names_sep = "_") %>% 
    unnest_wider(m7_vx, names_sep = "_") %>% 
    unnest_wider(m7_ob, names_sep = "_")
)

tabPCA_binary = 
  tabPCA %>%
  rowwise() %>% 
  mutate(across(matches("m6|m7"), ~ as.numeric(.x < 0.05))) 

list(tabPCA = tabPCA, tabPCA_binary = tabPCA_binary) %>% 
  set_names() %>% 
  saveRDS("rds/2008.7.8_pca.rds")

example3 %>% hoist(out, p = list("result" )) %>% unnest(p) %>% unnest(contains("m")) %>% filter(names(m6_vx) == "other") %>% unnest(contains("m")) %>% unnest(contains("m"))

############################################################
# Unpack completely
############################################################

# total effect
tab1a = 
  example3 %>% 
  hoist(out, "result") %>% 
  hoist(result, !!!funcs)  %>% 
  unnest(!!funcs) %>% 
  hoist(m1, pm1 = "p") %>% 
  hoist(m2, pm2 = "p") %>% 
  hoist(m3, pm3 = "p") %>% 
  hoist(m5, pm5 = "p")  %>% 
  hoist(m6, pm6 = "p")  %>% 
  discard(is.list)

# mediation
tab1b = 
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

(
  tab1a %>% left_join(tab1b)
)

# Top 10 rotated PCs 
(
  tabPCA = 
    example3 %>% 
    hoist(out, "result") %>% 
    hoist(result, !!!funcs)  %>% 
    unnest(!!funcs) %>% 
    hoist(m7, "p") %>% 
    unnest_wider(p) %>%  
    select(treatment, gene_set_name, controls, starts_with("R")) %>% 
    discard(is.list)
)

# Top 10 unrotated PCs
(
  tabPCA = 
    example3 %>% 
    hoist(out, "result") %>% 
    hoist(result, !!!funcs)  %>% 
    unnest(!!funcs) %>% 
    hoist(m98, "p") %>% 
    unnest_wider(p) %>%  
    select(treatment, gene_set_name, controls, num_range("Dim.", 1:10)) %>% 
    discard(is.list)
)



