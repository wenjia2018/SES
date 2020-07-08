library(tidyverse)
library(furrr)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase) 
library(dbr) # my package
walk(dir(path = "R",full.names = TRUE), source)

############################################################
# LOAD DATA (reconciled = ideal setting)
############################################################

load_data(reconciled = FALSE) %>% 
  list2env(.GlobalEnv)

recode_variables_in_dat() 

source("user_ms/define_treatments_controls_outcomes.R")

############################################################
# EXAMPLE 1: ESTIMATE WHOLE GENOME DE AND TFBM
############################################################

example1 = 
  args %>% 
  filter(treatment == "ses_sss_composite", 
         gene_set_name == "whole_genome_and_tfbm", 
         names(controls) == "basic") %>% 
  get_results(permutation = FALSE)

# INSPECT DE 
example1 %>% pluck("out", "result", "ttT")  
example1 %>% pluck("out", "result", "ttT") %>%  filter(adj.P.Val < 0.05) %>% pluck("gene") # which genes have corrected p-value < 0.05

# INSPECT TFBM: uncorrected p values
example1 %>% pluck("out", "result", "tfbm")  # TFBM de Novo
example1 %>% pluck("out", "result", "tfbm")  %>% filter(tfbm %in% immune_tfbms) # Inflammation

############################################################
# EXAMPLE 2: ESTIMATE UNIVARIATE/MULTIVARIATE MODELS FOR TABLE 1
############################################################

example2 = 
  args %>% 
  filter(treatment == "ses_sss_composite", 
         is.element(gene_set_name, table1), 
         names(controls) == "basic") %>%
  sample_n(1) %>% 
  get_results(permutation = FALSE) 

# INSPECT
# what do the column abbreviations of the following table mean?
example2 %>% pluck("out", "result", "table")
example2 %>% 
  filter(out_id == "result") %>%
  unnest(out) %>% 
  unnest(contains("m")) %>% 
  mutate(nm = names(m1)) %>% # any m would do
  filter(nm == "p") %>% 
  unnest(contains("m")) 

############################################################
# diagnose errors in any model?
############################################################

(errors =
   example2 %>% 
   filter(out_id == "error") %>%
   filter(map_lgl(out, negate(is.null)))) 

errors %>% pluck("out")

############################################################
# EXAMPLE 3: SIMILAR TO EXAMPLE 2 BUT WITH ALL CONTROLS AND CONSTITUENT TREATMENTS
############################################################

example3 = 
  args %>% 
  filter(is.element(gene_set_name, table1)) %>%  
  sample_n(10) %>% 
  get_results(permutation = TRUE, n_perm = 1000)


# scalar summaries
example3 %>%   
  filter(out_id == "result") %>%
  unnest(out) %>% 
  unnest(contains("m")) %>% 
  mutate(nm = names(m1)) %>% # any m would do
  filter(nm == "p") %>% 
  unnest(contains("m")) 

# details for m1
example3 %>%   
  filter(out_id == "result") %>%
  unnest(out) %>% 
  unnest(contains("m")) %>% 
  mutate(nm = names(m1)) %>% # any m would do
  filter(nm == "detail") %>% 
  unnest(m1)

# Pillai from m1
example3 %>%   
  filter(out_id == "result") %>%
  unnest(out) %>% 
  unnest(contains("m")) %>% 
  mutate(nm = names(m1)) %>% 
  filter(nm == "other") %>% # hack to get summaries which have length 1
  unnest(m1) %>% 
  mutate(Pillai = map_dbl(m1, ~.x))

saveRDS(example3, "rds/example3.rds")