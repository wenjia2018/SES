library(tidyverse)
example <- readRDS("~/ses-1/user_wx/m12_with1k_aging_sccont_NonBstrata_bespoke.rds")
example <- readRDS("~/ses-1/user_wx/m12_with1k_aging_scbinary_NonBstrata_bespoke.rds")
example <- readRDS("~/ses-1/user_wx/m12_with1k_aging_sc4levels_NonBstrata_bespoke.rds")

example <- readRDS("~/ses-1/user_wx/m12_with1k_disease_sc4levels_NonBstrata_bespoke.rds")
example <- readRDS("~/ses-1/user_wx/m12_with1k_disease_scbinary_NonBstrata_bespoke.rds")
example <- readRDS("~/ses-1/user_wx/m12_with1k_disease_sccont_NonBstrata_bespoke.rds")

example <- readRDS("~/ses-1/user_wx/m13_aging_with1k_sc5levels_allsample_bespoke_boot.rds")

example = example %>% 
  hoist(out, out1 = list("result", "example0")) %>%
  select(out1) %>% 
  unnest(out1) %>%
  hoist(out, m = list("result", "m12_fdr", 1, "other", "m")) %>% 
  filter(!map_lgl(m, ~is.null(.x))) %>% 
  mutate(p = m %>% map(~ dplyr::slice(.x, which.min(adj.p.withinunion))) %>% map_dbl(~ .x %>% pluck("adj.p.withinunion"))) %>% 
  arrange(p) %>% 
  dplyr::select(treatment, gene_set_name, controls, p) %>% 
  filter(p<0.05)


example <- readRDS("~/ses-1/user_wx/m7nn_with1k_aging_scbinary_NonBstrata_bespoke.rds")
example<- readRDS("~/ses-1/user_wx/m7nn_with1k_aging_sccont_NonBstrata_bespoke.rds")
example<- readRDS("~/ses-1/user_wx/m7nn_with1k_aging_sc4levels_NonBstrata_bespoke.rds")

example <- readRDS("~/ses-1/user_wx/m7nn_with1k_disease_sccont_NonBstrata_bespoke.rds")
example <- readRDS("~/ses-1/user_wx/m7nn_with1k_disease_scbinary_NonBstrata_bespoke.rds")
example <- readRDS("~/ses-1/user_wx/m7nn_with1k_disease_sc4levels_NonBstrata_bespoke.rds")
temp = 
  example %>% 
  hoist(out, out1 = list("result", "example0")) %>%
  select(out1) %>% 
  unnest(out1) %>%
  hoist(out, p = list("result", "m7_nn", 1, "p")) %>% 
  unnest_longer(p) %>% 
  group_by(treatment) %>% 
  arrange(p) %>% 
  mutate(p.adj = p.adjust(p, method = "fdr")) %>% 
  select(-out) %>% 
  filter(p.adj <0.05)
