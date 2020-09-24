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

example0 = readRDS("/home/share/scratch/fig1A_ses4binary_with1KI.rds")

with1KI = example0 %>%
  hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
  filter(controls=="all") %>% 
  dplyr::select(treatment, gene_set_name, p) %>% 
  rename(p_with1KI = p)




example0 = readRDS("/home/share/scratch/fig1A_ses4binary_without1KI.rds")

without1KI = example0 %>%
       hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
       filter(controls=="all") %>% 
       dplyr::select(treatment, gene_set_name, p) %>% 
  rename(p_without1KI = p)

a = with1KI %>% left_join(without1KI, by = c("treatment","gene_set_name"))

a %>% kableExtra::kable() %>% kableExtra::kable_styling()
