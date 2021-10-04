library(tidyverse)
example0_with1k <- readRDS("~/ses-1/user_wx/example_tmm_m12_withinflame.rds")
tp =
  example0_with1k %>%
  hoist(out, p = list("result", "m12_fdr", 1, "p")) %>% 
  hoist(out, t = list("result", "m12_fdr", 1, "detail","t")) %>% 
  filter(p<0.05) %>% 
  select(-out)
