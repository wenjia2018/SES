
with1KI <- example0_with1k %>%
  hoist(out, t = list("result", "m7_ob", 1, "detail","t")) %>% 
  unnest_longer("t") %>% 
  unnest("t") %>% 
  dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  filter(p.value < 0.001) %>% 
  select(1, 2, 7, 8, 9) %>% 
  filter(gene_set_name %in% table1[1:11])
  

without1KI <- example0_without1k %>%
  hoist(out, t = list("result", "m7_ob", 1, "detail","t")) %>% 
  unnest_longer("t") %>% 
  unnest("t") %>% 
  dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  filter(p.value < 0.001) %>% 
  select(1, 2, 7, 8, 9) %>% 
  filter(gene_set_name %in% table1[1:11])



