example0 = readRDS("/home/share/scratch/example0.rds")

example0 %>%
  filter(treatment %in% c("ses_composite_pp1","ses_sss_composite","income_hh_ff5","edu_max","SEI_ff5","sss_5")) %>% 
  filter(controls=="all") %>%
  hoist(out, "result") %>%
  hoist(result, m1 = "m1", m5 = "m5", m5b = "m5b") %>% 
  hoist("m1", pval = list(1,"p"), Pillai = list(1,"other","Pillai")) %>% 
  hoist(m5, FDR = list(1, "p")) %>% 
  hoist(m5b, FWER = list(1, "p")) %>% 
  select(1:5,7,9) %>%
  saveRDS("/home/xu/ses-1/user_wx/table1.rds")
    

example0 %>%
  filter(treatment %in% c("ses_composite_pp1","ses_sss_composite","income_hh_ff5","edu_max","SEI_ff5","sss_5")) %>% 
  filter(controls=="all") %>%
  hoist(out, "result") %>%
  hoist(result, m6_nn = "m6_nn", m7_nn = "m7_nn") %>% 
  hoist("m6_nn", pval_mancova = list(1,"p"), Pillai = list(1,"other","Pillai")) %>% 
  hoist("m7_nn",pval_uni = list(1,"p")) %>% 
  unnest_wider(pval_uni) %>% 
  select(1:5,7:15)%>%
  saveRDS("/home/xu/ses-1/user_wx/pca_nn.rds")

example0 %>%
  filter(treatment %in% c("ses_composite_pp1","ses_sss_composite","income_hh_ff5","edu_max","SEI_ff5","sss_5")) %>% 
  filter(controls=="all") %>%
  hoist(out, "result") %>%
  hoist(result, m6_vx = "m6_vx", m7_vx = "m7_vx") %>% 
  hoist("m6_vx", pval_mancova = list(1,"p"), Pillai = list(1,"other","Pillai")) %>% 
  hoist("m7_vx",pval_uni = list(1,"p")) %>% 
  unnest_wider(pval_uni) %>% 
  select(1:5,7:15)%>%
  saveRDS("/home/xu/ses-1/user_wx/pca_vx.rds")

example0 %>%
  filter(treatment %in% c("ses_composite_pp1","ses_sss_composite","income_hh_ff5","edu_max","SEI_ff5","sss_5")) %>% 
  filter(controls=="all") %>%
  hoist(out, "result") %>%
  hoist(result, m6_ob = "m6_ob", m7_ob = "m7_ob") %>% 
  hoist("m6_ob", pval_mancova = list(1,"p"), Pillai = list(1,"other","Pillai")) %>% 
  hoist("m7_ob",pval_uni = list(1,"p")) %>% 
  unnest_wider(pval_uni) %>% 
  select(1:5,7:15) %>%
  saveRDS("/home/xu/ses-1/user_wx/pca_ob.rds")


example0 %>% unnest(out) %>% `[[`("out") %>% `[[`("result")
