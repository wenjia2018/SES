bespoke_m7nnfun = function(example) {
  example = example%>% 
    hoist(out, out1 = list("result", "example0")) %>%
    select(out1) %>% 
    unnest(out1) %>%
    hoist(out, p = list("result", "m7_nn", 1, "p")) %>% 
    unnest_longer(p) %>% 
    group_by(treatment) %>% 
    arrange(p) %>% 
    mutate(p.adj = p.adjust(p, method = "fdr")) %>% 
    filter(p.adj <0.05)
  
}
example <- readRDS("~/ses-1/user_wx/m7nn_with1k_aging_sc5levels_allsample_bespoke.rds")
with = bespoke_m7nnfun(example)
source("/home/xu/ses-1/user_wx/extract_v2.R")
# in extract_v2 now is m7_nn could be changed to m7_ob or other rotation
result_med = outm7med(p = 0.05, mediators, control = "ancestryPC", example) 
temp = with %>% left_join(result_med, by = c("treatment", "gene_set_name")) %>% select(-out)
