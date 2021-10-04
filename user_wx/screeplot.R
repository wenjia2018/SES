threshold = 0.05
# rle
# example0_with1k <- readRDS("~/ses-1/user_wx/example_RLE_pca_nomed_withinflame.rds")
# example0_without1k<- readRDS("~/ses-1/user_wx/example_RLE_pca_nomed_noinflame.rds")
# tmm
example0_with1k <- readRDS("~/ses-1/user_wx/example_tmm_m7_withinflame_nonerotationonlypca.rds")

example0_with1k %>%
  hoist(out, varexplained = list("result", "m7_nn", 1, "other", "varexplained")) %>% 
  select(gene_set_name, varexplained) %>% 
  unique() %>% 
  unnest_longer(varexplained) %>% 
  mutate(varexplained_id = rep(str_c("d",1:10), 11)) %>% 
  ggplot(aes(x = factor(varexplained_id, level = str_c("d",1:10)), y = varexplained)) +
  geom_point() +
  facet_wrap(facets = "gene_set_name") +
  scale_y_continuous(labels = scales::percent)

example0_with1k <-readRDS("~/ses-1/user_wx/example_tmm_m7_withinflame.rds")
example0_without1k <- readRDS("~/ses-1/user_wx/example_tmm_m7_noinflame0209.rds")
example0_with1k %>%
  hoist(out, varexplained = list("result", "m7_ob", 1, "other", "varexplained")) %>% 
  select(gene_set_name, varexplained) %>% 
  unique() %>% 
  unnest_longer(varexplained) %>% 
  mutate(varexplained_id = rep(str_c("d",1:10), 11)) %>% 
  ggplot(aes(x = factor(varexplained_id, level = str_c("d",1:10)), y = varexplained)) +
  geom_point() +
  facet_wrap(facets = "gene_set_name")+
  scale_y_continuous(labels = scales::percent)






example0_without1k %>%
  hoist(out, varexplained = list("result", "m7_ob", 1, "other", "varexplained")) %>% 
  select(gene_set_name, varexplained) %>% 
  unique() %>% 
  filter(gene_set_name!="inflam1k_mRNA") %>% 
  unnest_longer(varexplained) %>% 
  mutate(varexplained_id = rep(str_c("d",1:10), 10)) %>% 
  ggplot(aes(x = factor(varexplained_id, level = str_c("d",1:10)), y = varexplained)) +
  geom_point() +
  facet_wrap(facets = "gene_set_name")+
  scale_y_continuous(labels = scales::percent)
