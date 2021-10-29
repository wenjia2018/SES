mediate_cluster <- readRDS("~/ses-1/user_wx/mediate_cluster_with1k_drinkmed.rds")
mediate_cluster <- readRDS("~/ses-1/user_wx/mediate_cluster_without1k_drinkmed.rds")
mediate_cluster <- readRDS("~/ses-1/user_wx/mediate_cluster_Denovo_drinkmed.rds")

mediators = 
  c(
"drink_category"
  )

extract_pp = function(m, mediator){
  m %>%
    hoist(out, result = list("result", mediator, "result")) %>% 
    filter(result!="NULL") %>% 
    unnest_longer(result) %>% 
    hoist(result, p = "p")%>% 
    hoist(result, detail = "other") %>% 
    unnest_wider("detail") %>% 
    dplyr::select(-out, -result) %>% 
    mutate(p = p %>% str_remove("<") %>% as.numeric() %>% format(digits = 3, scientific =T),
           mediator = mediator)
  # %>% 
  # adjP = ((p %>% str_remove("<") %>% as.numeric())) %>% p.adjust("fdr")) %>% 
  # filter(adjP<0.05) %>% 
  # mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric
}
med_out_all = function(example_data) {
  out = mediators %>% set_names() %>%  map(~ extract_pp(example_data, .x)) %>% bind_rows() 
  # %>% 
    # mutate(adjP = ((p %>% str_remove("<") %>% as.numeric())) %>% p.adjust("fdr")) %>% 
    # filter(adjP<0.05) %>%
    # mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric
  # out %>% select(-id, -controls) %>% kableExtra::kable() %>% kableExtra::kable_styling()
}


exB_data_with = 
  mediate_cluster %>%
  med_out_all() %>%
  mutate("1KI Genes" = "With 1KI Genes" %>% as.factor(),
         proportion = NA)
