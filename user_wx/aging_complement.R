load_data(reconciled = FALSE, remove_inflam = FALSE)
cluster_union = reduce(names(signatures$outcome_set) %>% 
                         str_detect('cl') %>%
                         keep(signatures$outcome_set, .), union)
aging_cluster_complement = setdiff(signatures$outcome_set$aging_mRNA, cluster_union)
signatures$outcome_set$aging_cluster_complement =  aging_cluster_complement
#  check sig in aging cluster complement

temp1 = outm8_allsig (p = p_eqtl, control = "ancestryPC_ses", example_race) 
temp2 = outm8_allsig (p = p_eqtl, control = "ancestryPC_ses", example_skincolor3)
aging_all = temp1 %>% 
  rbind(temp2) %>%
  filter(gene_set_name == "aging_mRNA") %>% 
  mutate(no_sig = gene_sig %>% map_int(~ length(.x))) %>% 
  filter(no_sig>0)
cluster_complement_sig =
  aging_all %>% 
  pull(gene_sig) %>% 
  `[[`(1)
cluster_complement_sig = cluster_complement_sig[cluster_complement_sig %in% aging_cluster_complement]
aging_complement_min = aging_all$ttT[[1]] %>%
  filter(gene %in% cluster_complement_sig) %>%
  slice(which.min(adj.P.Val)) %>% 
  mutate(p_omnibus = adj.P.Val,
         gene_sig = gene) %>% 
  select(p_omnibus, logFC, gene_sig)
aging_c = aging_all %>% 
  select(p_eqtl, treatment, gene_set_name, control_set) %>%
  mutate(gene_set_name = "aging_cluster_complement_mRNA") %>% 
  cbind(aging_complement_min) %>% 
  relocate(control_set, .after = last_col()) 
