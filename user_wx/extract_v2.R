
f1 = function(p, data){
  ex1 <- data %>%
    hoist(out, out1 = list("result","example1")) %>%
    select(1:3) %>%
    filter(p_eqtl == p) %>% 
    rename(out=out1) %>% 
    unnest(out) 
  # %>% 
  #   filter(control_set =="ancestryPC") 
  return(ex1)
}

outttT = function(p, control, data){
  out = f1(p, data) %>% 
    hoist(out, ttT = list("result", "ttT")) %>% 
    filter(map_lgl(.$ttT, ~dim(.)[1]!=0)) %>% 
    mutate(gene_sig = map(ttT, ~ filter(., adj.P.Val<0.05) %>% pull(gene))) %>% 
    filter(control_set == control) %>% 
    select(-out, -gene_set_name)
  return(out)
}

f0 = function(p, data){
  ex0 <- data %>%
    hoist(out, out1 = list("result","example0")) %>%
    dplyr::select(1:3) %>%
    filter(p_eqtl == p) %>% 
    dplyr::rename(out=out1) %>%
    unnest(out) 
  # %>% 
  # filter(control_set =="ancestryPC")
  return(ex0)
}



m8_present = function(p, control, data){
  
  out = f0(p, data) %>% 
    hoist(out, p_omnibus = list("result", "m8_fdr", 1, "p")) %>% 
    hoist(out, logFC = list("result", "m8_fdr", 1, "detail", "logFC")) %>% 
    hoist(out, gene_sig = list("result", "m8_fdr", 1, "detail", "gene")) %>% 
    dplyr::select(p_eqtl, treatment, gene_set_name, p_omnibus, logFC, gene_sig, control_set) %>% 
    # dplyr::filter(p<0.05) %>%
    filter(control_set == control)
  return(out)
}

# 
# outm8 = function(p, control, data){
#   m8_present(f0(p, data), control) %>% mutate(p_eqtl = p)
# }


outm8_allsig = function(p, control, data){
  out = f0(p, data) %>% 
    hoist(out, ttT = list("result", "m8_fdr", 1, "other", "m")) %>% 
    filter(ttT != "NULL") %>% 
    mutate(gene_sig = map(ttT, ~ dplyr::filter(., adj.P.Val<0.05) %>% pull(gene))) %>% 
    # select(-out) %>% 
    filter(control_set == control)
}

# med_extact_m8 = function (control, ex0){
#   out = ex0 %>%
#     hoist(out, med=list("result", "m8_fdr", 1, "mediation_mean")) %>% 
#     # mutate(mediator = focal) %>% 
#     filter(control_set== control) %>% 
#     dplyr::select(-out)
# }

outm8med = function(p, control, data){
  out = f0(p, data) %>%
    hoist(out, med=list("result", "m8_fdr", 1, "mediation_mean")) %>% 
    filter(control_set== control) %>% 
    dplyr::select(-out) %>% 
    filter(med!="NULL") %>%
    unnest_longer(med) %>%
    hoist(med, p = list("result","p")) %>%
    filter(!is.na(p)) %>% 
    # filter(p < threshold_med) %>% 
    dplyr::rename(p_med= p)
}

m7_present = function(ex0, control){
  
  var = ex0 %>%
    hoist(out, var_explained = list("result", "m7_ob", 1, "other", "varexplained")) %>%
    dplyr::select(treatment, gene_set_name, var_explained, control_set) %>% 
    filter(var_explained!="NULL")
  # filter out the non estimable results
  var$var_explained = var$var_explained %>% map(~ set_names(.x, str_c("d", 1:length(.x))))
  
  var = var %>% unnest_longer(var_explained)
  
  
  gene_list = ex0 %>%
    hoist(out, well_loaded = list("result", "m7_ob", 1, "other", "well_loaded")) %>%
    dplyr::select(treatment, gene_set_name, well_loaded, control_set) %>% 
    filter(well_loaded!="NULL")
  gene_list$well_loaded = gene_list$well_loaded %>% map(~ set_names(.x, str_c("d", 1:length(.x))))
  
  gene_list = gene_list %>% unnest_longer(well_loaded)
  # 
  out = ex0 %>%
    hoist(out, estimate = list("result", "m7_ob", 1, "detail", "t")) %>% 
    unnest_longer(estimate) %>% 
    hoist(estimate, p = "p.value") %>% 
    hoist(estimate, coef = "estimate") %>% 
    rename(p_id = estimate_id) %>% 
    dplyr::select(treatment, gene_set_name, p, coef, p_id, control_set) %>% 
    # dplyr::filter(p < threshold) %>%
    left_join(var, by = c("treatment", "gene_set_name", "p_id"= "var_explained_id", "control_set")) %>%
    left_join(gene_list, by = c("treatment", "gene_set_name", "p_id"= "well_loaded_id", "control_set")) %>%
    filter(control_set == control) %>%
    dplyr::select(1:7) %>%
    rename(p_pca = p) 
  return(out)
  
}

outm7pca = function(p, control, data){
  m7_present(f0(p, data), control) %>% mutate(p_eqtl = p)
}

med_extact_m7 = function (focal, control, ex0){
  out = ex0 %>%
    hoist(out, med=list("result", "m7_ob", 1, "mediation", focal, "result")) %>% 
    mutate(mediator = focal) %>% 
    filter(control_set== control) %>% 
    select(-out)
}

outm7med = function(p, mediators, control, data){
  out = mediators %>% 
    map_df(med_extact_m7, control, f0(p, data)) %>% 
    filter(med!="NULL") %>%
    unnest_longer(med) %>% 
    hoist(med, p = "p") %>%
    filter(p < threshold_med) %>% 
    rename(p_med= p)
}
if(0){
  outomni = p_eqtl %>% map(m8_present, control, race_bespoke_12.02.2021)
  outomni = p_eqtl %>% map(outm7pca, control, race_bespoke_12.02.2021)
  outomni = p_eqtl %>% map(outttT, control, race_bespoke_12.02.2021)
  outomni = p_eqtl %>% map(outm8_allsig, control, race_bespoke_12.02.2021)

  temp = outomni %>%
    bind_rows() %>%
    filter(gene_sig %>% map_dfc(~ length(.))>0) %>%
    unnest_longer(gene_sig)

  temp %>%
    hoist(out, p_omnibus = list("result", "m8_fdr", 1, "p")) %>%
    hoist(out, logFC = list("result", "m8_fdr", 1, "detail", "logFC")) %>%
    hoist(out, med_single=list("result", "m8_fdr", 1, "mediation_single")) %>%
    unnest_longer(med_single) %>%
    hoist(med_single, med_single_p = list("result", "p")) %>%
    filter(!is.na(med_single_p)) %>%
    hoist(med_single, med_single_prop = list("result", "other", "med_prop")) %>%
    hoist(med_single, med_single_beta = list("result", "other", "med_beta")) %>%
    dplyr::select(-controls, -ttT, -out, -med_single, -table1) %>%
    mutate(p_eqtl = p_eqtl %>% format(scientific =T)) %>%
    filter(med_single_p<0.05) %>%
    kableExtra::kable() %>%
    kableExtra::kable_styling()

  temp %>%
    hoist(out, p_omnibus = list("result", "m8_fdr", 1, "p")) %>%
    hoist(out, logFC = list("result", "m8_fdr", 1, "detail", "logFC")) %>%
    hoist(out, med_single=list("result", "m8_fdr", 1, "mediation_signle")) %>%
    unnest_longer(med_single) %>%
    hoist(med_single, med_single_p = list("result", "p")) %>%
    filter(!is.na(med_single_p)) %>%
    hoist(med_single, med_single_prop = list("result", "other", "med_prop")) %>%
    hoist(med_single, med_single_beta = list("result", "other", "med_beta")) %>%
    dplyr::select(-controls, -ttT, -out, -med_single, -table1) %>%
    mutate(p_eqtl = p_eqtl %>% format(scientific =T)) %>%
    filter(med_single_p<0.05) %>%
    kableExtra::kable() %>%
    kableExtra::kable_styling()

  temp %>%
    hoist(out, med_mean=list("result", "m8_fdr", 1, "mediation_mean")) %>%
    unnest_longer(med_mean) %>%
    hoist(med_mean, med_mean_p = list("result", "p")) %>%
    filter(!is.na(med_mean_p)) %>%
    hoist(med_mean, med_mean_prop = list("result", "other", "med_prop")) %>%
    hoist(med_mean, med_mean_beta = list("result", "other", "med_beta")) %>%
    dplyr::select(-controls, -ttT, -out, -med_mean, -table1) %>%
    mutate(p_eqtl = p_eqtl %>% format(scientific =T)) %>%
    filter(med_mean_p<0.05) %>%
    kableExtra::kable() %>%
    kableExtra::kable_styling()
  
}
