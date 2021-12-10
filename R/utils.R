mediate_spca = function(mediator, gene_set,pca_out){
  datt_m = 
    prepro(gene_set, treatment, c(controls, mediator)) %>%
    rename(treatment = treatment
           # mediator = mediator
    ) 
  if(remove_diseased_subjects) {
    datt_m =
      datt_m %>% 
      remove_diseased_subjects_from_datt(gene_set_name, controls)
  }

  pca_sparse <- sparsepca::spca(datt_m[gene_set], k = ncomp, scale = TRUE)
  pca_sparse$scores = pca_sparse$scores %>%  as_tibble() %>%  set_names(str_c("d", 1:ncomp))
  
  datt_m = dplyr::select(datt_m, -gene_set) 
  outcome = colnames(pca_sparse$scores) %>% set_names()
  
  convert_outcome = function(outcome, datt_m){
    datt_pca = bind_cols(datt_m, !!outcome := pca_sparse$scores[, outcome])
  }
  
  datt_pca = outcome %>% map(convert_outcome, datt = datt_m)
  # only do mediation for significant PCs, p values are corrected by bonferonni
  # actually should divided by 4 the treatment numbers, but Mike sometimes keep some
  # signatures when their p value is close and a little bigger to the cut off, so here
  # choose correcting pca component numbers only so mediational analysis will be done to less cases, save time.
  threshold = 0.05 / length(pca_out$p)
  sig = pca_out$p < 0.01
  print("please be aware: p value correction is bonferonni for now")
  # map2(datt_pca[sig], outcome[sig], fit_m99) %>% map(extract_m99)
  pmap(list(datt_pca[sig], outcome[sig], mediator), fit_m99) %>% map(extract_m99)
}

mediate_multiple_mediators = function(gene_set, rotate, pca_out){
  datt_m = 
    prepro(gene_set, treatment, c(controls, all_of(mediators))) %>%
    rename(treatment = treatment
           # mediator = mediator
    ) %>% 
    select(all_of(mediators), everything())
  
  if(remove_diseased_subjects) {
    datt_m =
      datt_m %>% 
      remove_diseased_subjects_from_datt(gene_set_name, controls)
  }
  pca_rotated = fit_pca_util(datt_m, gene_set, rotate) 
  
  datt_m = dplyr::select(datt_m, -gene_set) 
  outcome = colnames(pca_rotated$scores) %>% set_names()
  
  convert_outcome = function(outcome, datt_m){
    datt_pca = bind_cols(datt_m, !!outcome := pca_rotated$scores[, outcome])
  }
  
  datt_pca = outcome %>% map(convert_outcome, datt = datt_m)
  # only do mediation for significant PCs, p values are corrected by bonferonni
  # actually should divided by 4 the treatment numbers, but Mike sometimes keep some
  # signatures when their p value is close and a little bigger to the cut off, so here
  # choose correcting pca component numbers only so mediational analysis will be done to less cases, save time.
  threshold = 0.05 / length(pca_out$p)
  sig = pca_out$p < 0.01
  print("please be aware: p value correction is bonferonni for now")
  # map2(datt_pca[sig], outcome[sig], fit_m99) %>% map(extract_m99)
  pmap(list(datt_pca[sig], outcome[sig]), mediate_multimediate) %>% map(~.x %>% summary(opt = "avg"))
}


# mediate_multimediate for multiple mediators, need to change the mediators and functions accordingly
mediate_multimediate = function(datt, gene_set) {
  datt = datt %>% na.omit()
  M1reg = lm(as.formula(str_c("w5bmi_lm ~ treatment +", controls %>% str_c(collapse = " + "))), datt)
  M2reg = lm(as.formula(str_c("stress_perceived_lm ~ treatment +", controls %>% str_c(collapse = " + "))), datt)
  if(COR){
    datt = datt %>% mutate_at(.vars = vars(c("bills_binary", "currentsmoke_binary","insurance_lack_binary")),
                              .funs = list(~ .x %>% as.numeric()))
    M3reg = lm(as.formula(str_c("bills_binary ~ treatment +", controls %>% str_c(collapse = " + "))), datt)
    M4reg = lm(as.formula(str_c("currentsmoke_binary ~ treatment +", controls %>% str_c(collapse = " + "))), datt)
    M5reg = lm(as.formula(str_c("insurance_lack_binary ~ treatment +", controls %>% str_c(collapse = " + "))), datt)
    
  } else{
    M3reg = glm(as.formula(str_c("bills_binary ~ treatment +", controls %>% str_c(collapse = " + "))), datt, family = binomial("logit"))
    M4reg = glm(as.formula(str_c("currentsmoke_binary ~ treatment +", controls %>% str_c(collapse = " + "))), datt, family = binomial("logit"))
    M5reg = glm(as.formula(str_c("insurance_lack_binary ~ treatment +", controls %>% str_c(collapse = " + "))), datt, family = binomial("logit"))
    
  }

  Yreg = lm(as.formula(str_c(gene_set, "~ .")), datt)
  med_multimediate = multimediate::multimediate(lmodel.m = list(M1reg, M2reg, M3reg, M4reg, M5reg),
                                            correlated = COR,
                                            model.y = Yreg,
                                            treat = "treatment",
                                            treat.value = 1,
                                            control.value = 0,
                                            J = 1000,
                                            conf.level = 0.95)
}


mediate_mma = function(datt, gene_set) {
  x = datt %>% select(all_of(mediators), all_of(controls))
  contmed =  which(colnames(x)%like%"_lm")
  binmed = which(colnames(x)%like%"_binary")
  med_index = 1:length(mediators)
  y = datt %>% select(!!gene_set)
  m1 <- mma::mma(x = x, y = y, pred = treatment, contmed = contmed, binmed = binmed, jointm = list(n = 1, j1 = med_index),
            n2 = 1000, alpha = 0.05 )
  
}

cal_evalue = function(x) {
  ols_result = x %>% summary %>% broom::tidy() %>% filter(str_detect(term, "treatment"))
  
  evalues = EValue::evalues.OLS(est = ols_result %>% pull("estimate"),
                                se = ols_result %>% pull("std.error"),
                                sd = summary(x)$sigma)
  
}


DDA = function(mediator, gene_set, rotate, pca_out) {
  datt_m = 
    prepro(gene_set, treatment, c(controls, mediator)) %>%
    rename(treatment = treatment
           # mediator = mediator
    ) 
  
  if(remove_diseased_subjects) {
    datt_m =
      datt_m %>% 
      remove_diseased_subjects_from_datt(gene_set_name, controls)
  }
  pca_rotated = fit_pca_util(datt_m, gene_set, rotate) 
  
  datt_m = dplyr::select(datt_m, -gene_set) 
  outcome = colnames(pca_rotated$scores) %>% set_names()
  
  convert_outcome = function(outcome, datt_m){
    datt_pca = bind_cols(datt_m, !!outcome := pca_rotated$scores[, outcome])
  }
  
  datt_pca = outcome %>% map(convert_outcome, datt = datt_m)
  # only do mediation for significant PCs, p values are corrected by bonferonni
  # actually should divided by 4 the treatment numbers, but Mike sometimes keep some
  # signatures when their p value is close and a little bigger to the cut off, so here
  # choose correcting pca component numbers only so mediational analysis will be done to less cases, save time.
  threshold = 0.05 / length(pca_out$p)
  sig = pca_out$p < threshold
  print("please be aware: p value correction is bonferonni for now")
  
  covariates = colnames(datt_m) %>% str_c(collapse= " + ")
  formula_lm = names(datt_pca) %>% map(~ str_c(.x,  "~", covariates)) %>% map(as.formula)

  aa = map2(datt_pca[sig], formula_lm[sig], ~ DDA_component(model_formula = .y, predictor = mediator, data = .x, n_boot = n_boot)) 
  
  
}


DDA_component = function(model_formula, predictor, data, n_boot = n_boot) {
  
  vardist = dda.vardist(model_formula, pred = predictor, data = data, B = n_boot)
  
  # resdist = dda.resdist(model_formula, pred = predictor, data = data, B = n_boot)
  
  indep1 = dda.indep(model_formula, pred = predictor, data = data)
  # 
  indep2 = dda.indep(model_formula, pred = predictor, nlfun = tanh, data = data)
  
  # indep3 = dda.indep(model_formula, pred = predictor, data = data, hsic.method = "gamma", nlfun = 2) 
  
  indep4 = dda.indep(model_formula, pred = predictor, data = data, hsic.method = "boot", B = n_boot, nlfun = 2, parallelize = TRUE, cores = 4)
  
  return(list(
    vardist = vardist,
    # resdist = resdist, 
    indep1 = indep1,
    indep2 = indep2,
    # indep3 = indep3
    indep4 = indep4
    ))
  
}



mediate = function(mediator, gene_set, controls, treatment){
  datt_m = 
    prepro(gene_set, treatment, c(controls, mediator)) %>%
    rename(treatment = treatment,
           mediator = mediator) 
  
  if(remove_diseased_subjects) {
    datt_m =
      datt_m %>% 
      remove_diseased_subjects_from_datt(gene_set_name, controls)
  }
  ok = complete.cases(datt_m$mediator)
  datt_m = datt_m[ok, ]
  fit_m99(datt_m, gene_set, mediator) %>% extract_m99()  
}

# mediation for each gene

mediate_multiple = function(controls, treatment, gene_set){ 
  out = NULL
  args_m97 = 
    crossing(mediators, gene_set) %>% 
    rename(mediator = mediators) %>% 
    mutate(treatment = treatment,
           controls = list(controls)) %>% 
    mutate(names = str_c(mediator, gene_set, sep = "_"))
  out = pmap(list(args_m97$mediator,
       args_m97$gene_set,
       args_m97$controls,
      args_m97$treatment),
       safely(mediate)) %>% 
    set_names(args_m97$names)
  return(out)
  
}

mediate_pca = function(mediator, gene_set, rotate, pca_out){
  datt_m = 
    prepro(gene_set, treatment, c(controls, mediator)) %>%
    rename(treatment = treatment
           # mediator = mediator
           ) 
  if(remove_diseased_subjects) {
    datt_m =
      datt_m %>% 
      remove_diseased_subjects_from_datt(gene_set_name, controls)
  }
  pca_rotated = fit_pca_util(datt_m, gene_set, rotate) 
  
  datt_m = dplyr::select(datt_m, -gene_set) 
  outcome = colnames(pca_rotated$scores) %>% set_names()
  
  convert_outcome = function(outcome, datt_m){
    datt_pca = bind_cols(datt_m, !!outcome := pca_rotated$scores[, outcome])
  }
  
  datt_pca = outcome %>% map(convert_outcome, datt = datt_m)
  # only do mediation for significant PCs, p values are corrected by bonferonni
  # actually should divided by 4 the treatment numbers, but Mike sometimes keep some
  # signatures when their p value is close and a little bigger to the cut off, so here
  # choose correcting pca component numbers only so mediational analysis will be done to less cases, save time.
  threshold = 0.05 / length(pca_out$p)
  sig = pca_out$p < 0.01
  print("please be aware: p value correction is bonferonni for now")
 # map2(datt_pca[sig], outcome[sig], fit_m99) %>% map(extract_m99)
 pmap(list(datt_pca[sig], outcome[sig], mediator), fit_m99) %>% map(extract_m99)
}

get_table1 = function(example){ 
  
  tab1a = 
    example %>% 
    hoist(result, !!!funcs)  %>% 
    unnest(!!funcs) %>% 
    # hoist(m1, pm1 = "p") %>% 
    # hoist(m2, pm2 = "p") %>% 
    # hoist(m3, pm3 = "p") %>% 
    # hoist(m5, pm5 = "p") %>%
    # hoist(m5b, pm5b = "p") %>%
    hoist(m8_fwer, m8_fwer_p = "p") %>% 
    hoist(m8_fdr, m8_fdr_p = "p") %>% 
    discard(is.list)
  
  if(0){ 
    
    # mediation
    tab1b = 
      example %>% 
      hoist(result, !!!funcs)  %>% 
      unnest(m99) %>%
      unnest_wider(m99) %>% 
      hoist(w5bmi, w5bmi_p = c("result", "p"))  %>% 
      hoist(bingedrink, bingedrink_p = c("result", "p"))  %>% 
      hoist(currentsmoke, currentsmoke_p = c("result", "p"))  %>% 
      hoist(phys_activ_ff5, phys_activ_ff5_p = c("result", "p"))  %>%   
      discard(is.list)
    
    tab1a %>% left_join(tab1b)
    
  } else {
    
    tab1a
    
  }
}

detach_package <- function(abcpkg, character.only = FALSE)
  # https://stackoverflow.com/questions/6979917/how-to-unload-a-package-without-restarting-r
{
  if(!character.only)
  {
    abcpkg <- deparse(substitute(abcpkg))
  }
  search_item <- paste("package", abcpkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}


remove_errors = function(example){
  
  # ERRORS:
  example %>%
    hoist(out, "error") %>%
    mutate(error = map(error, as.character)) %>%
    unnest(error) %>%
    group_by(error) %>%
    slice(1) %>% 
    print()
  
  # WHAT CAUSES ERROR? RELATE NA TO ARGS OF model_fit()
  example %>%
    hoist(out, p = list("result", "m5", 1, "p")) %>%
    with(table(gene_set_name, is.na(p))) 
  
  # REMOVE MODELS THAT ERR
  example = example %>% hoist(out, "result") %>% drop_na()
}

# https://stats.stackexchange.com/questions/108007/correlations-with-unordered-categorical-variables

catcorrm <- function(vars, dat) sapply(vars, function(y) sapply(vars, function(x) vcd::assocstats(table(dat[,x], dat[,y]))$cramer))

# https://stackoverflow.com/questions/41582486/how-to-convert-r-mediation-summary-to-data-frame

# broom::tidy does a similar job but without total effect and prop 
mediate_summary = function(med.out) {
  summary(med.out) %>%
    capture.output() %>%
    discard(`==`, "") -> lines
  
  lines[which(grepl("^ ", lines)):(which(grepl("^Sample", lines))-3)] %>%
    sub("^       ", "med.out", .) %>%
    gsub("95% ", " 95%", .) %>%
    gsub("\\*", "",.) %>%
    gsub("\\.$","",.) %>% 
    gsub("CI ", "CI_", .) %>%
    sub(" \\(", "_(", .) %>%
    sub("p-", "p_", .) %>%
    sub(" ", "_", ., fixed=TRUE)  %>%
    textConnection() %>%
    read.table(header=TRUE) %>%
    setNames(sub("_$", "", colnames(.))) %>%
    dplyr::mutate(med.out=sub("\\.|_$", "", med.out),
                  med.out=gsub("_", " ", med.out))
  
  
}
