mediate = function(mediator, gene_set){
  
  datt_m = 
    prepro(gene_set, treatment, c(controls, mediator)) %>%
    rename(treatment = treatment,
           mediator = mediator) %>% 
    remove_diseased_subjects_from_datt(gene_set_name, controls)
  
  ok = complete.cases(datt_m$mediator)
  datt_m = datt_m[ok, ]
  fit_m99(datt_m, gene_set) %>% extract_m99()  
}

mediate_multiple = function(gene_set){ 
  
  args_m97 = 
    crossing(mediators, gene_set) %>% 
    rename(mediator = mediators) %>% 
    mutate(names = str_c(mediator, gene_set, sep = "_"))
  map2(args_m97$mediator, args_m97$gene_set, safely(mediate)) %>% 
    set_names(args_m97$names)
  
}

mediate_pca = function(mediator, gene_set, rotate){
  
  datt_m = 
    prepro(gene_set, treatment, c(controls, mediator)) %>%
    rename(treatment = treatment,
           mediator = mediator) %>% 
    remove_diseased_subjects_from_datt(gene_set_name, controls)
  
  pca_rotated = fit_pca_util(datt_m, gene_set, rotate, ncomp = 9) 
  
  datt_m = dplyr::select(datt_m, -gene_set) 
  gene_set = colnames(pca_rotated$scores)
  outcome = gene_set %>% set_names()
  
  convert_outcome = function(outcome, datt_m){
      datt_pca = bind_cols(datt_m, !!outcome := pca_rotated$scores[, outcome])
    }
  
  datt_pca = outcome %>% map(convert_outcome, datt = datt_m)

  out = map2(datt_pca[1:2], outcome[1:2], fit_m99) %>% map(extract_m99)  
  return(out)
}
# for each of the pca, do mediational analysis for each pc
out$m7_nn$mediation = mediate_pca(mediator, gene_set, "none")
out$m7_vx$mediation = mediate_pca(mediator, gene_set, "varimax")
out$m7_ob$mediation = mediate_pca(mediator, gene_set, "oblimin")


out$m7_nn$mediation = fit_m7(datt, gene_set, "none")$pca$scores %>% mediate_multiple(gene_set)
out$m7_vx$mediation = fit_m7(datt, gene_set, "varimax")$pca$scores %>% mediate_multiple(gene_set)
out$m7_ob$mediation = fit_m7(datt, gene_set, "oblimin")$pca$scores %>% mediate_multiple(gene_set)

out$varexplained = pca_rotated$Vaccounted[4,] %>% as.list()

out$pca = pca_rotated

out

############################################################
# cibersort cell type compositional analysis
############################################################

fit_m96 =function(treatment, controls){
  
  cell_types <- c(
    "B.cells.naive", "B.cells.memory", "Plasma.cells",
    "T.cells.CD8", "T.cells.CD4.naive", "T.cells.CD4.memory.resting",
    "T.cells.CD4.memory.activated",
    "T.cells.follicular.helper",
    "T.cells.regulatory..Tregs.", "T.cells.gamma.delta",
    "NK.cells.resting", "NK.cells.activated", "Monocytes", "Macrophages.M0",
    "Macrophages.M1",
    "Macrophages.M2", "Dendritic.cells.resting",
    "Dendritic.cells.activated", "Mast.cells.resting",
    # "Mast.cells.activated", # all zeros
    "Eosinophils", "Neutrophils"
  )
  
  term <- c(controls$basic, treatment)
  term <- c(controls, treatment)  # inside function
  rhs <- str_c(term, collapse = " + ")
  keep <- pData(dat) %>%
    dplyr::select(all_of(term)) %>%
    complete.cases() # keep only the complete cases
  nsubs <- sum(keep) # number of complete cases
  
  ########################################################
  # DEFINE DEPENDENT VARIABLE SETS AND LOAD INTO THE GLOBAL ENV
  ########################################################
  
  outcomes <- pData(dat) %>%
    dplyr::select(all_of(cell_types)) %>%
    compositions::acomp() %>%
    compositions::ilr()
  
  m <- lm(str_c("outcomes[keep, ] ~ ", rhs), data = pData(dat)[keep, ])
  
  # partial effect
  b <- coef(m)[treatment, ] %>% compositions::ilrInv()
  names(b) <- cell_types
  
  out <- list(b = coef(m)[treatment, ] %>% compositions::ilrInv(),
              Manova = m %>% car::Manova()
  )
  
  return(out)
}




extract_m96 = function(m, out = NULL){
  
  extract_anova = function(x) anova(x) %>% tidy %>% filter(str_detect(term, "treatment"))
  
  out$detail = extract_anova(m$Manova)
  out$p = extract_anova(m$Manova)$p.value
  out$other$b = m$b
  
  return(out = out)
}
if(is.element("m96", funcs))  out$m2  = fit_m96(treatment, controls) %>% extract_m96()