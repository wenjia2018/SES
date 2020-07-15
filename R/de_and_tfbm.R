de_and_tfbm <- function(treatment, controls, de_only = FALSE) {
  
  # Specify whole-genome regression of rna on design
  y <- dat %>% Biobase::exprs()
  X <- dat %>% Biobase::pData() %>% select(all_of(controls), all_of(treatment))
  
  keep = X %>% complete.cases()
  X = X[keep, ]
  y = y[, keep]
  
  X = model.matrix(~ ., data = X)
  
  # needed for limma::topTables() to work with factors
  treatment = X %>% colnames() %>% str_detect(treatment) %>% which()
  
  # Estimate DE using standard limmma/edger pipeline. 
  ttT <-
    lmFit(y, X) %>%
    eBayes %>%
    tidy_topTable(of_in = treatment)
  
  if(de_only) return(ttT)
  
  fig1 = DE_enrichplot(ttT)
  
  # genes whose uncorrected p-values below 0.05 (not an inference):
  ttT_sub = filter(ttT, P.Value <= 0.05)
  
  tfbm = 
    ttT %>%
    infer_db_adapted(ttT_sub = ttT_sub) %>%
    # extract_db  %>% 
    pluck("telis", "par", "p_over") %>% # tfbms over represented in DNA of genes labeled with high mRNA relation to treatment
    enframe(name = "tfbm", value = "tellis p.value for over-represented tfbms")
  
  return(list(ttT = ttT, tfbm = tfbm, fig1 = fig1))
}