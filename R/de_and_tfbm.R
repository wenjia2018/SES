de_and_tfbm <- function(treatment, controls, de_only = FALSE) {
  
  # Specify whole-genome regression of rna on design
  y <- dat %>% Biobase::exprs()
  X <- dat %>% Biobase::pData() %>% dplyr::select(all_of(controls), all_of(treatment))
  
  keep = X %>% complete.cases()
  X = X[keep, ]
  y = y[, keep]
  
  X = model.matrix(~ ., data = X)
  # drop columns to make sure X is full rank
  X = ordinal::drop.coef(X)
  
  # MatrixModels::model.Matrix can drop unused levels, but has error when the variable names contains some special characters 
  # https://stackoverflow.com/questions/44114391/error-in-sparse-model-matrix
  # X = MatrixModels::model.Matrix(~., data = X, sparse = TRUE, drop.unused.levels = TRUE)
  # needed for limma::topTables() to work with factors
  treatment = X %>% colnames() %>% str_detect(treatment) %>% which()
  

  # Estimate DE using standard limmma/edger pipeline. 
  ttT <-
    lmFit(y, X) %>%
    eBayes %>%
    tidy_topTable(of_in = treatment, confint = TRUE)
  
  if(de_only) return(list(ttT = ttT))
  # gene set tests self contained
  sets = signatures$outcome_set[table1]

  gene_set_test = mroast(y, index = sets, design = X, contrast = treatment) %>% rownames_to_column("gene_set_name")
  # optional gene.weights = ttT$logFC
  # ttT = ttT %>% mutate(weight = ifelse(logFC>0,1,-1)) gene.weights = ttT$weight
  
  # fig1 = DE_enrichplot(ttT)
  
  # genes whose uncorrected p-values below 0.05 (not an inference):
  ttT_sub = filter(ttT, P.Value <= 0.05)
  
  tfbm = 
    ttT %>%
    infer_db_adapted(ttT_sub = ttT_sub) %>%
    # extract_db  %>% 
    pluck("telis", "par", "p_over") %>% # tfbms over represented in DNA of genes labeled with high mRNA relation to treatment
    enframe(name = "tfbm", value = "tellis p.value for over-represented tfbms")
  
  return(list(ttT = ttT, tfbm = tfbm, gene_set_test = gene_set_test))
}