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
  # for whole genome regression without predefine any signature sets, sets will be empty,
  # add the following if else to accomendate this situation
  # now move gene set function to m10
  sets = signatures$outcome_set[table1 %>% str_subset("whole_genome",negate = T)]
  # check if any signature is not available, as mroast will not work if any of the elements in index are NULL
  ind = map_chr(sets, ~ is.null(.))
  if(TRUE %in% ind ){
    gene_set_test = NULL
  }else{
    gene_set_test = mroast(y, index = sets, design = X, contrast = treatment) %>% rownames_to_column("gene_set_name")
  }
 
  # optional gene.weights = ttT$logFC
  # ttT = ttT %>% mutate(weight = ifelse(logFC>0,1,-1)) gene.weights = ttT$weight
  
  # fig1 = DE_enrichplot(ttT)
  
  # genes whose uncorrected p-values below 0.05 (not an inference):
  ttT_sub = filter(ttT, P.Value <= 0.05)
  
  tfbm_all = 
    ttT %>%
    infer_db_adapted(ttT_sub = ttT_sub)
  
  tfbm_immue = 
    ttT %>% 
    infer_db_adapted(ttT_sub = ttT_sub, which_tfbms = immune_tfbms)
  
  tfbm = list(tfbm_all = ttT %>%
                infer_db_adapted(ttT_sub = ttT_sub),
              tfbm_immue = ttT %>% 
                infer_db_adapted(ttT_sub = ttT_sub, which_tfbms = immune_tfbms))
  tfbm_telis = tfbm %>% map(. %>%
                              pluck("telis", "par") %>% 
                              as.data.frame() %>% 
                              rownames_to_column(var = "tfbm"))
  tfbm_reg = tfbm %>% map(. %>%
                            pluck("m", "logFC") %>% 
                            as.data.frame()) 
  tfbm = map2(tfbm_telis, tfbm_reg, ~ left_join(.x, .y, by = c("tfbm" = "m_uni.term")) %>% select(tfbm, p_under, p_over, m_uni.p.value, m_cov.p.value ) )
  
     return(list(ttT = ttT, tfbm = tfbm))
                 # , gene_set_test = gene_set_test))
}