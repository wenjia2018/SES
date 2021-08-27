de_and_tfbm_normalization = function(treatment, controls, de_only = TRUE) {
  
  ph = Biobase::pData(raw)
  ex = Biobase::exprs(raw)

  ## Get pheno data for all subjects
  ph = pData(dat)
  all.equal(colnames(ex), ph$AID)
  
  # Specify whole-genome regression of rna on design
  X <- ph %>% dplyr::select(all_of(controls), all_of(treatment))
  keep = X %>% complete.cases()
  X = X[keep, ]
  # FIND GENES WITHOUT ENOUGH VARIATION
  # min.cout min.pro setting from Ravi's domain knowledge
  # Ravi observation population not experiment
  e_genes = edgeR::filterByExpr(ex[, keep], min.count = 10, min.pro = 0.10)
  e_genes = names(e_genes[e_genes == TRUE])
  # REMOVE THOSE GENES and missing subjects
  ex = ex[e_genes, keep]
  
  ## Convert gene names to HGNC
  if(full_reproducibility <- FALSE) {
    library('biomaRt')
    glist  = rownames(dat) # gene names
    mart   = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="useast.ensembl.org")
    G_list = getBM(filters= "ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values = glist, mart= mart)
    saveRDS(G_list, str_c(data_input,"genename_allbatches.rds"))
  }
  G_list = readRDS("/home/share/data_input/genename_allbatches.rds")
  # keep duplicated ensemble id suggested by Ravi
  G_list = G_list  %>% 
    dplyr::filter(hgnc_symbol != "")
    # dplyr::filter(!duplicated(ensembl_gene_id) & !duplicated(hgnc_symbol)) %>%
  
  
  # https://stackoverflow.com/questions/7547597/dictionary-style-replace-multiple-items/25790005#25790005
  
  dic = setNames(G_list$hgnc_symbol, G_list$ensembl_gene_id)
  rownames(ex) <- dic[rownames(ex)]
  nas = is.na(rownames(ex))
  ex = ex[!is.na(rownames(ex)), ]
  dge = edgeR::DGEList(counts = ex)
  
  dge_rle = edgeR::calcNormFactors(dge, method = "RLE")
  
  ## Design matrix
  X = model.matrix(~ ., data = X)
  # drop columns to make sure X is full rank
  X = ordinal::drop.coef(X)
  
  # Voom transformation
  v_rle = voom(dge_rle, X)
  
  # MatrixModels::model.Matrix can drop unused levels, but has error when the variable names contains some special characters 
  # https://stackoverflow.com/questions/44114391/error-in-sparse-model-matrix
  # X = MatrixModels::model.Matrix(~., data = X, sparse = TRUE, drop.unused.levels = TRUE)
  # needed for limma::topTables() to work with factors
  treatment = X %>% colnames() %>% str_detect(stringr::fixed(treatment)) %>% which()
  # treatment = X %>% colnames() %>% str_detect(treatment) %>% which()
  
  # Estimate DE using standard limmma/edger pipeline. 
  ttT_raw <-
    lmFit(v_rle, X) %>%
    eBayes
  #################  
  # extract F 
  # ttT_raw$F
  #################
  ttT <-
    ttT_raw %>%
    tidy_topTable(of_in = treatment, confint = TRUE)
  if(de_only) return(list(ttT = ttT))
}