de_and_tfbm_normalization_TMM = function(controls, treatment, gene_set_name, de_only = TRUE) {

  if(remove_diseased_subjects){
    keep_nondisease = 
      switch(
        gene_set_name,
        diabetes_mRNA = dat$diabetes != 1,
        CVD_mRNA = dat$heartatk != 1,
        Asthma_mRNA = dat$H5ID6F != 1 & dat$H5ID6FM != 1,
        Hypertension_mRNA = dat$H5ID6C != 1 & dat$H5ID6CM != 1,
        Aortic_Aneurysm_mRNA = dat$H5ID6Q != 1 & dat$H5ID6QM != 1,
        Melanoma_mRNA = dat$H5ID6A != 1 & dat$H5ID6AM != 1,
        rep(TRUE, dim(dat)[2])
      )
    # 
    keep_nondisease = ifelse(is.na(keep_nondisease), FALSE, keep_nondisease)
    # 
    dat = dat[, keep_nondisease]

  }
  rownames(dat@phenoData@data) = dat@phenoData@data$AID
  raw = raw[, sampleNames(dat)]

  # ph = Biobase::pData(raw)
  ex = Biobase::exprs(raw)
  
  ## Get pheno data for all subjects
  ph = Biobase::pData(dat)
  rownames(ph) = ph$AID
  all.equal(colnames(ex), ph$AID)
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
  
  # Specify whole-genome regression of rna on design
  data <- ph %>% dplyr::select(all_of(controls), all_of(treatment))
  keep = data %>% complete.cases()
  data = data[keep, ]
  batch_data = ph$batch[keep]
  # FIND GENES WITHOUT ENOUGH VARIATION
  # min.cout min.pro setting from Ravi's domain knowledge
  # Ravi observation population not experiment
  e_genes = edgeR::filterByExpr(ex[, keep], min.count = 10, min.pro = 0.10)
  e_genes = names(e_genes[e_genes == TRUE])
  # REMOVE THOSE GENES and missing subjects
  ex = ex[e_genes, keep]
  

  dge = edgeR::DGEList(counts = ex)
  
  dge_tmm = edgeR::calcNormFactors(dge, method = "TMM")
  
  ## Design matrix

  X = model.matrix(~ ., data = data)
  # drop columns to make sure X is full rank
  X = ordinal::drop.coef(X)
  # MatrixModels::model.Matrix can drop unused levels, but has error when the variable names contains some special characters 
  # https://stackoverflow.com/questions/44114391/error-in-sparse-model-matrix
  # X = MatrixModels::model.Matrix(~., data = X, sparse = TRUE, drop.unused.levels = TRUE)
  
  
  # Voom transformation
  v_tmm = voom(dge_tmm, X)
  # Batch correct
  v_tmm$E = sva::ComBat(v_tmm$E, batch_data)
  
  # Estimate DE using standard limmma/edger pipeline. 
  ttT_raw <-
    lmFit(v_tmm$E, X) %>%
    eBayes(trend = T)
 

  
   ttT_boot= NULL
  

  if(boot==TRUE) {
    pseudo_E = v_tmm$E - ttT_raw$coefficients[,treatment,drop=FALSE] %*% t(ttT_raw$design[,treatment,drop=FALSE])
    # library(foreach)
    # myCluster <- parallel::makeCluster(30)
    # doParallel::registerDoParallel(myCluster)
    
    # foreach(i = 1:N) %dopar% {
    for(i in 1:N){
      index = sample(1:dim(X)[1], replace = TRUE)
      X_boot = ordinal::drop.coef(X[index, ])
      ind = X_boot %>% colnames() %>% str_detect(treatment) %>% sum
      while(ind < 1){
        index = sample(1:dim(X)[1], replace = TRUE)
        X_boot = ordinal::drop.coef(X[index, ])
        ind = X_boot %>% colnames() %>% str_detect(treatment) %>% sum
      }
      treatment_ind = X_boot %>% colnames() %>% str_detect(stringr::fixed(treatment)) %>% which()
      ttT_boot[[i]] <-
        lmFit(pseudo_E[, index], X_boot) %>%
        eBayes(trend = T) %>%
        tidy_topTable(of_in = treatment_ind, confint = TRUE)
    } 
    # parallel::stopCluster(myCluster)
  }
  
  # needed for limma::topTables() to work with factors
  treatment = X %>% colnames() %>% str_detect(stringr::fixed(treatment)) %>% which()
  

  #################
  ttT <-
    ttT_raw %>%
    tidy_topTable(of_in = treatment, confint = TRUE)
  #################  
  # extract F and calculate R square
  # number of predictors in the model
  k = dim(ttT_raw$coefficients)[2] - 1
  # residual degrees of freedom
  df_residual = ttT_raw$df.residual %>% unique
  if(0) {
    # r square from the full model(including treatment and controls)
    rsquare_full = ttT_raw$F %>% map(~ .x/(.x + df_residual / k)) %>% set_names(ttT_raw$sigma %>% names)
    vars = ttT_raw$coefficients %>% colnames %>% `[`(- treatment)
    ttT_control = ttT_raw %>% tidy_topTable(of_in = vars, confint = TRUE)
    rsquare_ohnetreatment = ttT_control$F %>% map(~ .x/(.x + (df_residual+1) / (k-1))) %>% set_names(ttT_control$gene)
    
    rsquare_full = rsquare_full[order(names(rsquare_full))]
    rsquare_ohnetreatment = rsquare_ohnetreatment[order(names(rsquare_ohnetreatment))]
    
    rsq_diff = map2(rsquare_full, rsquare_ohnetreatment, `-`)
    f_statistic = map2(rsquare_full, rsquare_ohnetreatment, ~ (.x - .y) * (df_residual)/(1 - .x)) 
  }
  # F value
  f_statistic = ttT$t^2 
  F_pval = pf(f_statistic %>% unlist, 1, df_residual, lower.tail = FALSE) %>% p.adjust("fdr") %>% setNames(ttT$gene)
  if(de_only) return(list(ttT = ttT, ttT_boot = ttT_boot, F_pval = F_pval))
}