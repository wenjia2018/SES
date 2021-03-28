get_matrix_modified <- function(ttT = ttT, which_matrix = which_matrix, which_tfbms = which_tfbms, explicit_zeros = explicit_zeros){
  
  # LOAD TFBM MATRIX
  R = which_matrix # the gene x motif binding "R"egulation matrix
  
  # ATTEMPT TO REMOVE UNINFORMATIVE TFBMS. REMOVE ROWS WITH NO VARIATION.
  # COULD ALSO ADDRESS COLINEARITY HERE.
  R = R[, matrixStats::colSds(R) >= 0.1, drop = FALSE]
  
  # PERHAPS SUBSET COLS OF R
  # CASE 1: ONLY ONE NON NULL TFBM SPECIFIED BEWARE TYPE CHANGE
  # CASE 2: ANY OTHER NON-NULL, NON-SINGLETON PROPER SUBSET OF COLS OF R
  # CASE 3: EXAMINE ALL TFBMS
  in_frame  = rownames(R) %in% ttT$gene # which tfbm genes are in the sampling frame
  if((!is.null(which_tfbms)) & (length(which_tfbms) == 1)){
    Rt  = R[in_frame, which_tfbms, drop = FALSE]
    R   = matrix(Rt, ncol = 1) %>% `rownames<-`(names(Rt)) %>% `colnames<-`(which_tfbms)
  } else if((!is.null(which_tfbms)) & (length(which_tfbms) >= 2)) {
    R = R[in_frame, which_tfbms, drop = FALSE]        # restrict R to those genes which could in principle be differentially expressed
  } else if(is.null(which_tfbms)) {
    which_tfbms = colnames(R)
    R = R[in_frame, which_tfbms, drop = FALSE]
  }
  
  # The original base tfbm matrix only includes genes with at least one tfbm for at least one regulator.
  # Genes not mentioned in the base tfbm matrix implicitly have a count of zero for every tfbm.
  # The preceding code additionally excludes any genes in the base tfbm matrix but not in the sampling frame (ttT$gene).
  # The next option introduces explicit zeros to R, for all genes in the sampling frame but with (implicitly) zero tfbm count for any motif.
  if(explicit_zeros){
    print("Adding explicit zero counts to the tfbm matrix for genes in sampling frame, but not in tfbm matrix")
    not_in_R = setdiff(ttT$gene, rownames(R))
    nR = matrix(0,length(not_in_R), dim(R)[2])
    rownames(nR) = not_in_R
    R = rbind(R, nR)
  }
  return(R = R)
}
