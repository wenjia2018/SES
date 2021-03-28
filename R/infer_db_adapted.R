# GENERALIZE A FUNCTION FROM dbr, IN ORDER THAT TFBM ANALYSIS ACCOMODATE FACTORS
infer_db_adapted = 
  function (ttT, ttT_sub = NULL, which_matrix = NULL, which_tfbms = NULL, 
            explicit_zeros = TRUE, perm_telis = FALSE, n_sim = 1e+05){
    if (is.null(which_matrix)) 
      which_matrix = utr1
    # change get_matrix to get_matrix_modified to solve the problem of matrix subsetting with only 1 row, which becomes an array
    R = get_matrix_modified(ttT = ttT, which_matrix = which_matrix, which_tfbms = which_tfbms, 
                         explicit_zeros = explicit_zeros)
    ttT = dbr:::append_matrix(ttT = ttT, ttT_sub = ttT_sub, R = R)
    telis = NULL
    if (!is.null(ttT_sub)) 
      telis <- dbr:::get_telis(R = R, ttT_sub = ttT_sub, n_sim = n_sim, 
                               perm_telis = perm_telis)
    m = ttT %>% dplyr::select(logFC) %>% purrr::map(regress_db, 
                                                      X = dplyr::select(ttT, colnames(R)))
    return(out = list(ttT = ttT, telis = telis, m = m))
  }