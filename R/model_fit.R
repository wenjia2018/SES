model_fit = 
  function(gene_set_name, treatment, controls, funcs, out = NULL){  
    # e.g. ctra_mRNA
    # not a pure function (a and gene_set are mutable)
    
    ############################################################
    # CUT EXECUTION SHORT IF THE OUTCOME IS THE ENTIRE GENOME  
    ############################################################
    
    if(gene_set_name == "whole_genome_and_tfbm") return(de_and_tfbm(treatment, controls)) 
    if(funcs == "m96") return(celltype_cibersort(treatment, controls)) 
    ############################################################
    # OTHERWISE (FOR SMALLER GENE SETS OF INTEREST)
    ############################################################
    
    print("*******************")
    print(treatment)
    print(gene_set_name)
    print(controls)
    print("*******************") 
    
    gene_set = pluck(signatures, "outcome_set", gene_set_name)
    
    datt =
      prepro(gene_set, treatment, controls) %>% 
      rename(treatment = treatment) %>% 
      remove_diseased_subjects_from_datt(gene_set_name, controls)
    
    if(is.element("m1", funcs)){
      
      if(str_detect(gene_set_name, "inflam1k_mRNA")){ 
        
        out$m1 = NA # estimation fails because gene set is too big for lm
        
      } else {
        
        out$m1  = fit_m1(datt, gene_set) %>% extract_m1()      
        
      }
      
    }
    
    if(is.element("m2", funcs))  out$m2  = fit_m2(datt, gene_set) %>% extract_m2()
    if(is.element("m3", funcs))  out$m3  = fit_m3(datt, gene_set) %>% extract_m3()
    if(is.element("m4", funcs))  out$m4  = fit_m4(datt, gene_set) %>% extract_m4()
    if(is.element("m5", funcs))  out$m5  = fit_m5(controls, treatment, gene_set) %>% extract_m5()
    if(is.element("m5", funcs))  out$m5b = fit_m5(controls, treatment, gene_set) %>% extract_m5b()
    
    if(any(str_detect(funcs, "m6"))) {
      
      out$m6_nn = fit_m6(datt, gene_set, "none"   ) %>% extract_m6() 
      out$m6_vx = fit_m6(datt, gene_set, "varimax") %>% extract_m6()
      out$m6_ob = fit_m6(datt, gene_set, "oblimin") %>% extract_m6()
      
    }
    
    if(any(str_detect(funcs, "m7"))) {
      
      out$m7_nn  = fit_m7(datt, gene_set, "none"   ) %>% extract_m7()
      out$m7_vx  = fit_m7(datt, gene_set, "varimax") %>% extract_m7()
      out$m7_ob  = fit_m7(datt, gene_set, "oblimin") %>% extract_m7()
      
      # for each of the results just calculated above, append the genes loading high in this dimension
      out$m7_nn$other <- get_well_loaded_genes(datt, gene_set, "none")
      out$m7_vx$other <- get_well_loaded_genes(datt, gene_set, "varimax")
      out$m7_ob$other <- get_well_loaded_genes(datt, gene_set, "oblimin")
      
    } 
    
    if(any(str_detect(funcs, "m8"))){
      
      out$m8_fdr  = fit_m8(controls, treatment, gene_set) %>% extract_m8_fdr()
      out$m8_fwer = fit_m8(controls, treatment, gene_set) %>% extract_m8_fwer()
      
    }
   
    source("R/utils.R", local = TRUE)
    if(is.element("m97", funcs)) out$m97 = mediate_multiple(gene_set)
    if(is.element("m99", funcs)) out$m99 = mediators %>% set_names() %>% map(safely(mediate), gene_set = gene_set) 

    out = out %>% map(list) %>% as_tibble()
    return(out = out)
  }
