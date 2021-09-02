model_fit = 
  function(gene_set_name, treatment, controls, funcs, out = NULL){  
    # e.g. ctra_mRNA
    # not a pure function (a and gene_set are mutable)
    
    ############################################################
    # CUT EXECUTION SHORT IF THE OUTCOME IS THE ENTIRE GENOME  
    ############################################################
    # regression in the whole gemone controlling for ancestryPC of each signature
    if(gene_set_name == "whole_genome_and_tfbm") {
      if(normalization_bydesign == TRUE) return(de_and_tfbm_normalization_TMM(controls, treatment, gene_set_name)) 
      else return(de_and_tfbm(treatment, controls)) 
    }

    # regression in the whole gemone controlling for ancestryPC of whole genome
    if(gene_set_name == "whole_genome") return(fit_m10(treatment, controls, gene_set, ttT_within_genesets = FALSE, tfbm_genesets = FALSE, wholegenome = TRUE))
    if(funcs == "m96") return(celltype_cibersort(treatment, controls)) 
    # controls: NULL or controls + ses predictor
    if(funcs == "m95") return(AER_IV(gene_set_name, treatment, IV, controls)) 
    
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
      dplyr::rename(treatment = treatment) %>% 
      remove_diseased_subjects_from_datt(gene_set_name, controls)
    
    if(is.element("m1", funcs)){
      
      if(str_detect(gene_set_name, "inflam1k_mRNA|aging_cluster_complement_mRNA|aging_mRNA|aging_down_mRNA|aging_up_mRNA")){ 
        
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
    source("R/utils.R", local = TRUE)
    if(any(str_detect(funcs, "m7"))) {
      
      # out$m7_nn  = fit_m7(datt, gene_set, "none"   ) %>% extract_m7()
      # out$m7_vx  = fit_m7(datt, gene_set, "varimax") %>% extract_m7()
      out$m7_ob  = fit_m7(datt, gene_set, "oblimin") %>% extract_m7()
      # for each of the results just calculated above, append the genes loading high in this dimension
      # out$m7_nn$other$well_loaded <- get_well_loaded_genes(datt, gene_set, "none")
      # out$m7_vx$other$well_loaded <- get_well_loaded_genes(datt, gene_set, "varimax")
      out$m7_ob$other$well_loaded <- get_well_loaded_genes(datt, gene_set, "oblimin")
      
      # for each of the pca, do mediational analysis for each pc
      # out$m7_nn$mediation = mediators %>% set_names() %>% map(safely(mediate_pca), gene_set = gene_set, rotate = "none", out$m7_nn)
      # out$m7_vx$mediation = mediators %>% set_names() %>% map(safely(mediate_pca), gene_set = gene_set, rotate = "varimax", out$m7_vx)
      out$m7_ob$mediation = mediators %>% set_names() %>% map(safely(mediate_pca), gene_set = gene_set, rotate = "oblimin", out$m7_ob)
      
      out$m7_ob$dda = mediators %>% set_names() %>% map(safely(DDA), gene_set = gene_set, rotate = "oblimin", out$m7_ob)
    } 
    
    if(any(str_detect(funcs, "m8"))){
      
      out$m8_fdr = fit_m8(controls, treatment, gene_set) %>% extract_m8_fdr()
      
      if(mediation_mean & length(out$m8_fdr$sig_genes) > 0){
        # out$m8_fdr$mediation_single = mediate_multiple(controls, treatment, gene_set = out$m8_fdr$sig_genes)
        out$m8_fdr$mediation_mean = mediators %>% set_names() %>% map(safely(mediate), gene_set = out$m8_fdr$sig_genes, controls, treatment) 
      } else{
        # out$m8_fdr$mediation_single = NULL
        out$m8_fdr$mediation_mean = NULL
      }
      
      if(mediation_each_gene) {
        # out$m8_fdr$mediation_single = mediate_multiple(controls, treatment, gene_set = gene_set)
        plan(multicore, workers = 40)
        out$m8_fdr$mediation_single = furrr::future_map(gene_set, ~ mediate_multiple(controls, treatment, .x))
      }

      # out$m8_fwer = fit_m8(controls, treatment, gene_set) %>% extract_m8_fwer()
      
    }
    # do gene by gene mediation for all DE sig genes
    # if(any(str_detect(funcs, "m9"))){
    #   
    #   out$m9_fdr= de_and_tfbm(treatment, controls, de_only = TRUE) %>% extract_m8_fdr()
    #   
    #   if(length(out$m9_fdr$sig_genes) > 0){
    #     out$m9_fdr$mediation_single = mediate_multiple(controls, treatment, gene_set = out$m9_fdr$sig_genes)
    #     out$m9_fdr$mediation_mean = mediators %>% set_names() %>% map(safely(mediate), gene_set = out$m9_fdr$sig_genes, controls, treatment)
    #   } else{
    #     out$m9_fdr$mediation_single = NULL
    #     out$m9_fdr$mediation_mean = NULL
    #   }
    # }
    if(any(str_detect(funcs, "m12"))){
      
      out$m12_fdr = fit_m12(controls, treatment, gene_set_name) %>% extract_m8_fdr()
      
    }
    if(is.element("m10", funcs)) out$m10 = fit_m10(treatment, controls, gene_set, ttT_within_genesets = TRUE, tfbm_genesets = FALSE, wholegenome = FALSE)
    if(is.element("m11", funcs)) out$m11 = fit_m10(treatment, controls, gene_set, ttT_within_genesets = FALSE, tfbm_genesets = TRUE, wholegenome =FALSE)
    if(is.element("m97", funcs)) out$m97 = mediate_multiple(controls, treatment, gene_set)
    if(is.element("m99", funcs)) out$m99 = mediators %>% set_names() %>% map(safely(mediate), gene_set = gene_set, controls, treatment) 

    out = out %>% map(list) %>% as_tibble()
    return(out = out)
  }
