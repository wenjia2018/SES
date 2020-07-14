model_fit = 
  function(gene_set_name, treatment, controls, funcs, out = NULL){  
    # e.g. ctra_mRNA
    # not a pure function (a and gene_set are mutable)
    
    ############################################################
    # CUT EXECUTION SHORT IF THE OUTCOME IS THE ENTIRE GENOME  
    ############################################################
    
    if(gene_set_name == "whole_genome_and_tfbm") return(de_and_tfbm(treatment, controls)) 
    
    ############################################################
    # OTHERWISE (FOR SMALLER GENE SETS OF INTEREST)
    ############################################################
    
    print("*******************")
    print(treatment)
    print(gene_set_name)
    print("*******************") 
    
    gene_set = pluck(signatures, "outcome_set", gene_set_name)
    
    datt =
      prepro(gene_set, treatment, controls) %>% 
      rename(treatment = treatment) %>% 
      remove_diseased_subjects_from_datt(gene_set_name, controls)
    
    if(is.element("m1", funcs))  out$m1  = fit_m1(datt, gene_set) %>% extract_m1()
    if(is.element("m2", funcs))  out$m2  = fit_m2(datt, gene_set) %>% extract_m2()
    if(is.element("m3", funcs))  out$m3  = fit_m3(datt, gene_set) %>% extract_m3()
    if(is.element("m4", funcs))  out$m4  = fit_m4(datt, gene_set) %>% extract_m4()
    if(is.element("m5", funcs))  out$m5  = fit_m5(controls, treatment, gene_set) %>% extract_m5()
    if(any(str_detect(funcs, "m6"))) {
      
      out$m6_nn = fit_m6(datt, gene_set, "none"   ) %>% extract_m6() 
      out$m6_vx = fit_m6(datt, gene_set, "varimax") %>% extract_m6()
      out$m6_ob = fit_m6(datt, gene_set, "oblimin") %>% extract_m6()
      
      if(0){ 
        # ALTERNATIVELY, MORE GENERAL BUT PROBABLY TOO GENERAL UNLESS WE CONSIDER MANY ROTATIONS
        out = out %>% 
          append(c(m6_nn = "none", m6_vx = "varimax", m6_ob = "oblimin") %>% 
                   map(fit_m6, datt = datt, gene_set = gene_set) %>% 
                   map(extract_m6))
      }
      
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
    if(is.element("m98", funcs)) out$m98 = fit_m98(datt, gene_set) %>% extract_m98()
    source("R/utils.R", local = TRUE)
    if(is.element("m97", funcs)){
      args_m97 = crossing(mediators, gene_set) %>% rename(mediator = mediators) %>% 
        mutate(names = str_c(mediator, gene_set, sep = "_"))
      out$m97 = map2(args_m97$mediator, args_m97$gene_set, safely(mediate)) %>% set_names(args_m97$names)
    }
    if(is.element("m99", funcs)) out$m99 =  mediators %>% set_names() %>%map(safely(mediate), gene_set = gene_set) 
    if(0){ 
      # More elegant but requires more work to unify arguments across functions
      eval_me = function(x) compose(!!!str_subset(ls(.GlobalEnv), x)) 
      c("m1", "m3", "m4") %>% # generally should be funcs arg above
        set_names() %>% 
        map(eval_me) %>%
        map(exec, gene_set = gene_set, datt = datt)
    } 
    
    out = out %>% map(list) %>% as_tibble()
    return(out = out)
  }
