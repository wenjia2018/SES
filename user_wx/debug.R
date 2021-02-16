debugonce(fit_bespoke)
example_bespoke <- args_eqtl %>% mutate(out = pmap(list(gene_set_name = table1, p_eqtl = p_eqtl), safely(fit_bespoke)))


debugonce(model_fit)

out$m8_fdr = fit_m8(controls, treatment, gene_set) %>% extract_m8_fdr()
debugonce(mediate_multiple)
out$m8_fdr$mediation = mediate_multiple(controls, treatment, gene_set = out$m8_fdr$sig_genes)


 debugonce(mediate)
 pmap(list(args_m97$mediator,
           args_m97$gene_set,
           args_m97$controls,
           args_m97$treatment),
      safely(mediate)) %>% 
   set_names(args_m97$names)
 
 debugonce(fit_m99)
 fit_m99(datt_m, gene_set, mediator) %>% extract_m99()   
 