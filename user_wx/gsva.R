dat_sub = dat
pData(dat_sub) = pData(dat) %>% dplyr::select(ses_sss_composite)
pData(dat_sub) = pData(dat) %>% dplyr::select(income_hh_ff5)

rownames(pData(dat_sub)) = dat_sub@assayData[["exprs"]] %>% colnames
gsva_es <- GSVA::gsva(dat_sub, signatures$outcome_set[table1[1:2]])
# Estimating GSVA scores for 13 gene sets.
# Computing observed enrichment scores
# Estimating ECDFs with Gaussian kernels
# Using parallel with 56 cores

gsva_es %>% saveRDS("./user_wx/edu_gsva.rds")


tmp = gsva_es@assayData[["exprs"]]
scores = tmp %>% t %>% cbind(gsva_es@phenoData@data)




# scoreplot = function(set_name){
#   ggplot(scores, aes(x = edu_max, y = set_name, color = edu_max)) +
#     geom_boxplot(alpha=0.2) + 
#     xlab("edu")
#   
# }

scoreplot = function(of_in, set_name){
  ggpubr::ggboxplot(data = scores, x = of_in , y = set_name)
}


set_name = scores %>% colnames() %>% str_subset("mRNA")

pdf(file = "./user_wx/edu_scores_plot.pdf")

map(set_name, ~ scoreplot("edu_max", .x))

dev.off() 
 
 
