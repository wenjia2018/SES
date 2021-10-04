#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA

#' ### abbreviation
#' * ACME : average causal mediated effect
#' * PM : proportion mediated
#' * ADE : average direct effect
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA

library(tidyverse)
example <- readRDS("~/ses-1/user_wx/example_tmm_m17_withinflame_1000boots.rds")

# replace some cases with 0.008 threshold.
# recode edu as 0 1 and run the analysis again
multimediate.sescomposite <- readRDS("~/ses-1/user_wx/multimediate.sescomposite.rds")
multimediate.income <- readRDS("~/ses-1/user_wx/multimediate.income.rds")
multimediate.sei <- readRDS("~/ses-1/user_wx/multimediate.sei.rds")
multimediate.edu <- readRDS("~/ses-1/user_wx/multimediate.edu.rds")
multimediate.edu.0.008 <- readRDS("~/ses-1/user_wx/multimediate.edu.0.008.rds")


example$out[which(example$treatment =="income_hh_ff5"&example$gene_set_name == "CKD_mRNA")] =
  multimediate.income$out[which(multimediate.income$treatment =="income_hh_ff5"&multimediate.income$gene_set_name == "CKD_mRNA")]


example$out[which(example$treatment =="SEI_ff5"&example$gene_set_name %in% c("COPD_mRNA", "Asthma_mRNA", "CKD_mRNA"))] =
  multimediate.sei$out[which(multimediate.sei$treatment =="SEI_ff5"&multimediate.sei$gene_set_name %in% c("COPD_mRNA", "Asthma_mRNA", "CKD_mRNA"))]

example$out[which(example$treatment =="ses_sss_composite"&example$gene_set_name == "CVD_mRNA")] =
  multimediate.sescomposite$out[which(multimediate.sescomposite$treatment =="ses_sss_composite"&multimediate.sescomposite$gene_set_name == "CVD_mRNA")]

example$out[which(example$treatment =="edu_max"&example$gene_set_name %in% c("CKD_mRNA", "COPD_mRNA", "CVD_mRNA","Depression_mRNA", "Hypertension_mRNA", "inflam1k_mRNA","Rheumatoid_Arthritis_mRNA"))] =
  multimediate.edu$out[which(multimediate.edu$treatment =="edu_max"&multimediate.edu$gene_set_name %in% c("CKD_mRNA", "COPD_mRNA", "CVD_mRNA","Depression_mRNA", "Hypertension_mRNA", "inflam1k_mRNA","Rheumatoid_Arthritis_mRNA"))]

example$out[which(example$treatment =="edu_max"&example$gene_set_name %in% c("Hypertension_mRNA", "Rheumatoid_Arthritis_mRNA"))] =
  multimediate.edu.0.008$out[which(multimediate.edu.0.008$treatment =="edu_max"&multimediate.edu.0.008$gene_set_name %in% c("Hypertension_mRNA", "Rheumatoid_Arthritis_mRNA"))]




threshold = 0.05
saved = 
  example %>%  
  hoist(out, p = list("result", "m17_ob", 1, "p")) %>%
  unnest_longer(p) %>%
  dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>%
  mutate(p = p.adjust(p, method = "fdr")) %>%
  dplyr::filter(p < threshold) %>%
  mutate(pcmin = p_id %>% str_remove("d") %>% as.numeric()) %>%
  group_by(treatment, gene_set_name) %>%
  mutate(p_no = n()) %>%
  slice(which.min(pcmin)) %>%
  ungroup %>%
  hoist(out, proporp = list ("result", "m17_ob", 1, "mediation")) %>% 
  mutate(multimed = map2(.$p_id, .$proporp, ~ .y[[.x]]))


ses = saved %>% mutate(names = str_c(treatment,gene_set_name))
#' ### multiple mediation results 
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
for (i in 1 : dim(ses)[1]) {
  
  cat(" ############################################################","\n",
      "treatment is", ses[i, ]$treatment, "gene_set_name is", ses[i, ]$gene_set_name,"\n",
      "############################################################","\n")
  cat("multiple mediators","\n","\n")
  
  print(ses[i,]$multimed)
  cat("\n","\n")
  
}


