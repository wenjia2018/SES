
#' ---
#' title: Omnibus test results for skin color
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->


#' ###  regression of gene on several skin color indicators(5 levels, continuous)
#' * gene ~ skincolor + controls
#' * p values are FDR corrected within each signature
#' * Treatment:
#' *   1. skincolor as continuous 
#' * 2. Skincolor as 5 levels category variable(Black, Dark, Medium, Light, Light as reference, White is removed) 

#' * Control:
#' *   1. Basic: ses paper control removing Race ethnicity
#' * 2. All : basic + ancestry PC1:PC20
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
res_fun = function(example0_with1k, example0_without1k) {
  ex0_without1k <- example0_without1k %>%
    filter(gene_set_name!="inflam1k_mRNA" ) %>% 
    hoist(out, m = list("result", "m12_fdr", 1, "other", "m")) %>% 
    filter(!map_lgl(m, ~is.null(.x))) %>% 
    mutate(p = m %>% map(~ dplyr::slice(.x, which.min(adj.p.within))) %>% map_dbl(~ .x %>% pluck("adj.p.within"))) %>% 
    dplyr::select(treatment, gene_set_name, controls, p) %>% 
    filter(p<0.05) %>% 
    # dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
    # dplyr::filter(treatment %in% c("edu_p", "income_pp1_log", "SEI_max_p_w12", "ses_composite_pp1" )) %>% 
    mutate(pval=case_when(p<0.0001 ~ 0.0001,
                          p<0.001 ~ 0.001,
                          p<0.01 ~ 0.01,
                          p<0.05 ~ 0.05,
                          p>0.05 ~ 100),
           pval2=case_when(p<0.0001 ~ 100000,
                           p<0.001 ~ 25000,
                           p<0.01 ~ 15000,
                           p<0.05 ~ 10000,
                           p>0.05 ~ 0.0000001)
    )%>%
    mutate(
      # gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
      #      gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
      #      gene_set_name = gene_set_name %>% replace(gene_set_name == "Cvd", "CVD"),
      #      gene_set_name = gene_set_name %>% replace(gene_set_name == "Copd", "COPD"),
      #      gene_set_name = gene_set_name %>% replace(gene_set_name == "Ckd", "CKD"),
      #      gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI"),
      "1KI Genes" = "Without 1KI Genes"%>% as.factor())
  
  
  ex0_with1k <- example0_with1k %>%
    hoist(out, m = list("result", "m12_fdr", 1, "other", "m")) %>% 
    filter(!map_lgl(m, ~is.null(.x))) %>% 
    mutate(p = m %>% map(~ dplyr::slice(.x, which.min(adj.p.within))) %>% map_dbl(~ .x %>% pluck("adj.p.within"))) %>% 
    dplyr::select(treatment, gene_set_name, controls, p) %>% 
    filter(p<0.05) %>% 
    # dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
    # dplyr::filter(treatment %in% c("edu_p", "income_pp1_log", "SEI_max_p_w12", "ses_composite_pp1" )) %>% 
    mutate(pval=case_when(p<0.0001 ~ 0.0001,
                          p<0.001 ~ 0.001,
                          p<0.01 ~ 0.01,
                          p<0.05 ~ 0.05,
                          p>0.05 ~ 100),
           pval2=case_when(p<0.0001 ~ 100000,
                           p<0.001 ~ 25000,
                           p<0.01 ~ 15000,
                           p<0.05 ~ 10000,
                           p>0.05 ~ 0.0000001)
    )%>%
    mutate(
      # gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
      #      gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
      #      gene_set_name = gene_set_name %>% replace(gene_set_name == "Cvd", "CVD"),
      #      gene_set_name = gene_set_name %>% replace(gene_set_name == "Copd", "COPD"),
      #      gene_set_name = gene_set_name %>% replace(gene_set_name == "Ckd", "CKD"),
      #      gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI"),
      "1KI Genes" = "With 1KI Genes"%>% as.factor())
  
  bind_rows(ex0_with1k, ex0_without1k)
}


# nonhispanic black sample
example0_without1k <- readRDS("~/ses-1/user_wx/m12_without1k_scbinary_NonBstrata.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_scbinary_NonBstrata.rds")

res1 = res_fun(example0_with1k, example0_without1k)

res1 %>% 
  dplyr::select(-pval, - pval2)%>% 
  kableExtra::kable() %>% kableExtra::kable_styling()

example0_without1k <- readRDS("~/ses-1/user_wx/m12_without1k_sc4levels_NonBstrata.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_sc4levels_NonBstrata.rds")

res2 = res_fun(example0_with1k, example0_without1k)
res2 %>%  dplyr::select(-pval, - pval2)%>% 
  kableExtra::kable() %>% kableExtra::kable_styling()

example0_without1k <- readRDS("~/ses-1/user_wx/m12_without1k_sccont_NonBstrata.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_sccont_NonBstrata.rds")

res3 = res_fun(example0_with1k, example0_without1k)
res3 %>% 
  kableExtra::kable() %>% kableExtra::kable_styling()

res = bind_rows(res1,res2, res3) %>% 
  filter(gene_set_name != "inflame_mRNA") %>% 
  dplyr::select(-pval, - pval2)
# %>% 
#   filter(controls=="basic", gene_set_name %in% aging) %>% 
#   dplyr::select("treatment", "gene_set_name", "p", "1KI Genes")
res %>% kableExtra::kable() %>% kableExtra::kable_styling()

# nonhispanic black sample pc1-pc4
example0_without1k <- readRDS("~/ses-1/user_wx/m12_without1k_scbinary_NonBstrata_PC4.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_scbinary_NonBstrata_PC4.rds")

res1 = res_fun(example0_with1k, example0_without1k)

res1 %>% 
  dplyr::mutate(controls = "ancestry PC1-4") %>% 
  kableExtra::kable() %>% kableExtra::kable_styling()


example0_without1k <- readRDS("~/ses-1/user_wx/m12_without1k_sc4levels_NonBstrata_PC4.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_sc4levels_NonBstrata_PC4.rds")

res2 = res_fun(example0_with1k, example0_without1k)

res2 %>% 
  dplyr::mutate(controls = "ancestry PC1-4") %>% 
  kableExtra::kable() %>% kableExtra::kable_styling()


example0_without1k <- readRDS("~/ses-1/user_wx/m12_without1k_sccont_NonBstrata_PC4.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_sccont_NonBstrata_PC4.rds")

res3 = res_fun(example0_with1k, example0_without1k)


res3 %>% 
  dplyr::mutate(controls = "ancestry PC1-4") %>% 
  kableExtra::kable() %>% kableExtra::kable_styling()

res = bind_rows(res1,res2, res3) %>% 
  filter(gene_set_name != "inflame_mRNA") %>% 
  dplyr::select(-pval, - pval2)
# %>% 
#   filter(controls=="basic", gene_set_name %in% aging) %>% 
#   dplyr::select("treatment", "gene_set_name", "p", "1KI Genes")
res %>% 
  dplyr::mutate(controls = "ancestry PC1-4") %>% 
  kableExtra::kable() %>% kableExtra::kable_styling()
