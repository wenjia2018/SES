library(tidyverse)
threshold = 0.05
# rle
# example0_with1k <- readRDS("~/ses-1/user_wx/example_RLE_pca_nomed_withinflame.rds")
# example0_without1k<- readRDS("~/ses-1/user_wx/example_RLE_pca_nomed_noinflame.rds")
# tmm
example0_with1k <-readRDS("~/ses-1/user_wx/example_tmm_m7_withinflame.rds")
example0_without1k<- readRDS("~/ses-1/user_wx/example_tmm_m7_noinflame0209.rds")

example_tmm_m7_denovosesDE <- readRDS("~/ses-1/user_wx/example_tmm_m7_denovosesDE.rds") %>% 
  filter(gene_set_name == "ses4_de_novo")
loading_extract = function(data) {
  exB_with <- data %>%
    hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
    unnest_longer(p) %>% 
    # dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
    mutate(p = p.adjust(p, method = "fdr")) %>% 
    # dplyr::select(treatment, gene_set_name, p, p_id, pval2) %>% 
    dplyr::filter(p < threshold) %>% 
    group_by(treatment, gene_set_name) %>% 
    # slice(which.min(p)) %>% 
    mutate(pcmin = p_id %>% str_remove("d") %>% as.numeric()) %>%
    group_by(treatment, gene_set_name) %>%
    mutate(p_no = n()) %>% 
    slice(which.min(pcmin)) %>% 
    ungroup %>% 
    hoist(out, wellloaded = list("result", "m7_ob", 1, "other", "well_loaded")) %>%
    hoist(out, loading_matrix = list("result", "m7_ob", 1, "other", "loadings")) %>%
    mutate(sigpc_wellloadedgenes = wellloaded[[1]][pcmin]) %>% 
    mutate(sigpc_loading = loading_matrix %>% map(~.x[,pcmin])) %>% 
    select(-out, -wellloaded, -loading_matrix)
}



loading_extract(example0_with1k) %>% saveRDS("./user_wx/disease_signature_pca_with1k.rds")
loading_extract(example0_without1k) %>% saveRDS("./user_wx/disease_signature_pca_without1k.rds")
loading_extract(example_tmm_m7_denovosesDE)%>% saveRDS("./user_wx/denovo_pca_withoutallsignature.rds")


loading_extract_all = function(data) {
  exB_with <- data %>%
    hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
    unnest_longer(p) %>% 
    # dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
    mutate(p = p.adjust(p, method = "fdr")) %>% 
    # dplyr::select(treatment, gene_set_name, p, p_id, pval2) %>% 
    dplyr::filter(p < threshold) %>% 
    # group_by(treatment, gene_set_name) %>% 
    # slice(which.min(p)) %>% 
    mutate(pcmin = p_id %>% str_remove("d") %>% as.numeric()) %>%
    # group_by(treatment, gene_set_name) %>%
    # mutate(p_no = n()) %>% 
    # slice(which.min(pcmin)) %>% 
    ungroup %>% 
    hoist(out, wellloaded = list("result", "m7_ob", 1, "other", "well_loaded")) %>%
    hoist(out, loading_matrix = list("result", "m7_ob", 1, "other", "loadings")) %>%
    mutate(sigpc_wellloadedgenes = wellloaded[[1]][pcmin]) %>% 
    # mutate(sigpc_loading = loading_matrix %>% map(~.x[,pcmin])) %>% 
    select(-out, -wellloaded, -pcmin)
}



loading_extract_all(example0_with1k) %>% saveRDS("./user_wx/disease_signature_pca_with1k_allsig.rds")
loading_extract_all(example0_without1k) %>% saveRDS("./user_wx/disease_signature_pca_without1k_allsig.rds")
loading_extract_all(example_tmm_m7_denovosesDE)%>% saveRDS("./user_wx/denovo_pca_withoutallsignature_allsig.rds")



