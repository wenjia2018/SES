#' ---
#' title: var explained
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE, include=FALSE, comment=NA

#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA

library(tidyverse)
threshold = 0.05
# rle
# example0_with1k <- readRDS("~/ses-1/user_wx/example_RLE_pca_nomed_withinflame.rds")
# example0_without1k<- readRDS("~/ses-1/user_wx/example_RLE_pca_nomed_noinflame.rds")
# tmm
example0_with1k <-readRDS("~/ses-1/user_wx/example_tmm_m7_withinflame.rds")
example0_without1k<- readRDS("~/ses-1/user_wx/example_tmm_m7_noinflame0209.rds")

varexplained_extract = function(data) {
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
    hoist(out, varexplained = list("result", "m7_ob", 1, "other", "varexplained")) %>%
    mutate(sigpc_varexplained = varexplained[[1]][pcmin]) %>% 
    select(-out, -varexplained)
  return(exB_with)
}

#' ### with 1KI 
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
with = 
  varexplained_extract(example0_with1k) %>% 
  select(treatment, gene_set_name, sigpc_varexplained) 
varwith = with$sigpc_varexplained %>% map_dbl(~.x)

cat("mean of varexplained is", mean(varwith))
cat("range of varexplained is", "[", min(varwith), ",", max(varwith), "]")
with %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()
#' ### without 1KI 
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA 
without = 
  varexplained_extract(example0_without1k) %>% 
  select(treatment, gene_set_name, sigpc_varexplained) 
varwithout = without$sigpc_varexplained %>% map_dbl(~.x)

cat("mean of varexplained is", mean(varwithout))
cat("range of varexplained is", "[", min(varwithout), ",", max(varwithout), "]")

without%>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling() 
