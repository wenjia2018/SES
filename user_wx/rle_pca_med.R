
#' ---
#' title: RLE normalization PCA mediational analysis
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE


#' #### We do PCA first and then do mediational analysis only for the case where there is a significant total effect
#' *  for significant PCs using 5 mediators
#' *  p: p value for indirected effect / mediated effect;
#' *  result_id: which pc is significant;
#' *  adjP is adjusted for all the disease signatures(11) and 5 mediators using fdr
#' *  results are presented for adjP<0.05 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(here)
library(tidyverse)
# example0_with1k = readRDS("/home/share/scratch/xu/example0_new2signature_withinflam.rds")

example0 <- readRDS("~/ses-1/user_wx/example_RLE_pca_withinflame.rds")

threshold = 0.05

# we use fig1panleB to choose the mediation we wanted to do, that is only for significant total effect
fig1panelB <- example0 %>%
  hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
  unnest_longer(p) %>% 
  dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  mutate(p = p.adjust(p, method = "fdr")) %>% 
  dplyr::filter(p < threshold) %>% 
  dplyr::select(treatment, gene_set_name, p_id) %>% 
  mutate(id = str_c(treatment,"_",gene_set_name,"_",p_id))

extract_pp = function(m, mediator){
  m %>%
    hoist(out, "result") %>%
    hoist(result, "m7_ob") %>% 
    unnest(matches("^m7")) %>% 
    hoist(m7_ob, result = list("mediation", mediator, "result")) %>% 
    filter(result!="NULL") %>% 
    unnest_longer(result) %>% 
    hoist(result, p = "p")%>% 
    hoist(result, detail = "other") %>% 
    unnest_wider("detail") %>% 
    dplyr::select(-out, -result, -m7_ob) %>% 
    filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
    mutate(id = str_c(treatment,"_",gene_set_name,"_",result_id)) %>% 
    right_join(fig1panelB %>% dplyr::select(id)) %>%
    mutate(p = p %>% str_remove("<") %>% as.numeric() %>% format(digits = 3, scientific =T),
           mediator = mediator)
# %>% 
           # adjP = ((p %>% str_remove("<") %>% as.numeric())) %>% p.adjust("fdr")) %>% 
    # filter(adjP<0.05) %>% 
    # mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric
}

mediators = 
  c(
    "stress_perceived_lm",
    "bills_binary",
    "currentsmoke_binary",
    "w5bmi_lm",
    "insurance_lack_binary"
  )
out = mediators %>% set_names() %>%  map(~ extract_pp(example0, .x)) %>% bind_rows() %>% 
  mutate(adjP = ((p %>% str_remove("<") %>% as.numeric())) %>% p.adjust("fdr")) %>% 
  filter(adjP<0.05) %>%
  mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric
out %>% select(-id, -controls) %>% kableExtra::kable() %>% kableExtra::kable_styling()
