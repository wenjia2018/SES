
#' ---
#' title: oblique rotation PCA mediational analysis
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE


#' We do PCA for two new disease signature, which are constructed by DE ses, using oblique rotation and then do mediational analysis
#'  for significant PCs using 5 mediators
#'  
#' p: p value for indirected effect / mediated effect;
#' result_id: which pc is significant;
#' adjP is adjusted for 13 disease signatures and 5 mediators, ie. adjP = p /65;
#' results are presented for adjP<0.05 
library(here)
library(tidyverse)
# example0_with1k = readRDS("/home/share/scratch/xu/example0_new2signature_withinflam.rds")

example0 <- readRDS("/home/share/scratch/example0_without_1KI_0110.rds")
example0 <- readRDS("/home/share/scratch/example0_with_1KI_0110.rds")
example0 <- readRDS("/home/share/scratch/example0_without_1KI_de_novo_0110.rds")
example0 <- readRDS("/home/share/scratch/example0_without_1KI_de_novo2_0110.rds")
#' ## keep inflamation signatures
#' 
threshold = 0.05/10/5
fig1panelB <- example0 %>%
  hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
  unnest_longer(p) %>% 
  dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
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
    hoist(result, Prop_Mediated = "other") %>% 
    dplyr::select(1:5, 7) %>% 
    filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
    mutate(id = str_c(treatment,"_",gene_set_name,"_",result_id)) %>% 
    right_join(fig1panelB %>% dplyr::select(id)) %>% 
    mutate(p = p %>% str_remove("<") %>% as.numeric() %>% format(digits = 3, scientific =T),
           adjP = ((p %>% str_remove("<") %>% as.numeric()) * dim(.)[1]),
           adjP = ifelse(adjP>1, 1, adjP)) %>% 
    filter(adjP<0.05) %>% 
    mutate(adjP = adjP %>% format(digits = 3, scientific =T)) # convert from strings to numeric
}
#' ### bmi at wave 5

#+ echo=F, eval=T, warning=FALSE, message=FALSE
# w5bmi = extract_pp(example0, "w5bmi")
# w5bmi %>% kableExtra::kable() %>% kableExtra::kable_styling()
mediators = 
  c(
    "stress_perceived",
    "bills",
    "currentsmoke",
    "w5bmi",
    "insurance_lack"
  )
out = mediators %>% set_names() %>%  map(~ extract_pp(example0, .x))
out$w5bmi %>% kableExtra::kable() %>% kableExtra::kable_styling()

