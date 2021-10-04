
#' ---
#' title: Results of direction dependency
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
source("/home/xu/ses-1/user_wx/printdda.R")
#' ### perform DDAmodel selection
#' * evaluate independence properties of predictors and residuals
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
example0_with1k <- readRDS("~/ses-1/user_wx/example_tmm_m7_withinflame.rds")
example0_without1k <- readRDS("~/ses-1/user_wx/example_tmm_m7_noinflame0209.rds")
threshold = 0.05
nfactors = 10

dda = 
  example0_with1k %>% 
  hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
  unnest_longer(p) %>% 
  mutate(p = p.adjust(p, method = "fdr"),
         gene_set_name = str_c(gene_set_name,"_",p_id),
         "1KI Genes" = "With 1KI Genes") %>% 
  dplyr::filter(p < threshold, p_id %in% c(str_c("d",1:nfactors))) %>% 
  hoist(out, dda_all = list("result", "m7_ob", 1, "dda", "w5bmi_lm", "result")) %>% 
  mutate(dda = map2(.$p_id, .$dda_all, ~ .y[[.x]])) %>% 
  dplyr::select(treatment, gene_set_name, dda) 


for (i in 1 : dim(dda)[1]) {
  
  cat(" ############################################################","\n",
      "treatment is", dda[i, ]$treatment, "gene_set_name is", dda[i, ]$gene_set_name,"\n",
      "############################################################","\n")
  printresid(dda[i, ]$dda[[1]]$resdist)
  print(dda[i, ]$dda[[1]]$indep3)

}
