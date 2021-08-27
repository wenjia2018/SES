
#' ---
#' title: Evalue for PCA
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE, include=FALSE, comment=NA

#' ### Evalue
#' *  The E-value is defined as the minimum strength of association, on the risk ratio scale,
#'  that an unmeasured confounder would need to have with both the treatment and the outcome to fully explain away a specific
#'  treatment–outcome association, conditional on the measured covariates.
#' * Input required to calculate Evalue: OLS estimate and standard deviation, and standard deviation of the model (outcome)
#'
#' ### approximation or speical treatment made to calculate evalue for causal mediated effect
#' * estimate from mediation package is based on bootstrap. use bootstrap estimate and sd as OLS estimate and sd.  

#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
library(tidyverse)
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(Biobase)
library(ggformula)
library(ggpubr)
library(here)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(tidyverse)
library(stringi)
walk(dir(path = here("R"),full.names = TRUE), source)
example_bmi_m7_evalues_1208 <- readRDS("~/ses-1/user_wx/example_bmi_m7_evalues_1208.rds")

  
  

temp = 
  example_bmi_m7_evalues_1208 %>% 
  hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
  unnest_longer(p) %>% 
  filter(p<0.05/50) %>% 
  hoist(out, evals = list("result", "m7_ob", 1, "evalue")) %>% 
  mutate(evalue = map2(.$p_id, .$evals, ~ .y[[.x]])) %>% 
  dplyr::select(treatment, gene_set_name, p, p_id, evalue)

temp2 = 
  example_bmi_m7_evalues_1208 %>% 
  hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
  unnest_longer(p) %>% 
  filter(p<0.05/50) %>% 
  hoist(out, evals = list("result", "m7_ob", 1, "mediation", "w5bmi_lm", "result")) %>% 
  mutate(evalue = map2(.$p_id, .$evals, ~ .y[[.x]])) %>% 
  dplyr::select(treatment, gene_set_name, p, p_id, evalue)


#' ### Evalue for total effect
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
for (i in 1 : dim(temp)[1]) {
  
  cat(" ############################################################","\n",
      "treatment is", temp[i, ]$treatment, "gene_set_name is", temp[i, ]$gene_set_name,"\n",
      "############################################################","\n")
  cat("total effect","\n","\n")
  
   print(temp[i,]$evalue)
   
   cat("casual mediated effect","\n","\n")
   
  print(temp2[i,]$evalue[[1]]$evalue)
  cat("\n","\n")
  
}
