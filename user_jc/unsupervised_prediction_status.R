#' ---
#' title: Examples
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#' Set global options 
#+ setup, warning=FALSE, message=FALSE
# knitr::opts_chunk$set(echo = FALSE)

set.seed(1234) 
library(here)
library(tidyverse) 
library(tidymodels)
library(rlang)
library(Biobase) 
library(ggformula)
walk(dir(path = here("R"), full.names = TRUE), source) 

load_data(reconciled = FALSE, remove_inflam = FALSE)
define_treatments_and_controls()
sigs = signatures$outcome_set[c(table1, "ctra_mRNA")]

controls = controls %>% pluck("basic")
# treatment = "ses_sss_composite"

get_preds =
  function(treatment, residualized = "residualized") {
    # otherwise "not_residualized" wrt controls
    
    dt =
      prepro(reduce(sigs, union),
             treatment, controls) %>%
      rename(y = !!treatment) %>%
      na.omit()
    
    # possibly redefine  outcome as residualized outcome
    if(residualized == "residualized") dt$y = resid(lm(as.formula(str_c("y", " ~ ", str_c(controls, collapse = " + "))), dt))
    
    split    = initial_split(dt)
    training = training(split)
    testing  = testing(split)
    
    folds    = vfold_cv(training, v = 10)
    control  = control_resamples(save_pred = TRUE)
    
    args = list(training = training, folds = folds, control = control)
    preds = map(sigs[table1], safely(fit_m0), args)
    saveRDS(preds, here("user_jc", str_c("preds", "_", residualized , "_", treatment, ".rds")))
    
  }

############################################################
# COMPUTATIONALLY COSTLY!!! AND SAVES TO DISK
############################################################
list("ses_sss_composite",
     "sss_5",
     "SEI_ff5",
     # "edu_max",
     "income_hh_ff5") %>%
  walk(get_preds, residualized = "residualized") %>%
  walk(get_preds, residualized = "not_residualized")

############################################################
# EXAMINE
############################################################

get_mean_rsq =  
  function(x, y) {
    select(pluck(x, "result"), -fits, -preds) %>% 
      unnest(metrics) %>%
      filter(.metric == "rsq") %>% 
      mutate(gene_set = y)  %>% 
      select(gene_set, preprocessor, object, mean_rsq = mean) %>% 
      arrange(-mean_rsq)
  }

in_sample_svm_summary = 
  function(preds){ 
    accuracies =
      imap(preds, get_mean_rsq) %>% 
      map_df(filter, object == "svm") %>% 
      arrange(-mean_rsq)  
    return(accuracies = accuracies)
  }

svm_summaries = 
  list.files("user_jc",
             pattern = "^preds*",
             full.names = TRUE) %>% 
  set_names() %>% 
  map(compose(in_sample_svm_summary, readRDS))


if(0){
  
  
  experiment = 
    function(x) { 
      # for a big speed up
      select(x, preprocessor, object, preds) %>% 
        unnest(preds) %>% 
        gf_point(y...4 ~ .pred) %>% 
        gf_facet_wrap(~preprocessor + object, ncol = 2) %>% 
        print
    }
  debugonce(experiment)
  preds %>% map(experiment)
  
  
  
  if(0){
    
    # residual analysis: 
    # what is the expected prediction error (e.g. rmse), on average over all input space X?
    # how does the local prediction error (out-of-sample residual) vary over input-space?
    # how does prediction error vary? systematically with some variable?
    preds %>% 
      pluck("kidney_transplant_tolerance_mRNA") %>% 
      select(preprocessor, preds) %>% 
      unnest(preds) %>% 
      gather("gene", "mRNA", signatures$outcome_set[[1]]) %>%
      select(gene, mRNA, y, .pred) %>%
      gf_point(y-.pred ~ mRNA) %>% 
      gf_facet_wrap(~gene)
    
  }
  
}