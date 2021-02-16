
#' ### treatments: black and hispanic

#' ### controls: controls in ses paper + ses + bespoke ancestryPC
#' ### outcome:
#'c(
#'  "ctra_mRNA",
#'  "inflame_mRNA",
#'  "interferon_mRNA",
#'  "AntBIntF_mRNA", 
#'  "inflam1k_mRNA",
#' "aging_mRNA",
#' "aging_up_mRNA",
#'  "aging_down_mRNA"
#')
#'
#' ### analysis:
#' * omnibus P value: association between each disease set and race dummy, 
#' the smallest, whole-genome FDR corrected p-value via whole genome regressions
#'
#' * PCA p value: p value associated between each PC of the disease set and race dummy
#' 
#' * p_eqtl is a sequence of p (0.05, 0.01, 1e-3,...,1e-10) used to choose asscoiated snps
#' 
#' * omnibus p values are genowide corrected p values
#' 
#' * pca regression p value and mediation p values are unadjusted 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(here)
library(tidyverse)
bespoke <- readRDS("~/ses-1/user_wx/bespoke_v2.rds") 

f0 = function(p){
  example0 <- bespoke %>%
    hoist(out, out1 = list("result","example0")) %>%
    select(1:3) %>%
    filter(p_eqtl == p) %>% 
    rename(out=out1) %>%
    unnest(out) %>% 
    filter(control_set =="ancestryPC")
  return(example0)
}

f1 = function(p){
  example1 <- bespoke %>%
    hoist(out, out1 = list("result","example1")) %>%
    select(1:3) %>%
    filter(p_eqtl == p) %>% 
    rename(out=out1) %>% 
    unnest(out) %>% 
    filter(control_set =="ancestryPC") 
  return(example1)
}

mroast_present = function(example1){
  example1 %>%
    hoist(out, gsea = list("result", "gene_set_test")) %>% 
    dplyr::select(1,3,4) %>%
    unnest(gsea) %>% 
    as_tibble() %>% 
    filter(controls == "ancestryPC") %>% 
    dplyr::filter(FDR.Mixed<0.05) %>% 
    kableExtra::kable() %>%
    kableExtra::kable_styling() 
}

m8_present = function(example0){
  
  example0 %>% 
    hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
    dplyr::select(treatment, gene_set_name, p) %>% 
    dplyr::filter(p<0.05) %>%
    rename(p_omnibus = p) %>% 
    # filter(controls == "ancestryPC") %>% 
    kableExtra::kable() %>%
    kableExtra::kable_styling()
}

m7_present = function(example0){
  threshold = 0.05/10
  var = example0 %>%
    hoist(out, var_explained = list("result", "m7_ob", 1, "other", "varexplained")) %>% 
    dplyr::select(treatment, gene_set_name, var_explained)
  
  var$var_explained = var$var_explained %>% map(~ set_names(.x, str_c("d", 1:10)))
  
  var = var %>% unnest_longer(var_explained)
  
  
  gene_list = example0 %>%
    hoist(out, well_loaded = list("result", "m7_ob", 1, "other", "well_loaded")) %>% 
    dplyr::select(treatment, gene_set_name, well_loaded)
  
  gene_list$well_loaded = gene_list$well_loaded %>% map(~ set_names(.x, str_c("d", 1:10)))
  
  gene_list = gene_list %>% unnest_longer(well_loaded)
  
  example0 %>%
    hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
    unnest_longer(p) %>%
    dplyr::select(treatment, gene_set_name, p, p_id) %>% 
    dplyr::filter(p < threshold) %>%
    left_join(var, by = c("treatment", "gene_set_name", "p_id"= "var_explained_id")) %>% 
    left_join(gene_list, by = c("treatment", "gene_set_name", "p_id"= "well_loaded_id")) %>% 
    # filter(controls == "ancestryPC") %>% 
    dplyr::select(1:5) %>% 
    kableExtra::kable() %>%
    kableExtra::kable_styling() 
  
}


med_present = function (focal, example0){
  threshold_med = 0.05
  temp =  example0 %>%
    hoist(out, "result") %>%
    hoist(result, "m7_ob") %>% 
    unnest(matches("^m7")) %>% 
    hoist(m7_ob, result = list("mediation", focal, "result")) %>% 
    filter(result!="NULL")
  
  if(dim(temp)[1]==0){
    return(NULL)
  }else if(temp %>%  unnest_longer(result) %>% dim() %>% .[1] == 0){
    return(NULL)
  }else{
    out = temp %>% 
      unnest_longer(result) %>% 
      hoist(result,p = "p") %>% 
      dplyr::select(3, 4, 6, 8) %>%
      filter(p < threshold_med) %>%
      mutate(mediator = focal)
    return(out)
  }
}



mediators = 
  c(
    "stress_perceived_lm",
    "bills_binary",
    "currentsmoke_binary",
    "w5bmi_lm",
    "insurance_lack_binary",
    "lowbirthweight_binary",
    "high_lowbirth_binary",
    
    "totdiscrim2_binary",  # binary
    "discrim2_binary",   # binary
    "totdiscrim1_category",   # categorical of 4
    # special treatment in mediation
    "totdiscrim1_gamma", #exponential fit glm for mediation
    "totdiscrim2_gamma", #exponential fit glm for mediation
    "countdiscrimwhy_gamma",#exponential fit glm for mediation
    
    "totdiscrim1_pois", #poisson fit glm for mediation
    "totdiscrim2_pois", #poisson fit glm for mediation
    "countdiscrimwhy_pois",#poisson fit glm for mediation
    
    "totdiscrim1_category", #ordered logistic fit glm for mediation
    "totdiscrim2_category", #ordered logistic fit glm for mediation
    "countdiscrimwhy_category"#ordered logistic fit glm for mediation
  )
#' ## eqtl_p=0.05
#+ echo=F, eval=T, warning=FALSE, message=FALSE
# mroast_present(f0(0.05))
#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m8_present(f0(0.05))

#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m7_present(f0(0.05))

#' ### mediation 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
mediators %>%
  set_names() %>%
  map(med_present, f0(0.05)) %>% 
  keep(function(x) is_tibble(x)==TRUE) %>% 
  keep(function(x) dim(x)[1]>0) %>% 
  bind_rows() %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()


#' ## eqtl_p=0.01
#+ echo=F, eval=T, warning=FALSE, message=FALSE
# mroast_present(f0(0.05))
#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m8_present(f0(0.01))

#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m7_present(f0(0.01))

#' ### mediation 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
mediators %>%
  set_names() %>%
  map(med_present, f0(0.01)) %>% 
  keep(function(x) is_tibble(x)==TRUE) %>% 
  keep(function(x) dim(x)[1]>0) %>% 
  bind_rows() %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ## eqtl_p=1e-3
#+ echo=F, eval=T, warning=FALSE, message=FALSE
# mroast_present(f0(0.05))
#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m8_present(f0(1e-3))

#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m7_present(f0(1e-3))

#' ### mediation 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
mediators %>%
  set_names() %>%
  map(med_present, f0(1e-3)) %>% 
  keep(function(x) is_tibble(x)==TRUE) %>% 
  keep(function(x) dim(x)[1]>0) %>% 
  bind_rows() %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()
#' ## eqtl_p=1e-4
#+ echo=F, eval=T, warning=FALSE, message=FALSE
# mroast_present(f0(0.05))
#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m8_present(f0(1e-4))

#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m7_present(f0(1e-4))

#' ### mediation 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
mediators %>%
  set_names() %>%
  map(med_present, f0(1e-4)) %>% 
  keep(function(x) is_tibble(x)==TRUE) %>% 
  keep(function(x) dim(x)[1]>0) %>% 
  bind_rows() %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ## eqtl_p=1e-5
#+ echo=F, eval=T, warning=FALSE, message=FALSE
# mroast_present(f0(0.05))
#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m8_present(f0(1e-5))

#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m7_present(f0(1e-5))

#' ### mediation 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
mediators %>%
  set_names() %>%
  map(med_present, f0(1e-5)) %>% 
  keep(function(x) is_tibble(x)==TRUE) %>% 
  keep(function(x) dim(x)[1]>0) %>% 
  bind_rows() %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ## eqtl_p=1e-6
#+ echo=F, eval=T, warning=FALSE, message=FALSE
# mroast_present(f0(0.05))
#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m8_present(f0(1e-6))

#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m7_present(f0(1e-6))

#' ### mediation 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
mediators %>%
  set_names() %>%
  map(med_present, f0(1e-6)) %>% 
  keep(function(x) is_tibble(x)==TRUE) %>% 
  keep(function(x) dim(x)[1]>0) %>% 
  bind_rows() %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### eqtl_p=1e-7
#+ echo=F, eval=T, warning=FALSE, message=FALSE
# mroast_present(f0(0.05))
#' ## omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m8_present(f0(1e-7))

#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m7_present(f0(1e-7))

#' ### mediation 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
mediators %>%
  set_names() %>%
  map(med_present, f0(1e-7)) %>% 
  keep(function(x) is_tibble(x)==TRUE) %>% 
  keep(function(x) dim(x)[1]>0) %>% 
  bind_rows() %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ## eqtl_p=1e-8
#+ echo=F, eval=T, warning=FALSE, message=FALSE
# mroast_present(f0(0.05))
#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m8_present(f0(1e-8))

#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m7_present(f0(1e-8))

#' ### mediation 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
mediators %>%
  set_names() %>%
  map(med_present, f0(1e-8)) %>% 
  keep(function(x) is_tibble(x)==TRUE) %>% 
  keep(function(x) dim(x)[1]>0) %>% 
  bind_rows() %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ## eqtl_p=1e-9
#+ echo=F, eval=T, warning=FALSE, message=FALSE
# mroast_present(f0(0.05))
#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m8_present(f0(1e-9))

#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m7_present(f0(1e-9))

#' ### mediation 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
mediators %>%
  set_names() %>%
  map(med_present, f0(1e-9)) %>% 
  keep(function(x) is_tibble(x)==TRUE) %>% 
  keep(function(x) dim(x)[1]>0) %>% 
  bind_rows() %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ## eqtl_p=1e-10
#+ echo=F, eval=T, warning=FALSE, message=FALSE
# mroast_present(f0(0.05))
#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m8_present(f0(1e-10))

#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m7_present(f0(1e-10))

#' ### mediation 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
mediators %>%
  set_names() %>%
  map(med_present, f0(1e-10)) %>% 
  keep(function(x) is_tibble(x)==TRUE) %>% 
  keep(function(x) dim(x)[1]>0) %>% 
  bind_rows() %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()
