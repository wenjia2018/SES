#' ### treatments: black and hispanic

#' ### controls: controls in ses paper + ses + bespoke ancestryPC
#' * ancestryPC is obtained by snps which are in 2000 upstream and downstream region of the gene (in outcome set) transcription start site
#'
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

#' 
#' * omnibus P value: association between each disease set and race dummy, 
#' the smallest, whole-genome FDR corrected p-value via whole genome regressions
#'
#' * PCA p value: p value association between each PC of the disease set and race dummy
#' 
#' * omnibus p values are genowide corrected p values
#' 
#' * pca regression p value and mediation p values are unadjusted 

#+ echo=F, eval=T, warning=FALSE, message=FALSE

bespoke <- readRDS("~/ses-1/user_wx/bespoke_g2s.rds")
control = "ancestryPC"
threshold = 0.05/10
threshold_med = 0.05
  example0 <- bespoke %>%
    hoist(out, out1 = list("result","example0")) %>%
    select(1:2) %>%
    rename(out=out1) %>%
    unnest(out) 


  example1 <- bespoke %>%
    hoist(out, out1 = list("result","example1")) %>%
    select(1:2) %>%
    rename(out=out1) %>% 
    unnest(out) 

mroast_present = function(example1){
  example1 %>%
    hoist(out, gsea = list("result", "gene_set_test")) %>% 
    dplyr::select(1,3,4) %>%
    unnest(gsea) %>% 
    as_tibble() %>% 
    filter(control_set == control) %>% 
    dplyr::filter(FDR.Mixed<0.05) %>% 
    kableExtra::kable() %>%
    kableExtra::kable_styling() 
}

m8_present = function(example0, control){
  
  out = example0 %>% 
    hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
    dplyr::select(treatment, gene_set_name, p, control_set) %>% 
    dplyr::filter(p<0.05) %>%
    rename(p_omnibus = p) %>% 
    filter(control_set == control) 
  return(out)
}

m7_present = function(example0, control){

  var = example0 %>%
    hoist(out, var_explained = list("result", "m7_ob", 1, "other", "varexplained")) %>% 
    dplyr::select(treatment, gene_set_name, var_explained, control_set)
  
  var$var_explained = var$var_explained %>% map(~ set_names(.x, str_c("d", 1:10)))
  
  var = var %>% unnest_longer(var_explained)
  
  
  gene_list = example0 %>%
    hoist(out, well_loaded = list("result", "m7_ob", 1, "other", "well_loaded")) %>% 
    dplyr::select(treatment, gene_set_name, well_loaded, control_set)
  
  gene_list$well_loaded = gene_list$well_loaded %>% map(~ set_names(.x, str_c("d", 1:10)))
  
  gene_list = gene_list %>% unnest_longer(well_loaded)
  
  out = example0 %>%
    hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
    unnest_longer(p) %>%
    dplyr::select(treatment, gene_set_name, p, p_id, control_set) %>% 
    dplyr::filter(p < threshold) %>%
    left_join(var, by = c("treatment", "gene_set_name", "p_id"= "var_explained_id", "control_set")) %>% 
    left_join(gene_list, by = c("treatment", "gene_set_name", "p_id"= "well_loaded_id", "control_set")) %>% 
    filter(control_set == control) %>%
    dplyr::select(1:6) %>% 
    rename(p_pca = p) 
  return(out)
  
}


med_extact = function (focal, example0, control){
  out = example0 %>%
    hoist(out, med=list("result", "m7_ob", 1, "mediation", focal, "result")) %>% 
    mutate(mediator = focal) %>% 
    filter(control_set== control) %>% 
    select(-out)
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
    
    "totdiscrim2_category", #ordered logistic fit glm for mediation
    "countdiscrimwhy_category"#ordered logistic fit glm for mediation
  )

#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE

m8_present(example0, control) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE

m7_present(example0, control) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### mediation 
#+ echo=F, eval=T, warning=FALSE, message=FALSE

 mediators %>% 
  map_df(med_extact, example0, control) %>% 
  filter(med!="NULL") %>%
  unnest_longer(med) %>% 
  hoist(med, p = "p") %>%
  filter(p < threshold_med) %>% 
  rename(p_med= p) %>% 
   select(2,3,5,7,8,9)%>% 
   kableExtra::kable() %>%
   kableExtra::kable_styling()
