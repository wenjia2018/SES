
#' ---
#' title: skin color3 and race interaction model
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

#' ### treatments: 
#' * Among non hispanic white and hispanic race and Skin Color interaction model, use hispanic and white as reference
#' * construct the following dummies variables
#' * "color_byinterviewer3_DarkBlack", "color_byinterviewer3_LightMed","raceethnicity_NonHwhite",
#' "raceethnicity__color_byinterviewer3_NonHwhite|DarkBlack","raceethnicity__color_byinterviewer3_NonHwhite|LightMed"

#+ echo=F, eval=T, warning=FALSE, message=FALSE
load_data(reconciled = FALSE, remove_inflam = FALSE)
pData(dat) <- pData(dat) %>%
  # mutate_at(.vars = vars(matches("income|composite|SEI")),
  #           .funs = list(~ .x %>% ntile(4) %>% as.factor())
  # ) %>%
  mutate(raceethnicity = re %>%
           fct_recode(
             NonHwhite = "1",
             # white nonhispanic
             NonHblack = "2",
             # black nonhispanic
             NULL = "3",
             # asian nonhispanic
             NULL = "4",
             # other nonhispanic
             Hispanic = "5"
             # hispanic
           ) %>%
           relevel(ref = "NonHwhite"),
         color_byinterviewer3 = H3IR17 %>%
           as.character() %>% 
           as.factor %>% 
           fct_collapse(
             DarkBlack = c("1", "2"),
             LightMed = c("3", "4"),
             White = "5") %>%
           relevel(ref = "White"))


keep = (pData(dat)$raceethnicity !="NonHblack" & 
          !is.na(pData(dat)$color_byinterviewer3)) %>% ifelse(is.na(.), FALSE, .)
dat <- dat[, keep]

pData(dat) <- pData(dat) %>%
  finalfit::ff_interaction(raceethnicity, color_byinterviewer3, levels_sep = "|", var_sep = "__") %>% 
  mutate_at(.vars = vars(raceethnicity__color_byinterviewer3),
            .funs = list(~ fct_recode(., NULL = "NonHblack|White", NULL = "NonHblack|DarkBlack", NULL = "NonHblack|LightMed"))) %>% 
  fastDummies::dummy_cols(select_columns = c("raceethnicity", "raceethnicity__color_byinterviewer3")) 

descr::crosstab(dat$color_byinterviewer3,  dat$raceethnicity, prop.r = T, prop.c = T, prop.t = FALSE, plot = F)



control = "ancestryPC_ses"
threshold = 0.05/10
threshold_med = 0.05
#' ### controls: 
#' 
#' * ancestryPC_ses : basic + ses + bespoke ancestryPC
#' 
#' ### outcome:
#'c(
#'  "ctra",
#'  "inflamation1k",
#' "aging"
#')
#'
#' ### analysis:
#' * omnibus P value: association between each disease set and each treatment, 
#' the smallest, whole-genome FDR corrected p-value via whole genome regressions
#'
#' * PCA p value: p value associated between each PC of the disease set and each treatment
#' 
#' * p_eqtl is a sequence of p (0.05, 0.01, 1e-3,...,1e-10) used to choose asscoiated snps
#' 
#' * omnibus p values are genowide corrected p values
#' 
#' * pca regression p value and mediation p values are unadjusted 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(here)
library(tidyverse)
source("/home/xu/ses-1/user_wx/extract_v2.R")

color3_and_race3_NonHwhitehispanic_bespoke_31.03.2021 <- readRDS("/home/xu/ses-1/user_wx/color3_and_race3_NonHwhitehispanic_bespoke_31.03.2021.rds")
bespoke = color3_and_race3_NonHwhitehispanic_bespoke_31.03.2021 %>% filter(table1 %>% str_detect("aging|whole_genome"))
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
p_eqtl = c(0.05, 1e-2)

outomni = p_eqtl %>% map(m8_present, control, bespoke)

outomni %>%
  bind_rows() %>%
  dplyr::filter(p_omnibus < 0.05) %>%
  select(-gene_sig, -logFC) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_omnibus) %>%
  arrange(gene_set_name) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_omnibus=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### PCA regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE

# outm7pca = function(p, control){
#   m7_present(f0(p), control) %>% mutate(p_eqtl = p)
# }

outm71 = p_eqtl %>% map(outm7pca, control, bespoke)


outm71 %>%
  bind_rows() %>%
  dplyr::filter(p_pca < threshold) %>%
  select(-coef, -loadings, - well_loaded) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_pca) %>%
  arrange(gene_set_name) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_pca=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()
#' ### mediation 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
# outm7med = function(p, mediators, control){
#   out = mediators %>% 
#     map_df(med_extact, f0(p), control) %>% 
#     filter(med!="NULL") %>%
#     unnest_longer(med) %>% 
#     hoist(med, p = "p") %>%
#     filter(p < threshold_med) %>% 
#     rename(p_med= p)
# }

outm72 = p_eqtl %>% map_df(outm7med, mediators, control, bespoke)

outm72 %>% 
  dplyr::select(treatment, gene_set_name, med_id, mediator,p_med, p_eqtl) %>% 
  dplyr::filter(p_med < threshold_med) %>%
  mutate_at(.vars = vars(c("p_med")), .funs = funs(. %>%  str_remove("<") %>% as.numeric %>% format(scientific =TRUE))) %>% 
  pivot_wider(names_from = p_eqtl, values_from = p_med) %>%
  arrange(gene_set_name) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ format(.x, scientific = TRUE, digits = 2))) %>%
  mutate_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_med=",.x))) %>%
  rename_at(vars(starts_with("0.0")|starts_with("1e")), .funs = list( ~ str_c("p_eqtl=",.x))) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### pca linear hypothesis "raceethnicity__color_byinterviewer3_NonHwhite|DarkBlack" and "raceethnicity__color_byinterviewer3_NonHwhite|LightMed" 

outlh = p_eqtl %>% map_df(outm7_lh, control, bespoke)

outlh2 %>% 
  # filter(treatment == "color_byinterviewer3_LightMed") %>%
  unnest(lh) %>% 
  unnest(lh) %>% 
  filter(!is.na(Df)) %>% 
  filter(`Pr(>F)` <0.005)

#' ### omnibus linear hypothesis "raceethnicity__color_byinterviewer3_NonHwhite|DarkBlack" and "raceethnicity__color_byinterviewer3_NonHwhite|LightMed" 

outomni_lh = p_eqtl %>% map_df(outm10_lh, control, bespoke)

outomni_lh %>% 
  filter(!is.na(lh)) 


