# not finished!!!
#' ---
#' title: coefficient plot 
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE


#' ### stratification--raceethnicity: 
#' * "NonHwhite"
#' * "NonHblack",
#' * "Hispanic"
#' 
#' 
#' ### treatments: 
#' * "color_byinterviewer3_White"
#' * "color_byinterviewer3_LightMed",
#' * "color_byinterviewer3_DarkBlack"



control = "ancestryPC_ses"
threshold = 0.05/10
threshold_med = 0.05
p_eqtl = c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
#' ### controls: 
#' 
#' * ancestryPC_ses : basic + ses + bespoke ancestryPC
#' 
#' ### outcome: aging and clusters
#'
#' ### analysis:
#' * omnibus P value: association between each disease set and color related snps, 
#' the smallest, whole-genome FDR corrected p-value via whole genome regressions
#'
#' * PCA p value: p value associated between each PC of the disease set and color related snps
#' 
#' * p_eqtl is a sequence of p (0.05, 0.01, 1e-3,...,1e-10) used to choose asscoiated snps
#' 
#' * omnibus p values are genowide corrected p values
#' 
#' * pca regression p value and mediation p values are unadjusted 


#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(here)
library(tidyverse)
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase)
library(enrichplot)
library(dbr) # my package
walk(dir(path = here("R"), full.names = TRUE), source)
load_data(reconciled = FALSE, remove_inflam = FALSE)
source("/home/xu/ses-1/user_wx/extract_v2.R")
control = "ancestryPC_ses"
threshold = 0.05/10
threshold_med = 0.05
p_eqtl = c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
#' ## Strata: Nonhispanic White(white as reference)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

#' ### omnibus regression
#+ echo=F, eval=T, warning=FALSE, message=FALSE

bespoke <- readRDS("~/ses-1/user_wx/colorbinarycontinuous_bespoke_NonHwhite_strata_25.02.2021.rds")

outomni = p_eqtl %>% map(outm8_allsig, control, bespoke)

# outomni = p_eqtl %>% map(outm8_allsig, control)

temp = outomni %>%
  bind_rows() 
# mutate(gene= gene_sig %>% intersect(signatures$outcome_set[[.$table1]])) %>%
temp %>%
  mutate(gene= map2(table1, gene_sig, ~ intersect(signatures$outcome_set %>% pluck(.x), .y))) %>% 
  filter(gene %>% map_dfc(~ length(.))>0) %>% 
  unnest_longer(gene) %>% 
  mutate(adj.P = map2_dbl(gene, ttT, ~ filter(.y, gene==.x) %>% pull("adj.P.Val")),
         logFC = map2_dbl(gene, ttT, ~ filter(.y, gene==.x) %>% pull("logFC")),
         CI.L = map2_dbl(gene, ttT, ~ filter(.y, gene==.x) %>% pull("CI.L")),
         CI.R= map2_dbl(gene, ttT, ~ filter(.y, gene==.x) %>% pull("CI.R"))) %>% 
  select(-controls, -ttT, -gene_sig, -out) %>% 
  arrange(table1) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()

gene1 = temp %>%
  mutate(gene= map2(table1, gene_sig, ~ intersect(signatures$outcome_set %>% pluck(.x), .y))) %>% 
  filter(gene %>% map_dfc(~ length(.))>0)  %>%
  unnest_longer(gene) %>% 
  pull(gene) %>% 
  unique


dt = temp %>% 
  mutate(gene= map(table1, ~ intersect(signatures$outcome_set %>% pluck(.x), gene1))) %>%  
  unnest_longer(gene) %>% 
  filter(!is.na(gene)) %>% 
  mutate(adj.P = map2_dbl(gene, ttT, ~ filter(.y, gene==.x) %>% pull("adj.P.Val")),
         logFC = map2_dbl(gene, ttT, ~ filter(.y, gene==.x) %>% pull("logFC")),
         CI.L = map2_dbl(gene, ttT, ~ filter(.y, gene==.x) %>% pull("CI.L")),
         CI.R= map2_dbl(gene, ttT, ~ filter(.y, gene==.x) %>% pull("CI.R"))) %>% 
  select(-controls, -ttT, -gene_sig) %>% 
  arrange(table1) 
dt %>% 
  mutate(treatment = treatment %>% 
           str_remove("color_byinterviewer5_") %>% 
           factor(levels = c("Black", "Dark", "Medium", "Light"))) %>% 
  ggplot(aes(x = treatment, y = logFC, ymin = CI.L, ymax = CI.R)) +
  geom_pointrange()+
  facet_wrap(~ p_eqtl + table1, ncol = 4) +  
  labs(title = str_c("non hispanic white: ", gene1)) 


#' ## Strata: Nonhispanic black(Medium as reference)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
bespoke <- readRDS("~/ses-1/user_wx/color3_bespoke_NonHblack_strata_25.02.2021.rds")

outomni = p_eqtl %>% map(outttT,control, bespoke)
# outomni = p_eqtl %>% map(outm8_allsig, control)
temp = outomni %>%
  bind_rows() 
# mutate(gene= gene_sig %>% intersect(signatures$outcome_set[[.$table1]])) %>%
temp %>%
  mutate(gene= map2(table1, gene_sig, ~ intersect(signatures$outcome_set %>% pluck(.x), .y))) %>% 
  filter(gene %>% map_dfc(~ length(.))>0) %>% 
  unnest_longer(gene) %>% 
  mutate(adj.P = map2_dbl(gene, ttT, ~ filter(.y, gene==.x) %>% pull("adj.P.Val")),
         logFC = map2_dbl(gene, ttT, ~ filter(.y, gene==.x) %>% pull("logFC")),
         CI.L = map2_dbl(gene, ttT, ~ filter(.y, gene==.x) %>% pull("CI.L")),
         CI.R= map2_dbl(gene, ttT, ~ filter(.y, gene==.x) %>% pull("CI.R"))) %>% 
  select(-controls, -ttT, -gene_sig) %>% 
  arrange(gene) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()

gene1 = temp %>%
  mutate(gene= map2(table1, gene_sig, ~ intersect(signatures$outcome_set %>% pluck(.x), .y))) %>% 
  filter(gene %>% map_dfc(~ length(.))>0)  %>%
  unnest_longer(gene) %>% 
  pull(gene) %>% 
  unique
# for (i in 1: length(gene1)){
gene_i = gene1[1]
dt = temp %>% 
  mutate(gene= map(table1, ~ intersect(signatures$outcome_set %>% pluck(.x), gene_i))) %>%  
  unnest_longer(gene) %>% 
  filter(!is.na(gene)) %>% 
  mutate(adj.P = map2_dbl(gene, ttT, ~ filter(.y, gene==.x) %>% pull("adj.P.Val")),
         logFC = map2_dbl(gene, ttT, ~ filter(.y, gene==.x) %>% pull("logFC")),
         CI.L = map2_dbl(gene, ttT, ~ filter(.y, gene==.x) %>% pull("CI.L")),
         CI.R= map2_dbl(gene, ttT, ~ filter(.y, gene==.x) %>% pull("CI.R"))) %>% 
  select(-controls, -ttT, -gene_sig) %>% 
  arrange(gene) 
# f = 
dt %>% 
  filter(!is.na(adj.P)) %>% 
  mutate(treatment = treatment %>% 
           str_remove("color_byinterviewer3_") %>% 
           factor(levels = c("Light Medium"))) %>% 
  ggplot(aes(x = treatment, y = logFC, ymin = CI.L, ymax = CI.R)) +
  geom_pointrange() +
  facet_wrap(~ p_eqtl + table1, ncol = 4) +  
  labs(title = str_c("non hispanic black: ", gene_i)) 
# print(f)
# }



