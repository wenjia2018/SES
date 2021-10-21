library(here)
library(data.table)
# setDTthreads(threads = 20)
library(tidyverse)
library(EValue)
library(rlang)
library(skimr)
library(furrr)
library(limma)
# library(recipes)
# recipes and Evalue has conflict, after loading recipes, evalue package doesnot work with the error Error: $ operator is invalid for atomic vectors
library(parsnip)
library(workflows)
library(Biobase)


walk(dir(path = here("R"),full.names = TRUE), source)


############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################
# choose normalization methods for downstream analysis
tmm = TRUE
rle = FALSE
log2cpm = FALSE
# which PCA to perform
oblimin = FALSE
nn = TRUE
# explicitly assign ncomp as the smallest number of table signatures gene numbers
ncomp = 10
# for doing genowide DE analysis only
normalization_bydesign = FALSE
# specify if subjects with disease shall be removed
remove_diseased_subjects = TRUE
load_data(reconciled = FALSE, remove_inflam = FALSE)
define_treatments_and_controls()
recode_variables_in_dat()
print(abbreviations)
mediation_mean = FALSE
mediation_each_gene = FALSE
if(denovo <- FALSE){
  Significant <- readRDS("/home/share/scratch/Clustering/DeNovo/Significant.rds")
  full_clus_res <- readRDS("/home/share/scratch/Clustering/DeNovo/full_clus_res.rds")
  average_expr <- readRDS("/home/share/scratch/Clustering/DeNovo/average_expr.rds")
  gene_set_name = average_expr %>% names()
  args = crossing(treatment, gene_set_name, controls)
}
if(without1k <- TRUE){
  Significant <- readRDS("/home/share/scratch/Clustering/without_1k/Significant.rds")
  full_clus_res <- readRDS("/home/share/scratch/Clustering/without_1k/full_clus_res.rds")
  average_expr <- readRDS("/home/share/scratch/Clustering/without_1k/average_expr.rds")
  
}
if(with1k <- FALSE){
  Significant <- readRDS("/home/share/scratch/Clustering/1k/Significant.rds")
  full_clus_res <- readRDS("/home/share/scratch/Clustering/1k/full_clus_res.rds")
  average_expr <- readRDS("/home/share/scratch/Clustering/1k/average_expr.rds")
  
}
average_expr = average_expr %>% map(~ .x %>% t %>% as.data.frame)

get_sig_clus = function(treatment, Significant, full_clus_res){
  sig_treatment = Significant %>% map(~ pluck(.x, treatment))
  no_clus = full_clus_res %>% map(max)
  full_clus = no_clus %>% map(~ seq(1:.x))
  sig_clus = map2(sig_treatment, full_clus, ~ (.y %in%.x))
  
}

pheno = pData(dat)

subData = pheno %>% dplyr::select(all_of(treatment), batch, all_of(controls$basic))

non_missing = complete.cases(subData)

dat@phenoData@data = pheno[non_missing,]

cluster_reg = function(dt, outcome) {
  covariates = colnames(dt) %>% str_subset(outcome, negate = TRUE) %>% str_c(collapse= " + ")
  formula_lm = outcome %>%  str_c("~", covariates) %>% as.formula
  
  lm(formula_lm, dt)
}
extract_cluster_reg = function(m, out = NULL){
  extract_anova = function(x) anova(x) %>% tidy %>% filter(str_detect(term, "treatment"))
  extract_t = function(x) broom::tidy(x) %>% filter(str_detect(term, "treatment"))
  
  # Univariate parametric
  out$detail$anova = m %>% extract_anova
  out$detail$t = m %>% extract_t
  out$p = out$detail$anova %>% pluck("p.value")
  out$evalue = m %>% cal_evalue
}


fit_cluster_reg = function(treatment, gene_set_name, controls, out = NULL){
  clus = average_expr[[gene_set_name]] 
  colnames(clus) = str_c("d", 1:dim(clus)[2])
  clus = rownames_to_column(clus, var = "AID")
  datt_reg  = pData(dat) %>% 
    dplyr::select(AID = AID, all_of(controls), all_of(treatment)) %>% 
    rename(treatment = treatment) 
  
  if(remove_diseased_subjects) {
    datt_reg =
      datt_reg %>% 
      remove_diseased_subjects_from_datt(gene_set_name, controls)
  }
  
  outcome = colnames(clus %>% select(-AID)) %>% set_names() 
  
  convert_outcome_clus = function(outcome, datt_reg){
    datt_pca = left_join(datt_reg, clus %>% select(AID, !!outcome)) %>% select(-AID)
  }
  
  datt_pca = outcome %>% map(convert_outcome_clus, datt = datt_reg)
  
  sig = get_sig_clus(treatment, Significant, full_clus_res) %>% pluck(gene_set_name)
  
  pmap(list(datt_pca[sig], outcome[sig]), cluster_reg) %>% map(extract_cluster_reg)
}
# debugonce(fit_cluster_reg)

# plan(multicore, workers = 10)
example0 =
  args %>%
  filter(
    is.element(gene_set_name, table1) &
      # gene_set_name == "SES  Composite",
      names(controls) == "basic") %>%
  # filter(gene_set_name=="Depression_mRNA") %>% 
  mutate(out = pmap(., safely(fit_cluster_reg)),
         controls = names(controls))

example0 %>% saveRDS("./user_wx/Evalue_cluster_without1k_v2.rds")