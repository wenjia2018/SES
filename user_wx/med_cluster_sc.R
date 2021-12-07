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
normalization_bydesign = TRUE
# specify if subjects with disease shall be removed
remove_diseased_subjects = TRUE
load_data(reconciled = FALSE, remove_inflam = FALSE)
define_treatments_and_controls_sc()
# recode_variables_in_dat_racedummy()
gene_set_name = "aging_mRNA"
p_eqtl = 0.05
ancestryPC <- get_PC_dim(gene_set_name, p_eqtl)
define_treatments_and_controls_sc_bespoke(gene_set_name, ancestryPC)
custom_PCA <- readRDS(str_c("/home/share/dna_ancestry/dna/custom_PCA_", gene_set_name, "_", p_eqtl, ".rds"))
custom_PCA <- custom_PCA %>%
  select(-fid) %>%
  mutate(AID = AID %>% as.character())
recode_variables_in_dat_racedummy_bespoke(custom_PCA)
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
if(without1k <- FALSE){
  Significant <- readRDS("/home/share/scratch/Clustering/without_1k/Significant.rds")
  full_clus_res <- readRDS("/home/share/scratch/Clustering/without_1k/full_clus_res.rds")
  average_expr <- readRDS("/home/share/scratch/Clustering/without_1k/average_expr.rds")
  
}
if(with1k <- TRUE){
  # Significant <- readRDS("~/ses-1/Res/Clustering_aging/1k/1st/Significant.rds")
  Significant <- readRDS("~/ses-1/Res/Clustering_aging/1k_ancestry/Significant.rds")
  full_clus_res <- readRDS("~/ses-1/Res/Clustering_aging/1k_ancestry/full_clus_res.rds")
  average_expr <- readRDS("~/ses-1/Res/Clustering_aging/1k_ancestry/average_expr.rds")
  
}
average_expr = average_expr %>% map(~ .x %>% t %>% as.data.frame)

get_sig_clus = function(treatment, Significant, full_clus_res){
  sig_treatment = Significant %>% map(~ pluck(.x, treatment))
  no_clus = full_clus_res %>% map(max)
  full_clus = no_clus %>% map(~ seq(1:.x))
  sig_clus = map2(sig_treatment, full_clus, ~ (.y %in%.x))
  
}

pheno = pData(dat)
# plate has two unused levels, with 0 in these two groups, "PlateYear2Plate10" "PlateYear2Plate9"
# if not dropping them explicitly, when using MASS::polr causing a mismatched dimention problem, as 
# inside this function they drop it somewhere explicit and at some place using the original design matrix, 
# therefore causing the dimension mismatch.
pheno$Plate = droplevels(pheno$Plate) 
subData = pheno %>% dplyr::select(all_of(treatment), batch, all_of(controls$basic))

non_missing = complete.cases(subData)

dat@phenoData@data = pheno[non_missing,]

fit_mediate_cluster = function(treatment, gene_set_name, controls, out = NULL){
  out = mediators %>% 
    set_names() %>%
    map(safely(mediate_cluster), treatment, controls, gene_set_name)
}
mediate_cluster = function(mediator, treatment, controls, gene_set_name){
  clus = average_expr[[gene_set_name]] 
  colnames(clus) = str_c("d", 1:dim(clus)[2])
  clus = rownames_to_column(clus, var = "AID")

  datt_m  = pData(dat) %>% 
    dplyr::select(AID = AID, all_of(controls), all_of(treatment), all_of(mediator)) %>% 
    rename(treatment = treatment) 
  
  if(remove_diseased_subjects) {
    datt_m =
      datt_m %>% 
      remove_diseased_subjects_from_datt(gene_set_name, controls)
  }
  
  outcome = colnames(clus %>% dplyr::select(-AID)) %>% set_names() 
  
  convert_outcome_clus = function(outcome, datt_m){
    datt_pca = left_join(datt_m, clus %>% dplyr::select(AID, !!outcome)) %>% dplyr::select(-AID)
  }
  
  datt_pca = outcome %>% map(convert_outcome_clus, datt = datt_m)
  
  sig = get_sig_clus(treatment, Significant, full_clus_res) %>% pluck(gene_set_name)
  
  pmap(list(datt_pca[sig], outcome[sig], mediator), fit_m99) %>% map(extract_m99)
}
# debugonce(fit_mediate_cluster)

# plan(multicore, workers = 10)
example0 =
  args %>%
  filter(
    # is.element(gene_set_name, table1) &
    gene_set_name == "aging_mRNA",
    names(controls) == "ancestryPC") %>%
  # filter(gene_set_name=="Depression_mRNA") %>% 
  mutate(out = pmap(., safely(fit_mediate_cluster)),
         controls = names(controls))

example0 %>% saveRDS("./user_wx/mediate_cluster_sc5levels_allsample.rds")