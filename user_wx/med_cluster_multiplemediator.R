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

# 
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
controls = controls$basic
pheno = pData(dat)
pheno$Plate = droplevels(pheno$Plate) 
subData = pheno %>% dplyr::select(all_of(treatment), batch, all_of(controls))

non_missing = complete.cases(subData)

dat@phenoData@data = pheno[non_missing,]

fit_mediate_cluster = function(treatment, gene_set_name, controls, out = NULL){
  out = mediate_cluster_multiple(treatment, controls, gene_set_name)
}
mediate_cluster_multiple = function(treatment, controls, gene_set_name){
  print("*******************")
  print(treatment)
  print(gene_set_name)
  print(controls)
  print("*******************") 
  clus = average_expr[[gene_set_name]] 
  colnames(clus) = str_c("d", 1:dim(clus)[2])
  clus = rownames_to_column(clus, var = "AID")
  datt_m  = pData(dat) %>% 
    dplyr::select(AID = AID, all_of(controls), all_of(treatment), all_of(mediators)) %>% 
    rename(treatment = treatment) %>% 
    select(all_of(mediators), everything())
  
  if(remove_diseased_subjects) {
    datt_m =
      datt_m %>% 
      remove_diseased_subjects_from_datt(gene_set_name, controls)
  }
  
  outcome = colnames(clus %>% select(-AID)) %>% set_names() 
  
  convert_outcome_clus = function(outcome, datt_m){
    datt_pca = left_join(datt_m, clus %>% select(AID, !!outcome)) %>% select(-AID)
  }
  
  datt_pca = outcome %>% map(convert_outcome_clus, datt = datt_m)
  
  sig = get_sig_clus(treatment, Significant, full_clus_res) %>% pluck(gene_set_name)
  
  pmap(list(datt_pca[sig], outcome[sig]), mediate_multimediate) %>% map(~.x %>% summary(opt = "avg"))
}

if(1) {
  resetlevel = function(x) {
    levels(x) =  c(0:(length(levels(x))-1))
    x
  }
  temp_data = pData(dat)
  temp_data <- temp_data %>%
    # mutate_at(
    #   .vars = vars(matches("^edu_p$|^edu_max$")),
    #   .funs = list(~ .x %>%
    #                  fct_recode("0" = "high or less",
    #                             "1" = "more than high")
    #   ))  %>% 
    # mutate_at(
    #   .vars = vars(c("sex_interv")),
    #   .funs = list(~ .x %>%
    #                  fct_recode("0" = "f",
    #                             "1" = "m")
    #   ))  %>% 
    # mutate_at(
    #   .vars = vars(c("re")),
    #   .funs = list(~ .x %>%
    #                  fct_recode("0" = "1",
    #                             "1" = "2",
    #                             "2" = "3",
    #                             "3" = "4",
    #                             "4" = "5")
    #   ))   %>% 
    mutate_at(
      .vars = vars(c("BirthY", "time_biow5")),
      .funs = list(~ .x %>% as.factor()
      ))  %>% 
    mutate_at(
      .vars = vars(c("sex_interv", "re", "Plate", "BirthY", "W5REGION", "pregnant_biow5", 
                     "kit_biow5", "tube_biow5", "FastHrs",              
                     "travel_biow5", "months_biow5", "time_biow5", "edu_max",              
                     "bills_binary", "currentsmoke_binary", "insurance_lack_binary")),
      .funs = list(~  .x %>% resetlevel()
      )) 
  
  dat@phenoData@data = temp_data
  # 
  # levels(temp_data$Plate) = c(0:(length(levels(temp_data$Plate))-1))
  # levels(temp_data$BirthY) = c(0:(length(levels(temp_data$BirthY))-1))
  # levels(temp_data$W5REGION) = c(0:(length(levels(temp_data$W5REGION))-1))
  # levels(temp_data$pregnant_biow5) = c(0:(length(levels(temp_data$pregnant_biow5))-1))
  # levels(temp_data$kit_biow5) = c(0:(length(levels(temp_data$kit_biow5))-1))
  # levels(temp_data$tube_biow5) = c(0:(length(levels(temp_data$tube_biow5))-1))
  # levels(temp_data$FastHrs) = c(0:(length(levels(temp_data$FastHrs))-1))
  # 
  # 
  
}
# debugonce(fit_mediate_cluster)
COR = TRUE
# plan(multicore, workers = 10)
example0 =
  args %>%
  filter(
    is.element(gene_set_name, table1) &
      # gene_set_name == "SES  Composite",
      names(controls) == "basic") %>%
  # filter(gene_set_name=="Depression_mRNA") %>% 
  mutate(out = pmap(., safely(fit_mediate_cluster)),
         controls = names(controls))

example0 %>% saveRDS("./user_wx/mediate_cluster_multiple_with1k_withCOR.rds")