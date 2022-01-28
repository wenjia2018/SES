# main function to start running analysis
# for skin color with bespoke ancestry controls
# load packages
library(here)
# library(data.table)
# setDTthreads(threads = 20)
library(tidyverse)
library(EValue)
library(rlang)
library(skimr)
library(furrr)
library(limma)
# library(recipes)
# recipes and Evalue has conflict, after loading recipes,
# evalue package doesnot work with the error Error: $ operator is invalid for atomic vectors
library(parsnip)
library(workflows)
library(Biobase)
library(dbr) # my package

# source all the useful files
walk(dir(path = here("R"),full.names = TRUE), source)

############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################
# set some parameters values
# which PCA to perform
oblimin = FALSE
nn = TRUE
# explicitly assign ncomp for PCA analysis
ncomp = 10
fit_pca_util = partial(fit_pca_util, ncomp = ncomp) # specify n_perm
# type of mediation
mediation_mean = FALSE
mediation_each_gene = FALSE
# for doing genowide DE analysis only
normalization_bydesign = FALSE
# specify if subjects with disease shall be removed
remove_diseased_subjects = TRUE

funcs <- c("m1", "m2","m3", "m7", "m8", "m10")
funcs = "m13"
if(funcs == "m13"){
  boot = TRUE
  N = 10000

}else{
  boot = FALSE
}

fit_bespoke <- function(gene_set_name, p_eqtl) {
  load_data(reconciled = FALSE, remove_inflam = FALSE)
  
  # find the important PCs to be included for each analysis
  ancestryPC <- get_PC_dim(gene_set_name, p_eqtl)
  define_treatments_and_controls_sc_bespoke(gene_set_name, ancestryPC)
  custom_PCA <- readRDS(str_c("/home/share/dna_ancestry/dna/custom_PCA_", gene_set_name, "_", p_eqtl, ".rds"))
  custom_PCA <- custom_PCA %>%
    select(-fid) %>%
    mutate(AID = AID %>% as.character())
  recode_variables_in_dat_racedummy_bespoke(custom_PCA)
  if(gene_set_name == "whole_genome_and_tfbm"){
    example0 = NULL
    
    example1 <-
      args %>%
      filter(str_detect(names(controls), "ancestryPC")) %>%
      mutate(gene_set_name = "whole_genome_and_tfbm") %>%
      mutate(
        out = pmap(., safely(model_fit), funcs),
        control_set = names(controls)
      )
  }else{
    example0 <-
      args %>%
      filter(str_detect(names(controls), "ancestryPC")) %>%
      mutate(
        out = pmap(., safely(model_fit), funcs),
        control_set = names(controls))
    
    example1 <- NULL
  }
  
  
  return(list(
    example0 = example0,
    example1 = example1
  ))
}
# TABLE 1
table1 =
  c(
    # 'whole_genome_and_tfbm'
    # "CVD_mRNA",
    # "diabetes_mRNA",
    # "inflam1k_mRNA",
    # "Rheumatoid_Arthritis_mRNA", "Alzheimers_mRNA",
    # "Aortic_Aneurysm_mRNA", "COPD_mRNA",
    # "Asthma_mRNA","Hypertension_mRNA",
    # "Depression_mRNA",
    # "CKD_mRNA"
    # "ctra_mRNA",
    # "inflame_mRNA",
    # "interferon_mRNA",
    # "AntBIntF_mRNA",
    # "antibody_mRNA", #only 1 gene
    "aging_mRNA",
    "aging_up_mRNA",
    "aging_down_mRNA",
    "aging_cluster_complement_mRNA",
    "aging_down_cl1_mRNA",
    "aging_down_cl1a_mRNA",
    "aging_down_cl1b_mRNA",
    "aging_down_cl1c_mRNA",
    "aging_down_cl2_mRNA",
    "aging_down_cl3_mRNA",
    "aging_up_cl1_mRNA",
    "aging_up_cl2_mRNA",
    "aging_up_cl3_mRNA",
    "aging_up_cl4_mRNA"
  )

p_eqtl <- c(0.05)
# , 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)

args_eqtl <- crossing(table1, p_eqtl)
plan(multicore, workers = 14)
# debugonce(fit_bespoke)
example_bespoke <- args_eqtl %>% mutate(out = furrr::future_pmap(list(gene_set_name = table1, p_eqtl = p_eqtl), safely(fit_bespoke)))
example_bespoke %>% saveRDS("./user_wx/m13_aging_with1k_sc5levels_allsample_bespoke_bootxx.rds")
