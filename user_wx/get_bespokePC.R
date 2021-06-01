library(here)
library(tidyverse)
library(rlang)
library(skimr)
library(furrr)
library(tidyr)
source("./R/utils_snpPC.R")
# eqtl_pca <- function(gene_set_name, p_eqtl) {
# 
#   if (gene_set_name=="whole_genome"){
#     # for the whole genome
#     dat <- readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_17.11.2020.rds")
#     gene_sets <- Biobase::featureNames(dat)
#   }else{
#     # for a predefined signature 
#     signatures <- readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_22.03.2021_signature.rds")
#     gene_sets <- signatures$outcome_set %>% pluck(gene_set_name)
#   }
# 
#   G_list <- readRDS("/home/share/data_input/genename_15062020.rds")
#   # specific gens of our interest
#   gene <- G_list %>%
#     filter(hgnc_symbol %in% gene_sets) %>%
#     pull(ensembl_gene_id)
#   # eqtl catelog
#   catelog <- readRDS("/home/share/dna_ancestry/eqtl/whole.blood.rds")
#   # filter Snps according to some p value threshold p_eqtl
#   catelog <- catelog[catelog$pvalue < p_eqtl, ]
#   catelog <- catelog %>% filter(gene_id %in% gene)
#   
#   # save Snps which are assoicated with specific genes passing p_eqtl
#   data.table::fwrite(unique(catelog$rsid) %>% as.data.frame(),
#                      file = str_c("/home/share/dna_ancestry/dna/keep", gene_set_name, "_", p_eqtl, ".txt"),
#                      col.names = F
#   )
#   
#   # keep the snps from data set
#   # system("plink --bfile ./dna/omni_joined.freeze3.sharedMarkers --extract ./dna/keep.txt --make-bed --out ./dna/kept.omni")
#   system(str_c("plink --bfile /home/share/dna_ancestry/dna/omni_joined.freeze3.sharedMarkers --extract /home/share/dna_ancestry/dna/keep", gene_set_name, "_", p_eqtl, ".txt --make-bed --out /home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl, ".omni"))
#   # run PCAs from plink file
#   # system("plink --bfile ./dna/kept.omni --pca --out ./dna/kept.omni.pca")
#   system(str_c("plink --bfile /home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl, ".omni --pca --out /home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl, ".omni.pca"))
#   # output will contain a kept......pca.eigenvec file with family ID, AID, then 20 PCs
#   
#   eigenvalue <- data.table::fread(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl, ".omni.pca.eigenval"))
#   ancestryPC_dim = dim(eigenvalue)[1]
#   ancestryPC_names <- str_c("AncestryPC", c(1:ancestryPC_dim))
#   
#   custom_PCA = read.table(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl,".omni.pca.eigenvec"),
#                           col.names = c("fid", "AID", ancestryPC_names))
#   
#   
#   custom_PCA %>% saveRDS(str_c("/home/share/dna_ancestry/dna/custom_PCA_", gene_set_name, "_", p_eqtl,".rds"))
# }

# examples to get bespoke pca for new signatures and the whole genome:
if(0){
  gene_set_name = c(
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
  
  p_eqtl <- c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
  
  args =
    crossing(gene_set_name, p_eqtl)
  args %>%
    pmap(eqtl_pca)
}

if(0){
  gene_set_name = "whole_genome"
  p_eqtl <- c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
  args =
    crossing(gene_set_name, p_eqtl)
  plan(multicore, workers = 10)
  args %>%
    furrr::future_pmap(eqtl_pca)
  
  gene_set_name = "aging_cluster_complement_mRNA"
  p_eqtl <- c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
  args =
    crossing(gene_set_name, p_eqtl)
  args %>%
    pmap(eqtl_pca)
  plan(multicore, workers = 10)
  args %>%
    furrr::future_pmap(eqtl_pca)
}

# get bespoke pca for new signatures:

# gene_set_name = c(
#   "aging_down_cl1_mRNA",
#   "aging_down_cl1a_mRNA",   
#   "aging_down_cl1b_mRNA",
#   "aging_down_cl1c_mRNA",
#   "aging_down_cl2_mRNA",
#   "aging_down_cl3_mRNA",
#   "aging_up_cl1_mRNA",
#   "aging_up_cl2_mRNA",
#   "aging_up_cl3_mRNA",
#   "aging_up_cl4_mRNA"
# )
# 
# p_eqtl <- c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
# 
# args =
#   crossing(gene_set_name, p_eqtl)
# 
# plan(multicore, workers = 25)
# args %>% 
#   furrr::future_pmap(eqtl_pca)

# some cases where ancestry PC is less than default 20, has already revised the code to fix this bug
# gene_set_name = "aging_up_cl3_mRNA"
# p_eqtl <- c(1e-7, 1e-8, 1e-9, 1e-10)
# args =
#   crossing(gene_set_name, p_eqtl)
# 
# plan(multicore, workers = 4)
# args %>%
#   furrr::future_pmap(eqtl_pca)

if(1){
  gene_set_name = c(
    "CVD_mRNA",
    "diabetes_mRNA",
    "Rheumatoid_Arthritis_mRNA", 
    "Alzheimers_mRNA",
    "Aortic_Aneurysm_mRNA",
    "COPD_mRNA",
    "Asthma_mRNA",
    "Hypertension_mRNA",
    "Depression_mRNA",
    "CKD_mRNA"
  )
  p_eqtl <- c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
  args =
    crossing(gene_set_name, p_eqtl)
  plan(multicore, workers = 40)
  args %>%
    furrr::future_pmap(eqtl_pca)
  
}


