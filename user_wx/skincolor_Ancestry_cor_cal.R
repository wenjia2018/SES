# correlation of skincolor with ancestry PC
set.seed(123)

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
# recipes and Evalue has conflict, after loading recipes, evalue package doesnot work with the error Error: $ operator is invalid for atomic vectors
library(parsnip)
library(workflows)
library(Biobase)
library(dbr) # my package

walk(dir(path = here("R"),full.names = TRUE), source)
fit_m4 = partial(fit_m4, n_perm = 1000) # specify n_perm
# source("./dda_v0.1/dda_resdist.r")
# source("./dda_v0.1/dda_vardist.r")
# source("./dda_v0.1/dda_indep.r")
# source("./dda_v0.1/boot_hsic_test.R")
# source("./dda_v0.1/nlcor_test.r")
n_boot = 5000


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

mediation_mean = FALSE
mediation_each_gene = FALSE
funcs = str_subset(abbreviations$shorthand, "^m") 
funcs <- c("m1", "m2","m3", "m7", "m8", "m10")
funcs = "m7"

fit_pca_util = partial(fit_pca_util, ncomp = ncomp) # specify n_perm

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

corr_fun = function(gene_set_name, p_eqtl) {
  load_data(reconciled = FALSE, remove_inflam = FALSE)
  
  # find the important PCs to be included for each analysis
  ancestryPC_keep <- get_PC_dim(gene_set_name, p_eqtl)
  ancestryPC <- str_c("AncestryPC", c(1:20))
  define_treatments_and_controls_sc_bespoke(gene_set_name, ancestryPC)
  custom_PCA <- readRDS(str_c("/home/share/dna_ancestry/dna/custom_PCA_", gene_set_name, "_", p_eqtl, ".rds"))
  custom_PCA <- custom_PCA %>%
    select(-fid) %>%
    mutate(AID = AID %>% as.character())
  recode_variables_in_dat_racedummy_bespoke(custom_PCA)
  
  pheno = pData(dat) %>% dplyr::select(AID = AID, starts_with("AncestryPC"), all_of(treatment), color_byinterviewer5) 
  pheno$color_byinterviewer5 = factor(pheno$color_byinterviewer5, levels = c("White","Light",  "Medium", "Dark", "Black"))
  # brant::brant() to test proportional odds assumption in order to use ordinal logistic regression
  model_full5_ordered = MASS::polr(str_c("color_byinterviewer5 ~",ancestryPC %>% str_c(collapse = " + ")) , data = pheno, Hess=TRUE)
  parallel_test = brant::brant(model_full5_ordered)
  temp = car::Anova(model_full5_ordered)
  anova_ordinal = temp %>% as.data.frame() %>%  rownames_to_column(var = "AncestryPC") %>% filter(`Pr(>Chisq)`<0.05)
  
  model_full5_multinomial = VGAM::vglm(str_c("color_byinterviewer5 ~",ancestryPC %>% str_c(collapse = " + ")), family = VGAM::multinomial(refLevel="White"), model=TRUE, data=pheno)
  temp1 = car::Anova(model_full5_multinomial)
  anova_unordered = temp1 %>% as.data.frame() %>%  rownames_to_column(var = "AncestryPC") %>% filter(`Pr(>Chisq)`<0.05)
  
  cor_results = funcs::cor2(pheno %>% select(color_byinterviewer5, ancestryPC_keep)) 
  
  return(list(anova_ordinal = anova_ordinal, anova_unordered = anova_unordered, cor_results = cor_results, parallel_test = parallel_test))
}
p_eqtl = 0.05
args_eqtl <- crossing(table1, p_eqtl)
plan(multicore, workers = 14)
# debugonce(corr_fun)
example <- args_eqtl %>% mutate(out = furrr::future_pmap(list(gene_set_name = table1, p_eqtl = p_eqtl), safely(corr_fun)))
example %>% saveRDS("./user_wx/color_ancestry_cor.rds")

