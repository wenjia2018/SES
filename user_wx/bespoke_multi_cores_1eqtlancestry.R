# run analysis controlling for eqtl bespoke ancestry PCs
# eqtl associated snps are obtained given a specific gene set and a threshold p value
# !!! attention both define_treatment_controls_outcomes (define the treatment and control sets you want for the analysis)
# and recode_variables_in_dat(recode categorical variables with less levles for example) 
# has to be changed according to specific analysis!!!
# better check and redefine these two files for your specific analysis  
set.seed(123)
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
fit_m4 <- partial(fit_m4, n_perm = 1000) # specify n_perm

# explicitly assign ncomp as the smallest number of table signatures gene numbers
# attention for gene set smaller than 10
ncomp <- 10
fit_pca_util <- partial(fit_pca_util, ncomp = ncomp) # specify n_perm
print(abbreviations)
# linearhpo <- TRUE
if(FALSE) {
  ftest_v <- c("raceethnicity__color_byinterviewer3_NonHwhite|DarkBlack", "raceethnicity__color_byinterviewer3_NonHwhite|LightMed")
  
}
mediation_mean = FALSE
mediation_each_gene = FALSE
funcs <- str_subset(abbreviations$shorthand, "^m")
# select models to run
funcs <- c("m7", "m8", "m10" , "m11")
# funcs <- c( "m8")
# funcs <- c("m1","m2", "m3")
# funcs <- c( "m96")
# funcs <- NULL
# load snp related to skin color in case you want to perform snp regression
dt_color_snp <- readRDS("/home/share/dna_ancestry/dna/dt_color_snp.rds")
############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################

fit_bespoke <- function(gene_set_name, p_eqtl) {
  load_data(reconciled = FALSE, remove_inflam = FALSE)
  signatures$outcome_set$wholegenome <<-  featureNames(dat)
  # find the important PCs to be included for each analysis
  ancestryPC <- get_PC_dim("aging_mRNA", p_eqtl)
  define_treatments_and_controls_bespoke(gene_set_name, ancestryPC)
  # if(linearhpo){
  #   ftest_v <<- str_subset(c(treatment) %>% unique, "__") 
  # }
  
  # define_treatments_and_controls_snps(gene_set_name, ancestryPC)
  # custom_PCA <- readRDS(str_c("/home/share/dna_ancestry/dna/custom_PCA_", gene_set_name, "_", p_eqtl, ".rds"))
  custom_PCA <- readRDS(str_c("/home/share/dna_ancestry/dna/custom_PCA_", "aging_mRNA", "_", p_eqtl, ".rds"))
  custom_PCA <- custom_PCA %>%
    select(-fid) %>%
    mutate(AID = AID %>% as.character())
  recode_variables_in_dat_bespoke(custom_PCA)

    if(gene_set_name == "whole_genome"){
      example0 = NULL
      
      example1 <-
        args %>%
        filter(str_detect(names(controls), "ancestryPC_ses")) %>%
        mutate(gene_set_name = "whole_genome") %>%
        mutate(
          out = pmap(., safely(model_fit), funcs),
          control_set = names(controls)
        )
    }else{
      example0 <-
        args %>%
        filter(str_detect(names(controls), "ancestryPC_ses")) %>%
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
# table1 =
#   c(
#     "CVD_mRNA",
#     "inflam1k_mRNA",
#     "diabetes_mRNA",
#     "Rheumatoid_Arthritis_mRNA", 
#     "Alzheimers_mRNA",
#     "Aortic_Aneurysm_mRNA",
#     "COPD_mRNA",
#     "Asthma_mRNA",
#     "Hypertension_mRNA",
#     "Depression_mRNA",
#     "CKD_mRNA",
#     "whole_genome"
#   )
table1 <-
  c(
    # "ctra_mRNA",
    # "inflame_mRNA",
    # "interferon_mRNA",
    # "AntBIntF_mRNA",
    # "antibody_mRNA", #only 1 gene
    # "inflam1k_mRNA",
    "aging_mRNA",
    "aging_up_mRNA",
    "aging_down_mRNA",
    "aging_cluster_complement_mRNA",
    # "aging_down_cl1_mRNA",
    "aging_down_cl1a_mRNA",
    "aging_down_cl1b_mRNA",
    "aging_down_cl1c_mRNA",
    "aging_down_cl2_mRNA",
    "aging_down_cl3_mRNA",
    "aging_up_cl1_mRNA",
    "aging_up_cl2_mRNA",
    "aging_up_cl3_mRNA",
    "aging_up_cl4_mRNA",
    "whole_genome"
   # "wholegenome"# for whole genowide mediation 
  )
# table1 <-    c("aging_up_cl3_mRNA","aging_up_cl4_mRNA")
table1 <-"aging_cluster_complement_mRNA"
# table1 <-
#     c(
#       "darkblack005_mRNA",
#       "intersec001_mRNA",
#       "darkblack005intersec001_mRNA")
p_eqtl <- c(0.05)
# , 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)

args_eqtl <- crossing(table1, p_eqtl)
plan(multicore, workers = 14)
# debugonce(fit_bespoke)
# example_bespoke <- args_eqtl %>% mutate(out = pmap(list(gene_set_name = table1, p_eqtl = p_eqtl), safely(fit_bespoke)))

example_bespoke <- args_eqtl %>%
  mutate(out = furrr::future_pmap(list(gene_set_name = table1, p_eqtl = p_eqtl), safely(fit_bespoke)))

example_bespoke %>% saveRDS("./user_wx/skincolor_celltype_NonHwhite_eqtl005_aging_composite_ancestry_20.05.2021.rds")
# example_bespoke %>% saveRDS("./user_wx/bespoke_snps.rds")
# example_bespoke %>% saveRDS("./user_wx/bespoke_v4.rds")
# v2 :only ancestry controls
# v3 : all 3 sets of controls
# v4: corrected kept dim first row not as colnames
# v2 v3 and v4 results are same except mediation as it is based on bootstraping estimation
