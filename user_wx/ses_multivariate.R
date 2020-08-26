library(tidyverse)
library(limma)
library(tidymodels)
library(Biobase) 
library(dbr) # my package
source("R/load_data.R")
source("R/fit_multivariate_model.R")

############################################################
# LOAD DATA (reconciled = ideal setting)
############################################################

load_data(reconciled = TRUE) %>% 
  list2env(.GlobalEnv) 

############################################################
# DEFINE "TREATMENTS", CONTROLS, OUTCOMES (INCLUDING TFBMs) 
############################################################

# TABLE 1
table1 =
  c(
    "CVD_mRNA",
    "diabetes_mRNA",
    "inflam1k_mRNA", "breast_cancer_mRNA",
    "Lupus_mRNA", "Colorectal_mRNA",
    "Rheumatoid_Arthritis_mRNA", "Alzheimers_mRNA",
    "Aortic_Aneurysm_mRNA", "COPD_mRNA",
    "Asthma_mRNA","Hypertension_mRNA",
    "kidney_transplant_tolerance_mRNA"
  )

treatment = c("ses_sss_composite", "sss_5", "SEI_ff5", "edu_max", "income_hh_ff5")

controls = 
  list(
    basic = 
      c(
        "sex_interv", "re", "Plate", "AvgCorrelogram100" ,"age_w5",
        "BirthY", "W5REGION", "pregnant_biow5", 
        "kit_biow5", "tube_biow5",  "FastHrs",
        "travel_biow5",  "months_biow5", "time_biow5"
      ),
    biological = 
      c( 
        "B.cells.naive", "B.cells.memory", "Plasma.cells",
        "T.cells.CD8", "T.cells.CD4.naive", "T.cells.CD4.memory.resting",
        "T.cells.CD4.memory.activated",
        # "T.cells.follicular.helper",
        "T.cells.regulatory..Tregs.", "T.cells.gamma.delta",
        "NK.cells.resting", "NK.cells.activated", "Monocytes", "Macrophages.M0", 
        # "Macrophages.M1",
        "Macrophages.M2", "Dendritic.cells.resting",
        "Dendritic.cells.activated", "Mast.cells.resting",
        # "Mast.cells.activated", # not estimable in limma
        "Eosinophils", "Neutrophils" 
      ) 
  ) %>% 
  c(all = list(unique(unlist(.)))) 

gene_set_name = signature_names %>% append("whole_genome_and_tfbm") 

args = crossing(treatment, gene_set_name, controls)

immune_tfbms = 
  c(
    "CEBPG_CREB3L1", "CREB3", "CREB3L1", "IRF2", "IRF3",
    "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "JUN", "NFKB1", 
    "NFKB2", "NR3C1"
  )

############################################################
# EXAMPLE 1: ESTIMATE WHOLE GENOME DE AND TFBM
############################################################

example1 = 
  args %>% 
  filter(treatment == "ses_sss_composite", 
         gene_set_name == "whole_genome_and_tfbm", 
         names(controls) == "basic") %>% 
  get_results(permutation = FALSE)

# INSPECT DE
example1 %>% pluck("out", "result", "ttT")  
example1 %>% pluck("out", "result", "ttT") %>%  filter(adj.P.Val < 0.05) %>% pluck("gene") # which genes have corrected p-value < 0.05

# INSPECT TFBM: uncorrected p values
example1 %>% pluck("out", "result", "tfbm")  # TFBM de Novo
example1 %>% pluck("out", "result", "tfbm")  %>% filter(tfbm %in% immune_tfbms) # Inflammation

############################################################
# EXAMPLE 2: ESTIMATE UNIVARIATE/MULTIVARIATE MODELS FOR TABLE 1
############################################################

example2 = 
  args %>% 
  filter(treatment == "ses_sss_composite", 
         is.element(gene_set_name, table1), 
         names(controls) == "basic") %>% 
  get_results(permutation = FALSE) 

# INSPECT
# what do the column abbreviations of the following table mean?
example2 %>% pluck("out", "result", "table")
summaries = c("m1_detail")
summaries = c("m1_p", "m2_p" , "m3_p")
example2 %>% 
  filter(out_id == "result") %>%
  unnest(out) %>% 
  unnest(all_of(summaries)) %>% 
  select(treatment, controls, gene_set_name, all_of(summaries))

############################################################
# diagnose errors in any model?
############################################################

(errors =
    example2 %>% 
    filter(out_id == "error") %>%
    filter(map_lgl(out, negate(is.null)))) 

errors %>% pluck("out")

debugonce(fit_multivariate_model)
errors %>% 
  select(treatment, gene_set_name, controls) %>% 
  get_results(permutation = FALSE)

############################################################
# EXAMPLE 3: SIMILAR TO EXAMPLE 2 BUT WITH ALL CONTROLS AND CONSTITUENT TREATMENTS
############################################################

example3 = 
  args %>% 
  filter(is.element(gene_set_name, table1)) %>% 
  get_results(permutation = FALSE) 
