set.seed(123)

library(here)

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
source("./dda_v0.1/dda_resdist.r")
source("./dda_v0.1/dda_vardist.r")
source("./dda_v0.1/dda_indep.r")
source("./dda_v0.1/boot_hsic_test.R")
source("./dda_v0.1/nlcor_test.r")
n_boot = 5000


############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################
# choose normalization methods for downstream analysis
tmm = TRUE
rle = FALSE
log2cpm = FALSE
# which PCA to perform
oblimin = TRUE
nn = FALSE
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
funcs = str_subset(abbreviations$shorthand, "^m") 
funcs <- c("m1", "m2","m3", "m7", "m8", "m10")
funcs = "m7"

fit_pca_util = partial(fit_pca_util, ncomp = ncomp) # specify n_perm
# creat de novo
sig = Reduce(union, signatures$outcome_set[table1[1:11]])
example_TMM_DE <- readRDS("~/ses-1/user_wx/example_tmm_genowidebydesign.rds")
# parental SES
example_TMM_DE <- readRDS("~/ses-1/user_wx/example_parentlDE.rds")

data_DE = example_TMM_DE %>% 
  hoist(out, ttT = list("result", "ttT")) %>% 
  filter(map_lgl(.$ttT, ~dim(.)[1]!=0)) %>% 
  mutate(gene_sig = map(ttT, ~ filter(., adj.P.Val<0.05) %>% pull(gene)) %>% set_names(treatment)) %>% 
  dplyr::select(-out, -gene_set_name) 

data_DE_without1k = data_DE %>% 
  mutate(ttT = ttT %>% map(~ filter(.x, !(gene %in% sig)))) %>% 
  mutate(gene_sig = map(ttT, ~ filter(., adj.P.Val<0.05) %>% pull(gene)) %>% set_names(treatment)) 


DE_with1k = data_DE$gene_sig %>% 
  set_names(data_DE$treatment) %>% 
  map(~ intersect(.x, featureNames(dat)))
  


DE_remove1k = map(DE_with1k, ~ setdiff(.x, sig))

ses4_unique = Reduce(setdiff, list(DE_remove1k$ses_sss_composite,DE_remove1k$income_hh_ff5, DE_remove1k$sss_5, DE_remove1k$edu_max, DE_remove1k$SEI_ff5))
income_unique = Reduce(setdiff, list(DE_remove1k$income_hh_ff5, DE_remove1k$ses_sss_composite, DE_remove1k$sss_5, DE_remove1k$edu_max, DE_remove1k$SEI_ff5))

de_novo_signatures = list(
  # income_unique_de_novo = income_unique,
  # ses4_unique_de_novo = ses4_unique,
  # ses4income_de_novo = intersect(DE_remove1k$income_hh_ff5, DE_remove1k$ses_sss_composite),
  # ses4income_unique_de_novo = intersect(DE_remove1k$income_hh_ff5, DE_remove1k$ses_sss_composite) %>% setdiff(Reduce(union, list(DE_remove1k$sss_5, DE_remove1k$edu_max, DE_remove1k$SEI_ff5))),
  # sss_de_novo = DE_remove1k$sss_5,
  # income_de_novo = DE_remove1k$income_hh_ff5,
  # edu_de_novo = DE_remove1k$edu_max,
  # sei_de_novo = DE_remove1k$SEI_ff5,
  # ses_union_de_novo = Reduce(union,list(DE_remove1k$income_hh_ff5, DE_remove1k$ses_sss_composite, DE_remove1k$sss_5, DE_remove1k$edu_max, DE_remove1k$SEI_ff5)),
  ses4_de_novo = DE_remove1k$ses_sss_composite,
  ses4_with1k_de_novo = DE_with1k$ses_sss_composite
  ) 
signatures$outcome_set = c(signatures$outcome_set, de_novo_signatures)
de_novo_geneset = de_novo_signatures %>% names
gene_set_name = de_novo_geneset
args = crossing(treatment, gene_set_name, controls)
# debugonce(model_fit)

example0 =
  args %>%
  filter(
    # is.element(gene_set_name, de_novo_geneset) &
      # gene_set_name== "Aortic_Aneurysm_mRNA" &
      # treatment == "ses_sss_composite" &
      names(controls) == "basic") %>%
  mutate(out = pmap(., safely(model_fit), funcs),
         controls = names(controls))

# With controls used in SES paper: predicts the signatures from SES paper + Peters aging signature. 
example0 %>% saveRDS("./user_wx/example_tmm_m7_denovosesDE.rds") 

