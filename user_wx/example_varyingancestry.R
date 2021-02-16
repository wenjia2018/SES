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
library(MendelianRandomization)
walk(dir(path = here("R"),full.names = TRUE), source)
fit_m4 = partial(fit_m4, n_perm = 1000) # specify n_perm

# explicitly assign ncomp as the smallest number of table signatures gene numbers
# attention for gene set smaller than 10 
ncomp = 10 
fit_pca_util = partial(fit_pca_util, ncomp = ncomp) # specify n_perm
############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################
# load custom ancestry PC


load_data(reconciled = FALSE, remove_inflam = FALSE)

fit_bespoke = function(gene_set_name, p_eqtl){

  eigenvalue <- data.table::fread(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl,".omni.pca.eigenval"))
  eigenvec <- data.table::fread(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl,".omni.pca.eigenvec"))
  kept <- data.table::fread(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl,".omni.bim"))
  
  ag = PCDimension::AuerGervini(eigenvalue$V1, dd=c(dim(eigenvec)[1], dim(kept)[1]))
  d = PCDimension::agDimension(ag)
  ancestryPC_full = dput(str_c("AncestryPC", c(1:20)))
  ancestryPC = ancestryPC_full[1:d]
  
define_treatments_and_controls_bespoke(ancestryPC)

custom_PCA <- readRDS(str_c("/home/share/dna_ancestry/dna/custom_PCA_", gene_set_name, "_", p_eqtl,".rds"))
custom_PCA = custom_PCA %>% select(-fid) %>% mutate(AID = AID %>% as.character())
recode_variables_in_dat_bespoke(custom_PCA)
print(abbreviations)
funcs = str_subset(abbreviations$shorthand, "^m") 
funcs = funcs %>% str_subset("m[7-8]")

# debugonce(model_fit)


example0 =
  args %>% 
  filter(names(controls) =="ancestryPC") %>% 
  mutate(out = pmap(., safely(model_fit), funcs),
         control_set = names(controls))

example1 =
  args %>% 
  filter(names(controls) =="ancestryPC") %>% 
  mutate(gene_set_name = "whole_genome_and_tfbm") %>%
  mutate(out = pmap(., safely(model_fit), funcs),
         control_set = names(controls))

return(list(example0 = example0, example1 = example1))
}

gene_set_name =
  c(
    "ctra_mRNA",
    "inflame_mRNA",
    "interferon_mRNA",
    "AntBIntF_mRNA", 
    # "antibody_mRNA", #only 1 gene
    "inflam1k_mRNA",
    "aging_mRNA",
    "aging_up_mRNA",
    "aging_down_mRNA"
  )

p_eqtl = c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)

args = crossing(gene_set_name, p_eqtl)
# debugonce(fit_bespoke)
example_bespoke = args %>% 
  # slice(70) %>%
  mutate(out = pmap(., safely(fit_bespoke)))

example_bespoke %>% saveRDS("./user_wx/bespoke.rds")
