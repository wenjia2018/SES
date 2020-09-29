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

walk(dir(path = here("R"),full.names = TRUE), source)
fit_m4 = partial(fit_m4, n_perm = 1000) # specify n_perm

load_data(reconciled = FALSE, remove_inflam = FALSE)
define_treatments_and_controls()
recode_variables_in_dat()
print(abbreviations)
funcs = str_subset(abbreviations$shorthand, "^m") 
funcs = funcs %>% str_subset("m[7-8]")
signatures =  add_de_novo(signatures, treatment)
de_novo_geneset = str_c(treatment, "_de_novo")


de_novo_treatment = c("ses_sss_composite", "income_hh_ff5", "sss_5")
de_novo_geneset = de_novo_treatment %>% str_c("_","de_novo")
# signatures = get_de_novo(de_novo_treatment, control = "all", remove_inflam = TRUE)
signatures = readRDS("./user_wx/signatures.rds")
ncomp = signatures$outcome_set[de_novo_geneset] %>% map_dfc(length) %>% unlist() %>% min


fit_pca_util = partial(fit_pca_util, ncomp = ncomp) # specify n_perm
gene_set_name = gene_set_name %>% append(de_novo_geneset) 
args = crossing(treatment, gene_set_name, controls)

# debugonce(model_fit)
example0 =
  args %>%
  filter(is.element(gene_set_name, de_novo_geneset),
         treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5"), 
         names(controls) == "all") %>% 
  mutate(out = pmap(., safely(model_fit), funcs),
         controls = names(controls))



example0 %>% saveRDS("/home/share/scratch/example0_without_1KI_de_novo.rds")