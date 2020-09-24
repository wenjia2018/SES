
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


############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################

load_data(reconciled = FALSE, remove_inflam = TRUE)
define_treatments_and_controls()
recode_variables_in_dat()
print(abbreviations)

funcs = "m7"
# explicitly assign ncomp as the smallest number of table signatures gene numbers

ncomp = signatures$outcome_set[table1]%>% map_dfc(length) %>% unlist() %>%  min
ncomp = 10
fit_pca_util = partial(fit_pca_util, ncomp = ncomp) # specify n_perm


# debugonce(model_fit)
example4 =
  args %>% 
  filter(treatment == "ses_sss_composite",
         names(controls) == "all",
         gene_set_name %in% c("ses_de_up_mRNA","ses_de_down_mRNA")) %>%
  mutate(out = pmap(., safely(model_fit), funcs),
         controls = names(controls))

example4 %>% saveRDS("/home/share/scratch/example4.rds")