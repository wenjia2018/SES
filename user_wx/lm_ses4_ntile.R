# 
# In most cases, category 1 and 2 are very similar; there is a drop between 2 and 3; and 3, 4, 5 are not that different.
# 
# I would like to test this idea, please. 
# 
# A. For each of these 21 genes, could you run a regression equation for each please. 
# That is 21 regressions. Standard controls, including age and birth cohort. 
library(tidyverse)
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase)
library(ggformula)
library(here)
library(dbr)
library(venn)
walk(dir(path = here("R"),full.names = TRUE), source)
############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################

load_data(reconciled = FALSE, remove_inflam = FALSE)
define_treatments_and_controls()
recode_variables_in_dat()
treatment = "ses_sss_composite"
control = controls$all
gene_list = readxl::read_excel("./user_wx/ses4uniquedown_mean_ses4group.xlsx")$gene
# Specify whole-genome regression of rna on design
y <- dat[gene_list] %>% Biobase::exprs() %>%  t %>% as.data.frame() %>%  rownames_to_column(var = "AID") %>%
  mutate(AID = substring(AID, 2, 100) %>% as.character())
X <- dat %>% phenoData() %>% .@data %>% dplyr::select(all_of(control), all_of(treatment)) %>%
  # mutate(group = ses_sss_composite %>% ntile(5) %>% as.factor %>% relevel(ref = 4)) %>% 
  mutate(ses_sss_composite = ses_sss_composite %>% ntile(5) %>% as.factor,
         ses_binary = case_when(ses_sss_composite %in% c(1, 2) ~ 0,
                                ses_sss_composite %in% c(3, 4, 5) ~ 1) %>% as.factor()) %>% 
  dplyr::select(-ses_sss_composite)

keep = X %>% complete.cases()
X = X[keep, ]
y = y[keep, ]
# data = tibble(group = X$group, y %>% dplyr::select(-1)) 

gene = y %>% dplyr::select(-1) %>% as.list()
ses4model = model.matrix(~ ., data = X)
out <-  gene %>%
map(~ lm(.x ~ ses4model))

extract_lm = function(x) x %>% tidy %>% filter(str_detect(term, "ses4modelses_binary1"))

a = out %>% map(extract_lm)
sink("./user_wx/lm_21gene_binary.txt")
a %>% print 
sink()
