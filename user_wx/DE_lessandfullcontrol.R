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
# print(abbreviations)
# funcs = str_subset(abbreviations$shorthand, "^m") 
# funcs = funcs %>% str_subset("m[6-8]")
# explicitly assign ncomp as the smallest number of table signatures gene numbers
ncomp = signatures$outcome_set[table1]%>% map_dfc(length) %>% unlist() %>%  min

# full control

example1_fullcontrol =
  args %>%
  filter(treatment %in%c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5",
                         "ses_composite_pp1", "edu_p", "SEI_max_p_w12", "income_pp1_log"),
         gene_set_name == "whole_genome_and_tfbm",
         names(controls) == "all") %>%
  mutate(out = pmap(., safely(model_fit), funcs),
         controls = names(controls))

# # less control
# 
# controls$basic = c("sex_interv", "re", "Plate", "AvgCorrelogram100")
# 
#   example1_lesscontrol =
#     args %>%
#     filter(treatment %in%c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5"),
#            gene_set_name == "whole_genome_and_tfbm",
#            names(controls) == "basic") %>%
#     mutate(out = pmap(., safely(model_fit), funcs),
#            controls = names(controls))

  
  
# extract information 
  
  example1 %>%
    hoist(out, "error") %>%
    mutate(error = map(error, as.character)) %>%
    unnest(error)
  
  # RESULTS
  temp = example1 %>%
    hoist(out, "result")
  
  
  result = temp %>%
    pluck("result") %>%
    set_names(temp$treatment)
  
  #' #### ses DE genes
  #' 
  #+ echo=F, eval=T, warning=FALSE, message=FALSE
  DE = result %>% 
    map(pluck("ttT")) %>% 
    map(~ filter(.x,adj.P.Val<=0.05) %>% pull(gene))
  
  DE %>% openxlsx::write.xlsx("./user_wx/DE_xxx.xlsx")
  #' #### venn diagram for ses DE genes
  #' 
  #+ echo=F, eval=T, warning=FALSE, message=FALSE
  
  venn::venn(DE, ilabels = TRUE, zcolor = "style", ellipse = FALSE, opacity = 0.15)
  
