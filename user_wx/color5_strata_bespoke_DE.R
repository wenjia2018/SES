

#+ echo=F, eval=T, warning=FALSE, message=FALSE
NonHblack <- readRDS("~/ses-1/user_wx/color5_NonHblack_strata_bespoke_DE_16.02.2021.rds")
Hispanic <- readRDS("~/ses-1/user_wx/color5_Hispanic_strata_bespoke_DE_16.02.2021.rds")
NonHwhite <- readRDS("~/ses-1/user_wx/color5_NonHwhite_strata_bespoke_DE_16.02.2021.rds")
p_eqtl <- c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
control = "ancestryPC_ses"

library(tidyverse)
# functions to extract data
source("/home/xu/ses-1/user_wx/extract_v2.R")

#' ### sig genes for skincolor within nonhispanic black stratum
#+ echo=F, eval=T, warning=FALSE, message=FALSE
outblack = p_eqtl %>% map(outttT, control, NonHblack)
temp_black = outblack %>%
  bind_rows() %>%
  filter(gene_sig %>% map_dfc(~ length(.))>0) %>% 
  mutate(p_eqtl = p_eqtl %>% format(scientific =T)) %>%
  select(-ttT, -controls)
temp_black %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### sig genes for skincolor within nonhispanic white stratum
#+ echo=F, eval=T, warning=FALSE, message=FALSE
outwhite = p_eqtl %>% map(outttT, control, NonHwhite)
temp_white = outwhite %>%
  bind_rows() %>%
  filter(gene_sig %>% map_dfc(~ length(.))>0) %>% 
  mutate(p_eqtl = p_eqtl %>% format(scientific =T)) %>%
  select(-ttT, -controls)

temp_white %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### sig genes for skincolor within hispanic stratum
#+ echo=F, eval=T, warning=FALSE, message=FALSE
outhispanic = p_eqtl %>% map(outttT, control, Hispanic)
temp_hispanic = outhispanic %>%
  bind_rows() %>%
  filter(gene_sig %>% map_dfc(~ length(.))>0) %>% 
  mutate(p_eqtl = p_eqtl %>% format(scientific =T)) %>%
  select(-ttT, -controls)

temp_hispanic %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()



