#+ echo=F, eval=T
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
#+ echo=F, eval=T
a = readRDS("~/ses-1/user_wx/ses_gsea_webgestaltR_removefig1A.rds")

b = a %>% map(~ .x %>% filter(FDR<0.05))

#' ## ses4

b$ses4 %>%
  select(1,2,5,6,7,8,10) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()
#' 
#' ## income

b$income %>% kableExtra::kable() %>% kableExtra::kable_styling()
#' ## edu

b$edu %>% kableExtra::kable() %>% kableExtra::kable_styling()
#' ## SEI
#' 
b$SEI %>% kableExtra::kable() %>% kableExtra::kable_styling()
#' ## sss
#' 
b$sss %>% kableExtra::kable() %>% kableExtra::kable_styling()