#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
library(here)
library(data.table)
# setDTthreads(threads = 20)
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


walk(dir(path = here("R"),full.names = TRUE), source)


############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################
# choose normalization methods for downstream analysis
tmm = TRUE
rle = FALSE
log2cpm = FALSE
# which PCA to perform
oblimin = FALSE
nn = TRUE
# explicitly assign ncomp as the smallest number of table signatures gene numbers
ncomp = 10
# for doing genowide DE analysis only
normalization_bydesign = FALSE
# specify if subjects with disease shall be removed
remove_diseased_subjects = TRUE
load_data(reconciled = FALSE, remove_inflam = FALSE)
define_treatments_and_controls()
recode_variables_in_dat()

mediation_mean = FALSE
mediation_each_gene = FALSE

dt = pData(dat)
dt = dt %>% 
  dplyr::mutate(H5ID30 = case_when(H5ID30 ==0 ~ "NO",
                            H5ID30 ==1 ~ "YES"),
                H5TO12 = case_when(H5TO12==0 ~ "none",
                                   H5TO12==1 ~ "1 or 2 days in the past 12 months",
                                   H5TO12==2 ~ "once a month or less (3 to 12 days in the past 12 months)",
                                   H5TO12==3 ~ "2 or 3 days a month",
                                   H5TO12==4 ~ "1 or 2 days a week",
                                   H5TO12==5 ~ "3 to 5 days a week",
                                   H5TO12==6 ~ "every day or almost every day") %>% 
                  factor(levels = c("none","1 or 2 days in the past 12 months","once a month or less (3 to 12 days in the past 12 months)",
                                       "2 or 3 days a month","1 or 2 days a week","3 to 5 days a week","every day or almost every day")))
#' ### H5ID30 Were the past 7 days typical in terms of your physical activity?
#+ echo=F, eval=T, warning=FALSE, message=FALSE
b = dt$H5ID30 %>% table() %>% as.data.frame()
colnames(b) = c("physical activity typical in the past 7 days","Freq")
b %>% kableExtra::kable() %>% kableExtra::kable_styling()
#' ### H5TO12 During the past 12 months, on how many days did you drink alcohol (beer, wine, or liquor)?
#+ echo=F, eval=T, warning=FALSE, message=FALSE

a = dt$H5TO12 %>% table %>% as.data.frame()
colnames(a) = c("drink in the past 12 months","Freq")
a %>% kableExtra::kable() %>% kableExtra::kable_styling()
#' ###  wealth (H5EC2 + H5EC4) - (H5EC5A + H5EC5B + H5EC5C) = wealth
#' * if H5EC2 is NA use H5EC1 as single people only answer H5EC1, and H5EC2 is legitimate skip
#' * If H5EC3 is 0(means no gifts and inheritances), set H5EC4 to 0
#+ echo=F, eval=T, warning=FALSE, message=FALSE
dt$assets_household_net_w5 %>% summary
dt$assets_household_net_w5 %>% hist(breaks = 80)
#' ###  wealth (H5EC2 + H5EC4) - (H5EC5A + H5EC5B + H5EC5C) = wealth
#' * if H5EC2 is NA use H5EC1 as single people only answer H5EC1, and H5EC2 is legitimate skip
#' * If H5EC3 is 0(means no gifts and inheritances), set H5EC4 to 0
#' * remove all the missing records, does not change too much
#+ echo=F, eval=T, warning=FALSE, message=FALSE
temp_wealth = dt %>% select(AID, H5EC1, H5EC2, H5EC4, H5EC5A, H5EC5B, H5EC5C) %>% 
  mutate(H5EC2 = ifelse(is.na(H5EC2), H5EC1, H5EC2))
keep = complete.cases(temp_wealth)
keep %>% table
temp = temp_wealth[keep,]
temp = temp %>% mutate(wealth = (H5EC2 + H5EC4) - (H5EC5A + H5EC5B + H5EC5C))
temp$wealth %>% summary
temp$wealth %>% hist(breaks = 100)
table(temp$wealth<0)