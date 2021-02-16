
#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(here)
library(tidyverse)

walk(dir(path = here("R"), full.names = TRUE), source)
gene_set_name = "aging_mRNA"
p_eqtl = 0.05
load_data(reconciled = FALSE, remove_inflam = FALSE)
ancestryPC <- get_PC_dim(gene_set_name, p_eqtl)
define_treatments_and_controls_bespoke(gene_set_name, ancestryPC)
# define_treatments_and_controls_snps(gene_set_name, ancestryPC)
custom_PCA <- readRDS(str_c("/home/share/dna_ancestry/dna/custom_PCA_", gene_set_name, "_", p_eqtl, ".rds"))
custom_PCA <- custom_PCA %>%
  select(-fid) %>%
  mutate(AID = AID %>% as.character())
recode_variables_in_dat_bespoke(custom_PCA)
var_set = c("color_byinterviewer3_DarkBlack", "color_byinterviewer3_LightMed",
            "raceethnicity_NonHblack", "raceethnicity_Hispanic",
            "raceethnicity__color_byinterviewer3_NonHblack|DarkBlack",
            "raceethnicity__color_byinterviewer3_NonHblack|LightMed",
            "raceethnicity__color_byinterviewer3_Hispanic|DarkBlack",
            "raceethnicity__color_byinterviewer3_Hispanic|LightMed")
dt1 = pData(dat) %>% select(raceethnicity, color_byinterviewer3)
descr::crosstab(dt1$raceethnicity,dt1$color_byinterviewer3,prop.c = T, prop.t = T, plot = F)

dt = pData(dat) %>% mutate(outcome = ses_sss_composite) %>% select(outcome,all_of(var_set))
#+ echo=F, eval=T, warning=FALSE, message=FALSE, results = "asis"
fit <- lm(outcome ~ ., dt)
# show the theoretical model
equatiomatic::extract_eq(fit, wrap = TRUE, terms_per_line = 1)

#' ## effects
#' * effect of group NonHwhite and white : $\alpha$ 
#' * effect of group NonHwhite and LightMed : $\alpha + \beta_{2}$
#' * effect of group NonHwhite and DarkBlack : $\alpha + \beta_{1}$
#' * effect of group NonHblack and white : $\alpha + \beta_{3}$
#' * effect of group NonHblack and LightMed : $\alpha + \beta_{3} + \beta_{2} + \beta_{6}$
#' * effect of group NonHblack and DarkBlack : $\alpha + \beta_{3} + \beta_{1} + \beta_{5}$
#' * effect of group Hispanic and white : $\alpha + \beta_{4}$
#' * effect of group Hispanic and LightMed : $\alpha + \beta_{2} + \beta_{4} + \beta_{8}$
#' * effect of group Hispanic and DarkBlack : $\alpha + \beta_{1} + \beta_{4} + \beta_{7}$
#' 
#' 
#' 
