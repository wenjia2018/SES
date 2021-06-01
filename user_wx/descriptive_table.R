
#+ echo=F, eval=T, warning=FALSE, message=FALSE
############################################################
library(here)
library(tidyverse)
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(Biobase)


walk(dir(path = here("R"), full.names = TRUE), source)
load_data(reconciled = FALSE, remove_inflam = FALSE)
dt_color_snp <- readRDS("/home/share/dna_ancestry/dna/dt_color_snp.rds")
ancestryPC <- get_PC_dim("aging_mRNA", 0.05)
define_treatments_and_controls_bespoke("aging_mRNA", ancestryPC)
ancestryPC_ses = c(ancestryPC_keep, ses)
custom_PCA <- readRDS(str_c("/home/share/dna_ancestry/dna/custom_PCA_", "aging_mRNA", "_", 0.05, ".rds"))
custom_PCA <- custom_PCA %>%
  select(-fid) %>%
  mutate(AID = AID %>% as.character())
recode_variables_in_dat_bespoke(custom_PCA)

des = readRDS("/home/share/preprocessing/preprocessed_two_batches/waves_22.03.2021.rds")
des_genes =  dat
des_genes1 <- phenoData(des_genes)@data %>% as_tibble()


des_genes2 = des_genes1 %>% dplyr::select(AID, all_of(ancestryPC_ses))

keep = des_genes2 %>% complete.cases()

des_genes2 = des_genes2 %>% mutate(geni = keep) %>% select(AID, geni)

# if whole blood sample
# des_genes2 = des_genes1 %>% dplyr::mutate(geni=T) %>% dplyr::select(AID, geni)

dati = des %>%  left_join(des_genes2, by = "AID")
dati2 = dati %>%
  dplyr::filter(wave5==T) %>%
  dplyr::mutate(genedata = case_when(geni==TRUE ~ 1, is.na(geni)~ 0)) %>% 
  mutate(raceethnicity = re %>%
           fct_recode(
             NonHwhite = "1",
             # white nonhispanic
             NonHblack = "2",
             # black nonhispanic
             NULL = "3",
             # asian nonhispanic
             NULL = "4",
             # other nonhispanic
             Hispanic = "5"
             # hispanic
           ) %>%
           relevel(ref = "NonHwhite"),
         color_byinterviewer3 = H3IR17 %>%
           as.character() %>% 
           as.factor %>% 
           fct_collapse(
             DarkBlack = c("1", "2"),
             LightMed = c("3", "4"),
             White = "5") %>%
           relevel(ref = "White") %>% 
           factor(levels = c("White","LightMed","DarkBlack")))


library(table1)
table1::label(dati2$ses_sss_composite) <- "SES Composite"
table1::label(dati2$raceethnicity) <- "Race and Ethnicity"
table1::label(dati2$color_byinterviewer3) <- "Skin Color"
table1::label(dati2$edu_max) <- "Education"
table1::label(dati2$sss_5) <- "Subjective Social Status"
table1::label(dati2$sex_interv) <- "Sex"
table1::label(dati2$re) <- "Race and Ethnicity"
table1::label(dati2$age_w5) <- "Age at Wave V"
table1::label(dati2$BirthY) <- "Year of birth"
table1::label(dati2$W5REGION) <- "Region"
table1::label(dati2$kit_biow5) <- "Kit"
table1::label(dati2$tube_biow5) <- "Tube"
table1::label(dati2$FastHrs) <- "Fasting Hours"

table1::label(dati2$travel_biow5) <- "Travelling"
table1::label(dati2$months_biow5) <- "Blood Sample Month"
table1::label(dati2$time_biow5) <- "Blood Sample Time"
table1::label(dati2$stress_perceived) <- "Stress Perceived"
table1::label(dati2$bills) <- "Problems Paying Bills"
table1::label(dati2$currentsmoke) <- "Current Smoking"
table1::label(dati2$w5bmi) <- "BMI Wave V"
table1::label(dati2$insurance_lack) <- "Lack of Health Insurance"
table1::label(dati2$lowbirthweight) <- "Low Birthweight"

dati2 <-  dati2 %>% mutate(status= case_when(sex_interv=="m" & genedata==1 ~ "Male mRNA Sample",
                                             sex_interv=="f" & genedata==1 ~ "Female mRNA Sample",
                                             sex_interv=="m" & genedata==0 ~ "Male Wave V Sample",
                                             sex_interv=="f" & genedata==0 ~ "Female Wave V Sample"),
                           genedata = ifelse(genedata==1,"mRNA subsample", NA))

# table1::table1(~ses_sss_composite + color_byinterviewer3 + raceethnicity + sex_interv + age_w5 + BirthY +
#                  W5REGION + kit_biow5 + tube_biow5 + FastHrs + travel_biow5 +months_biow5 + time_biow5 | status, data = dati2)


table1::table1(~ses_sss_composite + color_byinterviewer3 + raceethnicity + sex_interv + age_w5 + BirthY +
                 W5REGION + kit_biow5 + tube_biow5 + FastHrs + travel_biow5 + months_biow5 + time_biow5 | genedata,
               overall = "Wave V",
               caption = "Supp. Table 1. Descriptive Information about the Wave V Total and mRNA Samples, Add Health",
               data = dati2)

# tab = arsenal::tableby(~ses_sss_composite + color_byinterviewer3 + raceethnicity + sex_interv + age_w5 + BirthY +
#                  W5REGION + kit_biow5 + tube_biow5 + FastHrs + travel_biow5 + months_biow5 + time_biow5,
#                  subset = dati2$geni,
#                # overall = "Wave V",
#                # caption = "Supp. Table 1. Descriptive Information about the Wave V Total and mRNA Samples, Add Health",
#                data = dati2)
# summary(tab)

dt = dati2 %>% filter(geni==TRUE)

sjPlot::tab_xtab(dt$raceethnicity, dt$color_byinterviewer3,
                 show.cell.prc = FALSE, show.summary = FALSE,
                 title =  "Supp. Table 2. Self-Described Ethnoracial Designation and Interviewer-Rated Skin Color (Trichotomized), Add Health",
                 var.labels = c("Race Ethnicity", "Skin Color"),
                 emph.total = TRUE,
                 show.row.prc = TRUE, 
                 show.col.prc = TRUE)
