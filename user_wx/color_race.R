
#' ---
#' title: race and skin color (by interviewer)
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

#' color_byinterviewer = H3IR17
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
library(tidyverse)
library(foreign)
library(descr)
# dat <- readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_17.11.2020.rds")
dat <- readRDS("/home/share/preprocessing/preprocessed_two_batches/all.batches.expression.set.tmm_waves_01.09.2021.rds")
AID_blood = dat@phenoData@data$AID
waves <- readRDS("/home/share/preprocessing/preprocessed_two_batches/waves_01.09.2021.rds")
# waves <- readRDS("/home/share/preprocessed_two_batches/waves_17.11.2020.rds")
waves_full <- readRDS("/home/share/preprocessing/preprocessed_two_batches/waves_full.rds")

sibling_w3 = read.xport("/home/share/data_input/sibling3.xpt") %>% as_tibble
ancestralPCA <- read.delim("/home/share/data_input/ancestralPCA.eigenStratPCs.allparticipants.082719")
ancestralPCA = ancestralPCA %>% mutate(AID=aid %>%as.character) %>% select(AID, fid)


sibling_full = sibling_w3 %>% 
  mutate_at(vars(c(SIB_AID1, SIB_AID2, SIB_AID3, SIB_AID4)), list(~ .x %>% as.character())) %>%
  mutate(SIB_AID1 = ifelse(SIB_REL1 %in% c("brother", "sister", "twin brother", "twin sister"), SIB_AID1, NA))%>%
  mutate(SIB_AID2 = ifelse(SIB_REL2 %in% c("brother", "sister", "twin brother", "twin sister"), SIB_AID2, NA))%>%
  mutate(SIB_AID3 = ifelse(SIB_REL3 %in% c("brother", "sister", "twin brother", "twin sister"), SIB_AID3, NA))%>%
  mutate(SIB_AID4 = ifelse(SIB_REL4 %in% c("brother", "sister", "twin brother", "twin sister"), SIB_AID4, NA)) %>% 
  filter(!((is.na(SIB_AID1) &is.na(SIB_AID2)&is.na(SIB_AID3)&is.na(SIB_AID4))))


temp_fid = sibling_full %>% left_join(ancestralPCA) 

sibling_AID = sibling_full$AID %>%
  union(sibling_full$SIB_AID1) %>% 
  union(sibling_full$SIB_AID2) %>% 
  union(sibling_full$SIB_AID3) %>% 
  union(sibling_full$SIB_AID4) %>%
  unique

waves = waves %>% mutate(
  raceethnicity = re %>%
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
    ))%>% 
  filter(!is.na(raceethnicity))
waves_full = waves_full %>% mutate(  
color_byinterviewer = H3IR17 %>%
    as.character() %>% 
    as.factor %>% 
    fct_recode(
    black = "1",
    dark_brown = "2",
    medium_brown = "3",
    light_brown = "4",
    white = "5"
  ))

#' ## blood sample
#+ echo=F, eval=T, warning=FALSE, message=FALSE

df = waves %>% left_join(waves_full %>% select(AID,color_byinterviewer)) %>% filter(AID %in% AID_blood)

crosstab(df$raceethnicity,  df$color_byinterviewer, prop.r = T, prop.c = T, prop.t = FALSE, plot = F)

#' ## full siblings in blood sample
#' * in the sibling dataset, there are brother, half-brother, half-sister, non-related, sister, twin brother,  twin sister 
#' * I keep only brother, sister, twin brother and twin sister
#+ echo=F, eval=T, warning=FALSE, message=FALSE
df_full_sibling = df %>% filter(AID %in% sibling_AID) %>% select(raceethnicity, color_byinterviewer)

crosstab(df_full_sibling$raceethnicity,  df_full_sibling$color_byinterviewer, prop.r = T, prop.c = T, prop.t = F, plot = F)



temp_fid = temp_fid%>% filter(AID %in% AID_blood)

temp_fid_ownfid = temp_fid %>% left_join(waves %>% select(AID,famid_fullsib))

