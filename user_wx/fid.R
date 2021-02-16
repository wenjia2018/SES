
#' ---
#' title: family id
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

#' 
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
library(tidyverse)
library(foreign)
library(descr)
dat <- readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_17.11.2020.rds")
AID_blood = dat@phenoData@data$AID
waves <- readRDS("/home/share/preprocessed_two_batches/waves_17.11.2020.rds")
otherrace = waves %>% filter(re %in%c(3,4)) %>% select(AID,re,H3IR17) %>%  pluck("AID")

sibling_w3 = read.xport("/home/share/data_input/sibling3.xpt") %>% as_tibble
ancestralPCA <- read.delim("/home/share/data_input/ancestralPCA.eigenStratPCs.allparticipants.082719")
ancestralPCA = ancestralPCA %>% mutate(AID=aid %>%as.character) %>% select(AID, fid)
# construct full sibling data 
sibling_full = sibling_w3 %>% 
  mutate_at(vars(c(SIB_AID1, SIB_AID2, SIB_AID3, SIB_AID4)), list(~ .x %>% as.character())) %>%
  mutate(SIB_AID1 = ifelse(SIB_REL1 %in% c("brother", "sister", "twin brother", "twin sister"), ifelse(SIB_AID1 != AID, SIB_AID1,NA), NA))%>%
  mutate(SIB_AID2 = ifelse(SIB_REL2 %in% c("brother", "sister", "twin brother", "twin sister"), ifelse(SIB_AID2 != AID, SIB_AID2,NA), NA))%>%
  mutate(SIB_AID3 = ifelse(SIB_REL3 %in% c("brother", "sister", "twin brother", "twin sister"), ifelse(SIB_AID3 != AID, SIB_AID3,NA), NA))%>%
  mutate(SIB_AID4 = ifelse(SIB_REL4 %in% c("brother", "sister", "twin brother", "twin sister"), ifelse(SIB_AID4 != AID, SIB_AID4,NA), NA)) %>% 
  filter(!((is.na(SIB_AID1) & is.na(SIB_AID2) & is.na(SIB_AID3) & is.na(SIB_AID4)))) 


fullsib_AID = sibling_full$AID %>%
  union(sibling_full$SIB_AID1) %>% 
  union(sibling_full$SIB_AID2) %>% 
  union(sibling_full$SIB_AID3) %>% 
  union(sibling_full$SIB_AID4) %>%
  unique %>%
  as_tibble_col(column_name = "AID")


sibling_w3 = sibling_w3 %>% left_join(ancestralPCA)

fid = list(t0 = sibling_w3 %>% select(AID, fid),
           t1 = sibling_w3 %>% select(SIB_AID1, fid) %>% dplyr::rename(AID = SIB_AID1),
           t2 = sibling_w3 %>% select(SIB_AID2, fid) %>% dplyr::rename(AID = SIB_AID2),
           t3 = sibling_w3 %>% select(SIB_AID3, fid) %>% dplyr::rename(AID = SIB_AID3),
           t4 = sibling_w3 %>% select(SIB_AID4, fid) %>% dplyr::rename(AID = SIB_AID4)) %>% 
  purrr::reduce(bind_rows) %>% 
  filter(AID!="") %>% 
  filter(!is.na(fid)) %>% 
  unique()


temp_fid = fullsib_AID %>% left_join(fid) %>% filter(AID %in% AID_blood) 
a=temp_fid$fid %>% table %>% as.data.frame()
a$Freq %>% table

temp_fid_oursample = temp_fid %>% filter(!(AID %in% otherrace))
a=temp_fid_oursample$fid %>% table %>% as.data.frame()
a$Freq %>% table

a

b = waves_fullsibling%>% filter(AID %in% AID_blood) %>% filter(!(AID %in% otherrace))
b=waves_fullsibling$famid_fullsib %>% table %>% as.data.frame()
b$Freq %>% table
