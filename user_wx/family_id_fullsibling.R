# create family id for full sibling: "brother", "sister", "twin brother", "twin sister"
library(tidyverse)
library(foreign)

sibling_w3 = read.xport("/home/share/data_input/sibling3.xpt") %>% as_tibble
dat <- readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_17.11.2020.rds")
AID_blood = dat@phenoData@data$AID


# keep only the AID who has full siblings, convert the rest to NA
sibling_full = sibling_w3 %>% 
  mutate_at(vars(c(SIB_AID1, SIB_AID2, SIB_AID3, SIB_AID4)), list(~ .x %>% as.character())) %>%
  mutate(SIB_AID1 = ifelse(SIB_REL1 %in% c("brother", "sister", "twin brother", "twin sister"), SIB_AID1, NA))%>%
  mutate(SIB_AID2 = ifelse(SIB_REL2 %in% c("brother", "sister", "twin brother", "twin sister"), SIB_AID2, NA))%>%
  mutate(SIB_AID3 = ifelse(SIB_REL3 %in% c("brother", "sister", "twin brother", "twin sister"), SIB_AID3, NA))%>%
  mutate(SIB_AID4 = ifelse(SIB_REL4 %in% c("brother", "sister", "twin brother", "twin sister"), SIB_AID4, NA)) %>% 
  filter(!((is.na(SIB_AID1) &is.na(SIB_AID2)&is.na(SIB_AID3)&is.na(SIB_AID4)))) %>% 
  mutate(famid = str_c("fam_",c(1:dim(.)[1])))

# reconstruct a dataframe for all the AID who has full siblings as one column, and their famid as another column
# there are some duplicated rows
t0 = sibling_full %>% select(AID, famid)
t1 = sibling_full %>% select(SIB_AID1, famid) %>% rename(AID = SIB_AID1)
t2 = sibling_full %>% select(SIB_AID2, famid) %>% rename(AID = SIB_AID2)
t3 = sibling_full %>% select(SIB_AID3, famid) %>% rename(AID = SIB_AID3)
t4 = sibling_full %>% select(SIB_AID4, famid) %>% rename(AID = SIB_AID4)

a = purrr::reduce(list(t0, t1, t2, t3, t4), bind_rows) %>% filter(!is.na(AID))

# remove the duplicated rows where the AIDs are same 
b = a %>% nest(AIDs = AID) %>% mutate(AIDs = AIDs %>% unlist(recursive = F))  %>% 
  mutate(AIDs =AIDs %>% map(str_replace_na) %>% map(as.numeric) %>% map(sort))%>% column_to_rownames(var = "famid") %>% unique() 

# unnest AIDs for next step 
cc = b %>% tidyr::separate(col = AIDs, sep = ",", into = c("AID1", "AID2", "AID3")) %>% rownames_to_column("famid") %>% 
  mutate_at(vars(c("AID1", "AID2", "AID3")), ~ gsub("[^0-9.-]", "",.))



t1 = cc %>% select(AID1, famid) %>% rename(AID = AID1)
t2 = cc %>% select(AID2, famid) %>% rename(AID = AID2)
t3 = cc %>% select(AID3, famid) %>% rename(AID = AID3)
# full siblings AIDs and their family id
aa = purrr::reduce(list(t1, t2, t3), bind_rows) 



