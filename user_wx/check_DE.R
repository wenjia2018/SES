# check DE genes from geonowide regression controlling for the ancestry PC for the whole genomoe with
# DE genes from genowide regression controlling for the ancestry PC for each signature

load_data(reconciled = FALSE, remove_inflam = FALSE)
p_eqtl = 0.01
# m8
example_race = readRDS("~/ses-1/user_wx/race_bespoke_15.03.2021.rds")
example_skincolor3 = readRDS("~/ses-1/user_wx/color3_bespoke_15.03.2021.rds")

race_m8 = f0(p = p_eqtl, example_race) %>% 
  filter(control_set == "ancestryPC_ses") %>%
  hoist(out, ttT = list("result", "m8_fdr", 1, "other", "m")) %>% 
  unnest(ttT) %>% 
  select(-out, - controls, table1) %>% 
  filter(adj.P.Val <0.05)

skincolor_m8 = f0(p = p_eqtl, example_skincolor3) %>% 
  filter(control_set == "ancestryPC_ses") %>%
  hoist(out, ttT = list("result", "m8_fdr", 1, "other", "m")) %>% 
  unnest(ttT) %>% 
  select(-out, - controls, table1) %>% 
  filter(adj.P.Val <0.05)

race_skincolor_m8 = race_m8 %>% rbind(skincolor_m8)

# genowide DE
skincolor3_bespoke_DE_22.03.2021 <- readRDS("~/ses-1/user_wx/skincolor3_bespoke_DE_22.03.2021.rds")

race_bespoke_DE_22.03.2021 <- readRDS("~/ses-1/user_wx/race_bespoke_DE_22.03.2021.rds")

race_DE = outttT(p = p_eqtl, control = "ancestryPC_ses", race_bespoke_DE_22.03.2021)
skincolor_DE = outttT(p = p_eqtl, control = "ancestryPC_ses", skincolor3_bespoke_DE_22.03.2021)

race_skincolor_DE = skincolor_DE %>%
  rbind(race_DE) %>% 
  unnest(ttT) %>% 
  filter(adj.P.Val <0.05) %>% 
  mutate()


your_func <- function(df, x){
  
  df %>% 
    # group_by(z) %>%
    transmute(!!paste(x) := ifelse(gene %in% signatures$outcome_set[[x]], TRUE, FALSE ))
    # ungroup() %>%
    # select(!!paste(x, "ratio", sep = "_") )
}

a=race_skincolor_DE %>% 
  bind_cols(map_dfc(focal, ~ your_func(race_skincolor_DE, .x)))


focal = c("aging_mRNA",
    "aging_up_mRNA",
    "aging_down_mRNA",
    "aging_down_cl1_mRNA",
    "aging_down_cl1a_mRNA",
    "aging_down_cl1b_mRNA",
    "aging_down_cl1c_mRNA",
    "aging_down_cl2_mRNA",
    "aging_down_cl3_mRNA",
    "aging_up_cl1_mRNA",
    "aging_up_cl2_mRNA",
    "aging_up_cl3_mRNA",
    "aging_up_cl4_mRNA",
    "aging_cluster_complement")


