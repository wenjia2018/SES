bespoke1 <- readRDS("~/ses-1/user_wx/color_conti_binary_bespoke.rds")
bespoke2 <- readRDS("~/ses-1/user_wx/color_dummy_bespoke.rds")
bespoke1$out[[1]]$result$example0
bespoke2$out[[1]]$result$example0

for(i in 1:dim(bespoke1)[1]){
bespoke1$out[[i]]$result$example1 =  bespoke1$out[[i]]$result$example1 %>% rbind(bespoke2$out[[i]]$result$example1)
bespoke1$out[[i]]$result$example0 =  bespoke1$out[[i]]$result$example0 %>% rbind(bespoke2$out[[i]]$result$example0)
}
bespoke1 %>% saveRDS("./user_wx/color_conti_binary_dummy_bespoke.rds")

bespoke1 <- readRDS("~/ses-1/user_wx/color_binarycontinuous_bespoke_1902.rds")
bespoke2 <- readRDS("~/ses-1/user_wx/color3_bespoke_18.02.2021.rds")
bespoke3 <- readRDS("~/ses-1/user_wx/color5_bespoke_2202.rds")


for(i in 1:dim(bespoke1)[1]){
  bespoke1$out[[i]]$result$example1 =  bespoke1$out[[i]]$result$example1 %>% 
    rbind(bespoke2$out[[i]]$result$example1) %>% 
    rbind(bespoke3$out[[i]]$result$example1)
  bespoke1$out[[i]]$result$example0 =  bespoke1$out[[i]]$result$example0 %>%
    rbind(bespoke2$out[[i]]$result$example0) %>% 
    rbind(bespoke3$out[[i]]$result$example0)
}
bespoke1 %>% saveRDS("./user_wx/color_conti_binary_dummy3and5_bespoke.rds")

# add aging cluster complement in race
bespoke1 = readRDS("~/ses-1/user_wx/race_bespoke_12.02.2021.rds")
bespoke2 <- readRDS("~/ses-1/user_wx/race_bespoke_aging_cluster_complement_15.03.2021.rds")
bespoke1 %>% rbind(bespoke2) %>% saveRDS("./user_wx/race_bespoke_15.03.2021.rds")

# add aging cluster complement in skincolor3
bespoke1 = readRDS("~/ses-1/user_wx/color3_bespoke_18.02.2021.rds")
bespoke2 <- readRDS("~/ses-1/user_wx/color3_bespoke_aging_cluster_complement_15.03.2021.rds")

bespoke1 %>% rbind(bespoke2) %>% saveRDS("./user_wx/color3_bespoke_15.03.2021.rds")

# add aging cluster complement in skincolor3 in race nonhispanic black strata
bespoke1 = readRDS("/home/xu/ses-1/user_wx/color3_bespoke_NonHblack_strata_25.02.2021.rds")
bespoke2 <- readRDS("~/ses-1/user_wx/color3_bespoke_NonHblack_strata_15.03.2021.rds")

bespoke1 %>% rbind(bespoke2) %>% saveRDS("~/ses-1/user_wx/color3_bespoke_NonHblack_strata_16.03.2021.rds")



color3_bespoke_28.03.2021 <- readRDS("~/ses-1/user_wx/color3_bespoke_28.03.2021.rds")
color3_bespoke_29.03.2021 <- readRDS("~/ses-1/user_wx/color3_bespoke_29.03.2021.rds")

for(i in 1:dim(color3_bespoke_28.03.2021)[1]){
  color3_bespoke_28.03.2021$out[[i]]$result$example0 =  color3_bespoke_28.03.2021$out[[i]]$result$example1 %>% rbind(bespoke2$out[[i]]$result$example1)
  
}
