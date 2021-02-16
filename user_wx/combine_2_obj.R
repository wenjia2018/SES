bespoke1 <- readRDS("~/ses-1/user_wx/color_conti_binary_bespoke.rds")
bespoke2 <- readRDS("~/ses-1/user_wx/color_dummy_bespoke.rds")
bespoke1$out[[1]]$result$example0
bespoke2$out[[1]]$result$example0

for(i in 1:dim(bespoke1)[1]){
bespoke1$out[[i]]$result$example1 =  bespoke1$out[[i]]$result$example1 %>% rbind(bespoke2$out[[i]]$result$example1)
bespoke1$out[[i]]$result$example0 =  bespoke1$out[[i]]$result$example0 %>% rbind(bespoke2$out[[i]]$result$example0)
}
bespoke1 %>% saveRDS("./user_wx/color_conti_binary_dummy_bespoke.rds")
