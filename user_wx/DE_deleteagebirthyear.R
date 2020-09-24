
DE_removeagebirthyear <- readRDS("~/ses-1/user_wx/DE_removeagebirthyear.rds")

error = DE_removeagebirthyear %>%
  hoist(out, "error") %>%
  mutate(error = map(error, as.character)) %>%
  unnest(error)

# RESULTS
temp = DE_removeagebirthyear %>%
  hoist(out, "result")


result = temp %>%
  pluck("result") %>%
  set_names(temp$treatment)


DE_full = result %>% 
  map(pluck("ttT")) %>% 
  map(~ filter(.x,adj.P.Val<=0.05))
#' #### ses DE genes
#' 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
DE_updown = result %>% 
  map(pluck("ttT")) %>% 
  map(~ filter(.x,adj.P.Val<=0.05) %>% mutate(updown = ifelse(logFC>0,"up","down")) %>% dplyr::select(1,8))

DE = result %>% 
  map(pluck("ttT")) %>% 
  map(~ filter(.x,adj.P.Val<=0.05) %>% pull(gene))

sig = Reduce(union, signatures$outcome_set[table1[1:11]])

DE_remove1k = map(DE, ~ setdiff(.x, sig))
venn::venn(DE_remove1k, ilabels = TRUE, zcolor = "style", ellipse = FALSE, opacity = 0.15)
# ses 4 and income unique genes after removing all genes in our table 1 disease signatures
ses4_unique = Reduce(setdiff, list(DE_remove1k$ses_sss_composite, DE_remove1k$income_hh_ff5, DE_remove1k$sss_5, DE_remove1k$edu_max, DE_remove1k$SEI_ff5))
income_unique = Reduce(setdiff, list(DE_remove1k$income_hh_ff5, DE_remove1k$ses_sss_composite, DE_remove1k$sss_5, DE_remove1k$edu_max, DE_remove1k$SEI_ff5))

# ses 4 up and down
ses4_up = DE_updown$ses_sss_composite %>% filter(updown =="up") %>% pull(gene) %>% intersect(DE_remove1k$ses_sss_composite)
ses4_down = DE_updown$ses_sss_composite %>% filter(updown =="down") %>% pull(gene) %>% intersect(DE_remove1k$ses_sss_composite)

# income up and down
income_up = DE_updown$income_hh_ff5 %>% filter(updown =="up") %>% pull(gene) %>% intersect(DE_remove1k$income_hh_ff5)
income_down = DE_updown$income_hh_ff5 %>% filter(updown =="down") %>% pull(gene) %>% intersect(DE_remove1k$income_hh_ff5)

# ses4 unique up and down
ses4_unique_up = DE_updown$ses_sss_composite %>% filter(updown =="up") %>% pull(gene) %>% intersect(ses4_unique)
ses4_unique_down = DE_updown$ses_sss_composite %>% filter(updown =="down") %>% pull(gene) %>% intersect(ses4_unique)


# income unique up and down
income_unique_up = DE_updown$income_hh_ff5 %>% filter(updown =="up") %>% pull(gene) %>% intersect(income_unique)
income_unique_down = DE_updown$income_hh_ff5 %>% filter(updown =="down") %>% pull(gene) %>% intersect(income_unique)


# ses 4 and income intersection up and down

ses4_income_intersection = intersect(DE_remove1k$ses_sss_composite, DE_remove1k$income_hh_ff5)
ses4_income_intersection_up = intersect(ses4_up, income_up)
ses4_income_intersection_down = intersect(ses4_down, income_down)

a = list(ses4 = DE_remove1k$ses_sss_composite,
         ses4_up = ses4_up,
         ses4_down = ses4_down,
         ses4_unique = ses4_unique, ses4_unique_up = ses4_unique_up, ses4_unique_down = ses4_unique_down,
         income = DE_remove1k$income_hh_ff5,
         income_up = income_up,
         income_down = income_down,
         income_unique = income_unique,
         income_unique_up = income_unique_up,
         income_unique_down =income_unique_down,
         ses4_income_intersection = ses4_income_intersection,
         ses4_income_intersection_up = ses4_income_intersection_up,
         ses4_income_intersection_down = ses4_income_intersection_down)


b = list(ses4_unique_down = tibble(gene = ses4_unique_down) %>% left_join(DE_full$ses_sss_composite),
         income_unique_down = tibble(gene = income_unique_down) %>% left_join(DE_full$income_hh_ff5),
         ses4_income_intersection_down = tibble(gene = ses4_income_intersection_down)%>% left_join(DE_full$ses_sss_composite))

b %>%  openxlsx::write.xlsx("./user_wx/DE_removetable1signatures_ses4income_fulllist_removeagebirthyear.xlsx")

