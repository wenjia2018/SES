# treatment = treatment[5]
# 
# controls = controls$ancestryPC_ses
ftest_v = c("raceethnicity__color_byinterviewer3_NonHwhite|DarkBlack", 
            "raceethnicity__color_byinterviewer3_NonHwhite|LightMed")
treatment = a[3]
# ftest_v = str_subset(c(treatment) %>% unique, "__")

if(treatment %in% ftest_v) {
  ftest_v = replace(ftest_v, ftest_v == treatment, "treatment")
}
ftest_v




a=m$fit$d1 %>% car::linearHypothesis(str_subset(names(coef(m$fit$d1)), ftest_v), verbose=F)

example