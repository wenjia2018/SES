
#' ---
#' title: fixed effect model detail information
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

# https://stackoverflow.com/questions/21667262/how-to-find-difference-between-values-in-two-rows-in-an-r-dataframe-using-dplyr
# https://stackoverflow.com/questions/29171403/comparing-rows-within-groups-in-a-dataframe
#' ## Controls, treatment and threshold
#' 
#' ### controls:
#'   basic = 
#' c(
#'  "sex_interv", "famid_fullsib",
#'   "Plate", "AvgCorrelogram100" ,"age_w5",
#'  "NK.cells.activated",
#'  "T.cells.CD8",
#'  "Macrophages.M0", 
#'  "Macrophages.M2",
#'  "B.cells.naive",
#'  "T.cells.CD4.memory.resting"
#' )
#' 

#' ### treatments:
#' 
#' * skin color 3 levels: white as reference, LightMed and DarkBlack as treatments
#' * fam_67 has Nonhispanic white and Hispanic race ethnicity diff
#' 
#' ## sample composition
#+ echo=F, eval=T, warning=FALSE, message=FALSE
tempfamFE <- readRDS("~/ses-1/tempfamFE.rds")

all = tempfamFE$famid_fullsib %>% table %>% as.data.frame()

totalsib = all %>% filter(Freq!=0) %>% dim %>% `[`(1) 

tempfamFE = 
  tempfamFE %>% 
  select(famid_fullsib, color_byinterviewer3, raceethnicity) %>%
  mutate(color_byinterviewer3 = as.character(color_byinterviewer3)) %>% 
  filter(!is.na(famid_fullsib)) %>% 
  group_by(famid_fullsib) %>% 
  mutate(colordiff = ifelse(color_byinterviewer3==lag(color_byinterviewer3), "same", "diff"))

totalsame = tempfamFE$colordiff %>% table %>% as.data.frame() %>% filter(. =="same") %>% pull(Freq)

totaldiff = tempfamFE$colordiff %>% table %>% as.data.frame() %>% filter(. =="diff") %>% pull(Freq)

black = tempfamFE %>% filter(raceethnicity=="NonHblack")

blacksib = black$famid_fullsib %>% table %>% as.data.frame() %>% filter(Freq!=0) %>% dim %>% `[`(1) 

blacksame = black$colordiff %>% table %>% as.data.frame() %>% filter(. =="same") %>% pull(Freq)
blackdiff = black$colordiff %>% table %>% as.data.frame() %>% filter(. =="diff") %>% pull(Freq)

white = tempfamFE %>% filter(raceethnicity=="NonHwhite")

whitesib = white$famid_fullsib %>% table %>% as.data.frame() %>% filter(Freq!=0) %>% dim %>% `[`(1) 

whitesame = white$colordiff %>% table %>% as.data.frame() %>% filter(. =="same") %>% pull(Freq)
whitediff = white$colordiff %>% table %>% as.data.frame() %>% filter(. =="diff") %>% pull(Freq)

his = tempfamFE %>% filter(raceethnicity=="Hispanic")

hissib = his$famid_fullsib %>% table %>% as.data.frame() %>% filter(Freq!=0) %>% dim %>% `[`(1) 

hissame = his$colordiff %>% table %>% as.data.frame() %>% filter(. =="same") %>% pull(Freq)
hisdiff = his$colordiff %>% table %>% as.data.frame() %>% filter(. =="diff") %>% pull(Freq)
tribble(
  ~ sample,  ~ NO.sib, ~ "color rated same sib", ~ "color rated different sib", ~ "percentage of rated different",
  "total", totalsib, totalsame, totaldiff, str_c(totaldiff/totalsib*100, "%"),
  "White", whitesib, whitesame, whitediff, str_c(whitediff/whitesib*100, "%"),
  "Black", blacksib, blacksame, blackdiff, str_c(blackdiff/blacksib*100, "%"),
  "Hispanic", hissib, hissame, hisdiff, str_c(hisdiff/hissib*100, "%")
  
) %>% kableExtra::kable() %>% 
  kableExtra::kable_styling()

#' ## intra-family correlation Cohen’s kappa 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
# https://www.r-bloggers.com/2018/04/k-is-for-cohens-kappa/
comp =
  tempfamFE %>% 
  select(famid_fullsib, color_byinterviewer3) %>% 
  pivot_wider(names_from = "famid_fullsib", values_from = "color_byinterviewer3") %>% 
  t %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "famid_fullsib") %>% 
  unnest_wider("V1")


kappa2(comp %>% select(2,3))
agree(comp %>% select(2,3))

comp_white = comp %>% 
  right_join(tempfamFE %>% select(famid_fullsib, raceethnicity) %>% unique %>% filter(raceethnicity =="NonHwhite"))

kappa2(comp_white %>% select(2,3))

comp_black = comp %>% 
  right_join(tempfamFE %>% select(famid_fullsib, raceethnicity) %>% unique %>% filter(raceethnicity =="NonHblack"))

kappa2(comp_black %>% select(2,3))

comp_his = comp %>% 
  right_join(tempfamFE %>% select(famid_fullsib, raceethnicity) %>% unique %>% filter(raceethnicity =="Hispanic"))

kappa2(comp_his %>% select(2,3))


#' ### treatments:
#' 
#' * skin color 5 levels: white, Light brown, Medium brown and Dark brown, and Black 
#' 
#' ## sample composition
#+ echo=F, eval=T, warning=FALSE, message=FALSE
full <- readRDS("~/ses-1/user_wx/tempfamFE.rds")
tempfamFE = full

all = tempfamFE$famid_fullsib %>% table %>% as.data.frame()

totalsib = all %>% filter(Freq!=0) %>% dim %>% `[`(1) 

tempfamFE = 
  tempfamFE %>% 
  select(famid_fullsib, color_byinterviewer5, raceethnicity) %>%
  mutate(color_byinterviewer5 = as.character(color_byinterviewer5)) %>% 
  filter(!is.na(famid_fullsib)) %>% 
  group_by(famid_fullsib) %>% 
  mutate(colordiff = ifelse(color_byinterviewer5==lag(color_byinterviewer5), "same", "diff"))

totalsame = tempfamFE$colordiff %>% table %>% as.data.frame() %>% filter(. =="same") %>% pull(Freq)

totaldiff = tempfamFE$colordiff %>% table %>% as.data.frame() %>% filter(. =="diff") %>% pull(Freq)

black = tempfamFE %>% filter(raceethnicity=="NonHblack")

blacksib = black$famid_fullsib %>% table %>% as.data.frame() %>% filter(Freq!=0) %>% dim %>% `[`(1) 

blacksame = black$colordiff %>% table %>% as.data.frame() %>% filter(. =="same") %>% pull(Freq)
blackdiff = black$colordiff %>% table %>% as.data.frame() %>% filter(. =="diff") %>% pull(Freq)

white = tempfamFE %>% filter(raceethnicity=="NonHwhite")

whitesib = white$famid_fullsib %>% table %>% as.data.frame() %>% filter(Freq!=0) %>% dim %>% `[`(1) 

whitesame = white$colordiff %>% table %>% as.data.frame() %>% filter(. =="same") %>% pull(Freq)
whitediff = white$colordiff %>% table %>% as.data.frame() %>% filter(. =="diff") %>% pull(Freq)

his = tempfamFE %>% filter(raceethnicity=="Hispanic")

hissib = his$famid_fullsib %>% table %>% as.data.frame() %>% filter(Freq!=0) %>% dim %>% `[`(1) 

hissame = his$colordiff %>% table %>% as.data.frame() %>% filter(. =="same") %>% pull(Freq)
hisdiff = his$colordiff %>% table %>% as.data.frame() %>% filter(. =="diff") %>% pull(Freq)
tribble(
  ~ sample,  ~ NO.sib, ~ "color rated same sib", ~ "color rated different sib", ~ "percentage of rated different", ~ "Cohen's K"
  "total", totalsib, totalsame, totaldiff, str_c(totaldiff/totalsib*100, "%"), 
  "White", whitesib, whitesame, whitediff, str_c(whitediff/whitesib*100, "%"),
  "Black", blacksib, blacksame, blackdiff, str_c(blackdiff/blacksib*100, "%"),
  "Hispanic", hissib, hissame, hisdiff, str_c(hisdiff/hissib*100, "%")
  
) %>% kableExtra::kable() %>% 
  kableExtra::kable_styling()

#' ## intra-family correlation Cohen’s kappa 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
# https://www.r-bloggers.com/2018/04/k-is-for-cohens-kappa/
comp =
  tempfamFE %>% 
  select(famid_fullsib, color_byinterviewer5) %>% 
  pivot_wider(names_from = "famid_fullsib", values_from = "color_byinterviewer5") %>% 
  t %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "famid_fullsib") %>% 
  unnest_wider("V1")


all_k = kappa2(comp %>% select(2,3))
agree(comp %>% select(2,3))

comp_white = comp %>% 
  right_join(tempfamFE %>% select(famid_fullsib, raceethnicity) %>% unique %>% filter(raceethnicity =="NonHwhite"))

white_k = kappa2(comp_white %>% select(2,3))

comp_black = comp %>% 
  right_join(tempfamFE %>% select(famid_fullsib, raceethnicity) %>% unique %>% filter(raceethnicity =="NonHblack"))

black_k = kappa2(comp_black %>% select(2,3))

comp_his = comp %>% 
  right_join(tempfamFE %>% select(famid_fullsib, raceethnicity) %>% unique %>% filter(raceethnicity =="Hispanic"))

his_k = kappa2(comp_his %>% select(2,3))


tribble(
  ~ sample,  ~ NO.sib, ~ "color rated same sib", ~ "color rated different sib", ~ "percentage of rated different", ~ "Cohen's K(Intraclass correlation)",
  "total", totalsib, totalsame, totaldiff, str_c(totaldiff/totalsib*100, "%"), all_k$value,
  "White", whitesib, whitesame, whitediff, str_c(whitediff/whitesib*100, "%"), white_k$value,
  "Black", blacksib, blacksame, blackdiff, str_c(blackdiff/blacksib*100, "%"), black_k$value,
  "Hispanic", hissib, hissame, hisdiff, str_c(hisdiff/hissib*100, "%"), his_k$value) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()
