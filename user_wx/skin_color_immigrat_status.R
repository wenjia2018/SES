
#' ---
#' title: skin color and immigration status summary
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

#' #### Skin color and Immigrate status 
#' * immigrate status is created by variable H1GI11(Were you born in the United States?)
#' * immigrate status = no if H1GI11 = 1, immigrate status = yes if H1GI11 = 0

#+ echo=F, eval=T, warning=FALSE, message=FALSE 
library(tidyverse)
# dat <- readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_22.03.2021.rds")
dat <- readRDS("/home/share/preprocessing/preprocessed_two_batches/all.batches.expression.set.tmm_waves_01.09.2021.rds")
AID_blood = dat@phenoData@data$AID
waves <- readRDS("/home/share/preprocessing/preprocessed_two_batches/waves_01.09.2021.rds")

waves = waves %>% mutate(
  color_byinterviewer5 = H3IR17 %>%
    as.character() %>% 
    as.factor %>% 
    fct_recode(
      White = "5",
      Light = "4",
      Medium = "3",
      Dark = "2",
      Black = "1"
    ),
  color_byinterviewer3 = H3IR17 %>%
    as.character() %>% 
    as.factor %>% 
    fct_collapse(
      DarkBlack = c("1", "2"),
      LightMed = c("3", "4"),
      White = "5"),
  immigrat = ifelse(immigrat==0, "Immigrate", "Nonimmigrate") %>% as.factor()) 
# %>% 
#   filter(!is.na(raceethnicity))

#' ## whole sample summary
#+ echo=F, eval=T, warning=FALSE, message=FALSE

data = waves %>%
  dplyr::select(color_byinterviewer5, color_byinterviewer3, immigrat, re)

#' ### skincolor 3 levles
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
descr::crosstab(data$color_byinterviewer3,  data$immigrat, prop.r = T, prop.c = T, prop.t = FALSE, plot = F)

#' ### skincolor 5 levles
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
descr::crosstab(data$color_byinterviewer5,  data$immigrat, prop.r = T, prop.c = T, prop.t = FALSE, plot = F)
#' ## blood  sample summary
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
data_blood = waves %>%
  filter(AID %in% AID_blood) %>% 
  select(color_byinterviewer3, color_byinterviewer5, immigrat, re)

#' ### skincolor 3 levles
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
descr::crosstab(data_blood$color_byinterviewer3,  data_blood$immigrat, prop.r = T, prop.c = T, prop.t = FALSE, plot = F)

#' ### skincolor 5 levles
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
descr::crosstab(data_blood$color_byinterviewer5,  data_blood$immigrat, prop.r = T, prop.c = T, prop.t = FALSE, plot = F)


#' ## blood sample nonhispanic black race group summary
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

nonhispanicblack = data_blood %>% filter(re ==2) %>% filter(color_byinterviewer3!="White")
#' ### skincolor 3 levles
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
descr::crosstab(nonhispanicblack$color_byinterviewer3,  nonhispanicblack$immigrat, prop.r = T, prop.c = T, prop.t = FALSE, plot = F)

#' ### skincolor 5 levles
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

descr::crosstab(nonhispanicblack$color_byinterviewer5,  nonhispanicblack$immigrat, prop.r = T, prop.c = T, prop.t = FALSE, plot = F)

