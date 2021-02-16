
#' ---
#' title: race ethnicity and skin color summary
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

#+ echo=F, eval=T, warning=FALSE, message=FALSE 
library(tidyverse)
dat <- readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_17.11.2020.rds")
AID_blood = dat@phenoData@data$AID

waves <- readRDS("/home/share/preprocessed_two_batches/waves_17.11.2020.rds")

# compete t-tests for these variables by white versus black, and also white versus hispanic. 
waves = waves %>% mutate(
  raceethnicity = re %>%
    fct_recode(
      White = "1",
      # white nonhispanic
      Black = "2",
      # black nonhispanic
      Asia = "3",
      # asian nonhispanic
      Othernonhisp = "4",
      # other nonhispanic
      Hispanic = "5"
    ),
  color_byinterviewer = H3IR17 %>%
    as.character() %>% 
    as.factor %>% 
    fct_recode(
      black = "1",
      dark_brown = "2",
      medium_brown = "3",
      light_brown = "4",
      white = "5"
    ),
  mRNAsample = ifelse(AID %in% AID_blood, "mRNA", "NonmRNA"))

#' ## whole sample summary
#+ echo=F, eval=T, warning=FALSE, message=FALSE, results="asis" 

data = waves %>%
  dplyr::select(mRNAsample, raceethnicity, color_byinterviewer) 



tab3 <- arsenal::tableby(~ raceethnicity + color_byinterviewer, data=data, test=FALSE)
summary(tab3)
#' ## mRNA  sample summary
#+ echo=F, eval=T, warning=FALSE, message=FALSE, results="asis"  
tab2 <- arsenal::tableby(~ raceethnicity + color_byinterviewer, data=data %>% filter(mRNAsample=="mRNA"), test=FALSE)
summary(tab2)


tab1 <- arsenal::tableby(mRNAsample ~ raceethnicity + color_byinterviewer, data=data, test=FALSE)
summary(tab1)
# reporttools::tableNominal(vars = data[,2:3], group = data[,1], cap = 
#                +     "Table of nominal variables.", lab = "tab: nominal")
# 
# stargazer::stargazer(data ~ AID_blood, type = "html")
