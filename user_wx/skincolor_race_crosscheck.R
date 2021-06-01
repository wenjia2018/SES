
#' ---
#' title: skin color race cross check
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE
#' ## skin color race cross check
#' 
#' * skin color by interviewer is only available at wave 3: H3IR17
#' * race by interviewer are available at wave 1(H1GI9), wave 3(H3IR4) wave 4(H4IR4)
#' * self reported race are available at wave 1(H1GI6[A-E]), wave 3(H3OD4[A-E]) wave 5(H5OD4[A-E])
#' * our raceethnicity variable is created by self reported and interviewer coded together, i.e. 
#' use wave 1 selfreported race, if NA, use wave 3 selfreported, if NA, wave 5 selfreported, if NA, then
#' interviewer coded race w1 w3 w4 sequentially to best minimize NA in this variables 
#' 
#' 
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
library(tidyverse)
dat <- readRDS("/home/share/preprocessing/preprocessed_two_batches/dt_batches1_2_steve_waves_29.10.2020.rds")
AID_blood = dat@phenoData@data$AID

waves <- readRDS("/home/share/preprocessing/preprocessed_two_batches/waves_17.11.2020.rds")

dt = waves %>%
  filter(AID %in% AID_blood) %>% 
  select(race_interv, H3IR17, H3IR4, H4IR4, re) %>%
  mutate_at(.vars = vars(race_interv), .funs = list(~ .x %>% 
                                                      fct_recode("white" = "w", 
                                                                 "black" = "b", 
                                                                 "native" = "n", 
                                                                 "asian"  = "a",  
                                                                 "other"  = "o",  
                                                                 ))) %>% 
  mutate_at(.vars = vars(c(H3IR4, H4IR4)), .funs = list(~ .x %>% 
                                                          factor %>% 
                                                          fct_recode("white" = "1", 
                                                                     "black" = "2", 
                                                                     "native" = "3",
                                                                     "asian"  = "4"))) %>% 
  rename(race_interv_w1 = race_interv,
         race_interv_w3 = H3IR4,
         race_interv_w4 = H4IR4) %>% 
  mutate(color_byinterviewer = H3IR17 %>%
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
             White = "5") %>%
           relevel(ref = "White"),
         raceethnicity = re %>%
             fct_recode(
               `Non Hispanic White` = "1",
               # white nonhispanic
               `Non Hispanic Black` = "2",
               # black nonhispanic
               NULL = "3",
               # asian nonhispanic
               NULL = "4",
               # other nonhispanic
               Hispanic = "5"
             ))
  

  
  

#' ## race by interviewer from wave 1 and wave 3
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
descr::crosstab(dt$race_interv_w1, dt$race_interv_w3,  prop.r = T, prop.c = T, prop.t = FALSE, plot = F)

#' ## race by interviewer from wave 1 and wave 4
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

descr::crosstab(dt$race_interv_w1, dt$race_interv_w4,  prop.r = T, prop.c = T, prop.t = FALSE, plot = F)

#' ## race by interviewer from wave 3 and wave 4
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
descr::crosstab(dt$race_interv_w3, dt$race_interv_w4,  prop.r = T, prop.c = T, prop.t = FALSE, plot = F)


#' ## race by interviewer from wave 1 and skin color by interviewer from wave 3
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
descr::crosstab(dt$race_interv_w1, dt$color_byinterviewer, prop.r = T, prop.c = T, prop.t = FALSE, plot = F)

#' ## race by interviewer from wave 3 and skin color by interviewer from wave 3
#+ echo=F, eval=T, warning=FALSE, message=FALSE 


descr::crosstab(dt$race_interv_w3, dt$color_byinterviewer, prop.r = T, prop.c = T, prop.t = FALSE, plot = F)


#' ## race by interviewer from wave 4 and skin color by interviewer from wave 3
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

descr::crosstab(dt$race_interv_w4, dt$color_byinterviewer, prop.r = T, prop.c = T, prop.t = FALSE, plot = F)


#' ## created raceethnicity and race by interviewer from wave 1
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

descr::crosstab(dt$raceethnicity, dt$race_interv_w1,  prop.r = T, prop.c = T, prop.t = FALSE, plot = F)

#' ## created raceethnicity and race by interviewer from wave 3
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

descr::crosstab(dt$raceethnicity, dt$race_interv_w3,  prop.r = T, prop.c = T, prop.t = FALSE, plot = F)




#' ## created raceethnicity and race by interviewer from wave 4
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

descr::crosstab(dt$raceethnicity, dt$race_interv_w4,  prop.r = T, prop.c = T, prop.t = FALSE, plot = F)



#' ## created raceethnicity and skin color by interviewer from wave 3
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

descr::crosstab(dt$raceethnicity, dt$color_byinterviewer,  prop.r = T, prop.c = T, prop.t = FALSE, plot = F)

#' ## created raceethnicity and skin color by interviewer(3 levels) from wave 3
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

ct = descr::crosstab(dt$raceethnicity, dt$color_byinterviewer3,  prop.r = T, prop.c = T, prop.t = FALSE, plot = F, row.labels = T)
ct
# ct %>% xtable %>% openxlsx::write.xlsx("./user_wx/testtable.xlsx")
# testtable <- read_excel("user_wx/testtable.xlsx")


sjPlot::tab_xtab(dt$raceethnicity, dt$color_byinterviewer3,
                 show.cell.prc = FALSE, show.summary = FALSE,
                 title = "Race ethnicity and skin color contingency table",
                 var.labels = c("Race Ethnicity", "Skin Color"),
                 emph.total = TRUE,
                 show.row.prc = TRUE, 
                 show.col.prc = TRUE)
