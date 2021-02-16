
#' ---
#' title: race and discrimination summary
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE


library(tidyverse)
dat <- readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_29.10.2020.rds")
AID_blood = dat@phenoData@data$AID

waves <- readRDS("/home/share/preprocessed_two_batches/waves_03.11.2020.rds")

# compete t-tests for these variables by white versus black, and also white versus hispanic. 
# not appropriate as these variables are not approximately normal distributed
 
# whole sample
white = waves %>% filter(re==1) %>% select(lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2) 
 black = waves %>% filter(re==2) %>% select(lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2) 
 hispanic = waves %>% filter(re==5) %>% select(lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2) 

#' ## white black 2 sample t test 
 map2(white, black, t.test) %>% map(~ .$p.value) %>% as.data.frame() %>% `rownames<-`("white_black")

#' ## white hispanic 2 sample t test 
 map2(white,hispanic, t.test) %>% map(~ .$p.value) %>% as.data.frame() %>%  `rownames<-`("white_hispanic")
# blood sample
 white_blood = waves %>% filter(re==1, AID %in% AID_blood) %>% select(lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2) 
 black_blood = waves %>% filter(re==2, AID %in% AID_blood) %>% select(lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2) 
 hispanic_blood = waves %>% filter(re==5, AID %in% AID_blood) %>% select(lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2) 
 
#' ## white black 2 sample t test 
 map2(white_blood, black_blood, t.test) %>% map(~ .$p.value) %>% as.data.frame() %>% `rownames<-`("white_black")
 #' ## white hispanic 2 sample t test
 map2(white_blood,hispanic_blood, t.test) %>% map(~ .$p.value) %>% as.data.frame() %>%  `rownames<-`("white_hispanic")
 # chi square test for whole sample
 waves = waves %>% mutate(
         raceethnicity = re %>%
         fct_recode(
                 NonHwhite = "1",
                 # white nonhispanic
                 NonHblack = "2",
                 # black nonhispanic
                 NULL = "3",
                 # asian nonhispanic
                 NULL = "4",
                 # other nonhispanic
                 Hispanic = "5"
         ))
 
#' ## whole sample summary
#+ echo=F, eval=T, warning=FALSE, message=FALSE
 
data = waves %>%
        dplyr::select(raceethnicity, lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2) %>% 
        dplyr::mutate_all(as.factor)
 
Hmisc::describe(data)

#' ## white black chisquare test (asymptotic)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
white_black = data %>% dplyr::filter(raceethnicity != "Hispanic")

 map(white_black %>% select(-1), chisq.test, white_black$raceethnicity) %>%
         map(~ .$p.value) %>%
         as.data.frame() %>% 
         `rownames<-`("white_black")


#' ## white black chisquare test (simulated p value)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

 map(white_black %>% select(-1), chisq.test, white_black$raceethnicity, simulate.p.value=TRUE) %>%
         map(~ .$p.value) %>%
         as.data.frame() %>% 
         `rownames<-`("white_black")


#' ## white hispanic chisquare test (asymptotic)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
white_hispanic = data %>% dplyr::filter(raceethnicity != "NonHblack")

 map(white_hispanic %>% select(-1), chisq.test, white_hispanic$raceethnicity) %>%
         map(~ .$p.value) %>%
         as.data.frame() %>% 
         `rownames<-`("white_hispanic")

#' ## white hispanic chisquare test (simulated p value)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

map(white_hispanic %>% select(-1), chisq.test, white_hispanic$raceethnicity, simulate.p.value=TRUE) %>%
        map(~ .$p.value) %>%
        as.data.frame() %>% 
        `rownames<-`("white_hispanic")
#' ## blood sample summary
#+ echo=F, eval=T, warning=FALSE, message=FALSE

data_blood = waves %>%
        filter(AID %in% AID_blood) %>% 
        select(raceethnicity, lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2) %>% 
        dplyr::mutate_all(as.factor)
Hmisc::describe(data_blood)


#' ## white black chisquare test (asymptotic)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
white_black = data_blood %>% dplyr::filter(raceethnicity != "Hispanic")

map(white_black %>% select(-1), chisq.test, white_black$raceethnicity) %>%
        map(~ .$p.value) %>%
        as.data.frame() %>% 
        `rownames<-`("white_black")


#' ## white black chisquare test (simulated p value)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

map(white_black %>% select(-1), chisq.test, white_black$raceethnicity, simulate.p.value=TRUE) %>%
        map(~ .$p.value) %>%
        as.data.frame() %>% 
        `rownames<-`("white_black")


#' ## white hispanic chisquare test (asymptotic)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
white_hispanic = data_blood %>% dplyr::filter(raceethnicity != "NonHblack")

map(white_hispanic %>% select(-1), chisq.test, white_hispanic$raceethnicity) %>%
        map(~ .$p.value) %>%
        as.data.frame() %>% 
        `rownames<-`("white_hispanic")

#' ## white hispanic chisquare test (simulated p value)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

map(white_hispanic %>% select(-1), chisq.test, white_hispanic$raceethnicity, simulate.p.value=TRUE) %>%
        map(~ .$p.value) %>%
        as.data.frame() %>% 
        `rownames<-`("white_hispanic")





