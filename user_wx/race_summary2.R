
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

#+ echo=F, eval=T, warning=FALSE, message=FALSE 
library(tidyverse)
dat <- readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_29.10.2020.rds")
AID_blood = dat@phenoData@data$AID

waves <- readRDS("/home/share/preprocessed_two_batches/waves_03.11.2020.rds")

# compete t-tests for these variables by white versus black, and also white versus hispanic. 
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
    )) %>% 
  filter(!is.na(raceethnicity))

#' ## whole sample summary
#+ echo=F, eval=T, warning=FALSE, message=FALSE

data = waves %>%
  dplyr::select(raceethnicity, lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2) %>% 
  dplyr::mutate_all(as.factor)

Hmisc::describe(data)

#' ## blood  sample summary
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
data_blood = waves %>%
  filter(AID %in% AID_blood) %>% 
  select(raceethnicity, lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2) %>% 
  dplyr::mutate_all(as.factor)
Hmisc::describe(data_blood)

#' ## whole sample and blood sample plot
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
A = ggplot(data %>% filter(raceethnicity != "NA" & !is.na(lowbirthweight)), aes(x = raceethnicity, fill = lowbirthweight)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

B = ggplot(data_blood %>% filter(raceethnicity != "NA" & !is.na(lowbirthweight)), aes(x = raceethnicity, fill = lowbirthweight)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

ggpubr::ggarrange(A, B, ncol = 2, labels = c("full_sample", "blood_sample"), common.legend = TRUE, legend="right")


A = ggplot(data %>% filter(raceethnicity != "NA" & !is.na(high_lowbirth)), aes(x = raceethnicity, fill = high_lowbirth)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

B = ggplot(data_blood %>% filter(raceethnicity != "NA" & !is.na(high_lowbirth)), aes(x = raceethnicity, fill = high_lowbirth)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

ggpubr::ggarrange(A, B, ncol = 2, labels = c("full_sample", "blood_sample"), common.legend = TRUE, legend="right")


A = ggplot(data %>% filter(raceethnicity != "NA" & !is.na(discrim1)), aes(x = discrim1, fill = raceethnicity)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

B = ggplot(data_blood %>% filter(raceethnicity != "NA" & !is.na(discrim1)), aes(x = discrim1, fill = raceethnicity)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

ggpubr::ggarrange(A, B, ncol = 2, labels = c("full_sample", "blood_sample"), common.legend = TRUE, legend="right")


A = ggplot(data %>% filter(raceethnicity != "NA" & !is.na(discrim2)), aes(x = discrim2, fill = raceethnicity)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

B = ggplot(data_blood %>% filter(raceethnicity != "NA" & !is.na(discrim2)), aes(x = discrim2, fill = raceethnicity)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

ggpubr::ggarrange(A, B, ncol = 2, labels = c("full_sample", "blood_sample"), common.legend = TRUE, legend="right")

A = ggplot(data %>% filter(raceethnicity != "NA" & !is.na(totdiscrim1)), aes(x = totdiscrim1, fill = raceethnicity)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

B = ggplot(data_blood %>% filter(raceethnicity != "NA" & !is.na(totdiscrim1)), aes(x = totdiscrim1, fill = raceethnicity)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

ggpubr::ggarrange(A, B, ncol = 2, labels = c("full_sample", "blood_sample"), common.legend = TRUE, legend="right")

A = ggplot(data %>% filter(raceethnicity != "NA" & !is.na(totdiscrim2)), aes(x = totdiscrim2, fill = raceethnicity)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

B = ggplot(data_blood %>% filter(raceethnicity != "NA" & !is.na(totdiscrim2)), aes(x = totdiscrim2, fill = raceethnicity)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

ggpubr::ggarrange(A, B, ncol = 2, labels = c("full_sample", "blood_sample"), common.legend = TRUE, legend="right")



# whole sample
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
white = waves %>% filter(re==1) %>% select(lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2) 
black = waves %>% filter(re==2) %>% select(lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2) 
hispanic = waves %>% filter(re==5) %>% select(lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2) 

#' ## whole sample white black and white hispanic 2 sample Wilcoxon test p value
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
a = map2(white, black, wilcox.test) %>% map(~ .$p.value) %>% as.data.frame() %>% `rownames<-`("white_black")
b = map2(white,hispanic, wilcox.test) %>% map(~ .$p.value) %>% as.data.frame() %>%  `rownames<-`("white_hispanic")

rbind(a,b)

#' ## whole sample white black and white hispanic 2 sample chisquare test p value
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
white_black = data %>% dplyr::filter(raceethnicity != "Hispanic")

map(white_black %>% select(-1), chisq.test, white_black$raceethnicity) %>%
  map(~ .$p.value) %>%
  as.data.frame() %>% 
  `rownames<-`("white_black")

rbind(a,b)



# blood sample
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

white_blood = waves %>% filter(re==1, AID %in% AID_blood) %>% select(lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2) 
black_blood = waves %>% filter(re==2, AID %in% AID_blood) %>% select(lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2) 
hispanic_blood = waves %>% filter(re==5, AID %in% AID_blood) %>% select(lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2) 

#' ## blood sample white black and white hispanic 2 sample Wilcoxon test p value
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
a = map2(white_blood, black_blood, wilcox.test) %>% map(~ .$p.value) %>% as.data.frame() %>% `rownames<-`("white_black")

b = map2(white_blood,hispanic_blood, wilcox.test) %>% map(~ .$p.value) %>% as.data.frame() %>%  `rownames<-`("white_hispanic")

rbind(a,b)


#' ##  freq distribution for countdiscrim in whole sample
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

countdiscrim = waves %>% select(AID, re, matches("H5MN5[A-E]")) %>% 
  mutate_at(vars(matches("H5MN5[A-E]")),
            .funs = list(new = ~ case_when(. %in% c(2, 3, 4) ~ 1,
                                           . %in% c(1) ~ 0))) %>% 
  rename( "less_respect" = H5MN5A_new,
          "poorer_service" = H5MN5B_new,
          "not_smart" = H5MN5C_new,
          "afraid_of" = H5MN5D_new,
          "threatened" = H5MN5E_new) %>%
  filter(re %in% c(1,2,5))


a = countdiscrim %>% select(8:12) %>% 
  gather(
    key = "discrim_experience",
    value = "discrim_count"
  ) %>%
  group_by(discrim_experience) %>%
  summarise_at("discrim_count", sum, na.rm = TRUE) %>% 
  mutate(freq = discrim_count/dim(countdiscrim)[1])
a
# countdiscrim  %>% filter(re %in% c(1,2,5)) %>% count(less_respect)


#' ###  countdiscrim freq distribution for white
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

w = countdiscrim  %>% filter(re == 1) %>% select(8:12) %>% 
  gather(
  key = "discrim_experience",
  value = "discrim_count"
) %>%
  group_by(discrim_experience) %>%
  summarise_at("discrim_count", sum, na.rm = TRUE) %>% 
  mutate(freq = discrim_count/dim(countdiscrim%>% filter(re == 1))[1])

w
#' ###  countdiscrim freq distribution for black
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

b = countdiscrim  %>% filter(re == 2) %>% select(8:12) %>% 
  gather(
    key = "discrim_experience",
    value = "discrim_count"
  ) %>%
  group_by(discrim_experience) %>%
  summarise_at("discrim_count", sum, na.rm = TRUE) %>% 
  mutate(freq = discrim_count/dim(countdiscrim%>% filter(re == 2))[1])
b
#' ###  countdiscrim freq distribution for hispanic
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
h = countdiscrim  %>% filter(re == 5) %>% select(8:12) %>% 
  gather(
    key = "discrim_experience",
    value = "discrim_count"
  ) %>%
  group_by(discrim_experience) %>%
  summarise_at("discrim_count", sum, na.rm = TRUE) %>% 
  mutate(freq = discrim_count/dim(countdiscrim%>% filter(re == 5))[1])
h

countdiscrim %>% select(2, 8:12) %>% 
  gather(
    key = "discrim_experience",
    value = "discrim_count",
    -1) %>%
  mutate(
    re = re %>%
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
      )   
  ) %>% 
  group_by(re, discrim_experience) %>%
  summarise_at("discrim_count", sum, na.rm = TRUE) %>% 
  mutate(freq = discrim_count/dim(countdiscrim)[1]) %>% 
  ggplot(aes(x = discrim_experience, y = discrim_count)) + 
  geom_bar(stat = "identity") +
  facet_wrap("re") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


#' ##  freq distribution for countdiscrim in blood sample
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

countdiscrim = waves %>% select(AID, re, matches("H5MN5[A-E]")) %>% 
  filter(AID %in% AID_blood) %>%
  mutate_at(vars(matches("H5MN5[A-E]")),
            .funs = list(new = ~ case_when(. %in% c(2, 3, 4) ~ 1,
                                           . %in% c(1) ~ 0))) %>% 
  rename( "less_respect" = H5MN5A_new,
          "poorer_service" = H5MN5B_new,
          "not_smart" = H5MN5C_new,
          "afraid_of" = H5MN5D_new,
          "threatened" = H5MN5E_new) %>%
  filter(re %in% c(1,2,5))


# countdiscrim %>% select(8:12) %>%
#   summarise_all(sum, na.rm = TRUE) %>% as.data.frame %>%
#   `rownames<-`("count")  %>% t %>% as.data.frame %>%  rownames_to_column("discrim") %>% 
#   mutate(freq = count/dim(countdiscrim)[1] )

# countdiscrim  %>% filter(re %in% c(1,2,5)) %>% count(less_respect)

a = countdiscrim %>% select(8:12) %>% 
  gather(
       key = "discrim_experience",
       value = "discrim_count"
       ) %>%
  group_by(discrim_experience) %>%
  summarise_at("discrim_count", sum, na.rm = TRUE) %>% 
  mutate(freq = discrim_count/dim(countdiscrim)[1])
 a 
# countdiscrim %>% select(8:12) %>% 
#   gather(
#     key = "discrim_experience",
#     value = "discrim_count"
#   ) %>%
#   group_by(discrim_experience) %>%
#   summarise_at("discrim_count", sum, na.rm = TRUE) %>% 
#   mutate(freq = discrim_count/dim(countdiscrim)[1])%>% 
#   ggplot(aes(x  = discrim_experience, y = discrim_count)) +
#   geom_bar(position = "stack", stat = "identity", na.rm = TRUE)

#' ###  countdiscrim freq distribution for white
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

w = countdiscrim  %>% filter(re == 1) %>% select(8:12) %>% 
  gather(
    key = "discrim_experience",
    value = "discrim_count"
  ) %>%
  group_by(discrim_experience) %>%
  summarise_at("discrim_count", sum, na.rm = TRUE) %>% 
  mutate(freq = discrim_count/dim(countdiscrim %>% filter(re == 1))[1])

w

#' ###  countdiscrim freq distribution for black
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

b = countdiscrim  %>% filter(re == 2) %>% select(8:12) %>% 
  gather(
    key = "discrim_experience",
    value = "discrim_count"
  ) %>%
  group_by(discrim_experience) %>%
  summarise_at("discrim_count", sum, na.rm = TRUE) %>% 
  mutate(freq = discrim_count/dim(countdiscrim%>% filter(re == 2))[1])
b
#' ###  countdiscrim freq distribution for hispanic
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
h = countdiscrim  %>% filter(re == 5) %>% select(8:12) %>% 
  gather(
    key = "discrim_experience",
    value = "discrim_count"
  ) %>%
  group_by(discrim_experience) %>%
  summarise_at("discrim_count", sum, na.rm = TRUE) %>% 
  mutate(freq = discrim_count/dim(countdiscrim %>% filter(re == 5))[1])
h

countdiscrim %>% select(2, 8:12) %>% 
  gather(
    key = "discrim_experience",
    value = "discrim_count",
    -1) %>%
  mutate(
    re = re %>%
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
      )   
  ) %>% 
  group_by(re, discrim_experience) %>%
  summarise_at("discrim_count", sum, na.rm = TRUE) %>% 
  mutate(freq = discrim_count/dim(countdiscrim)[1]) %>% 
  ggplot(aes(x = discrim_experience, y = discrim_count)) + 
  geom_bar(stat = "identity") +
  facet_wrap("re") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
