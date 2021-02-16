
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
  dplyr::select(raceethnicity, lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2)

Hmisc::describe(data %>% 
                  dplyr::mutate_all(as.factor))

#' ## blood  sample summary
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
data_blood = waves %>%
  filter(AID %in% AID_blood) %>% 
  select(raceethnicity, lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2)

Hmisc::describe(data_blood %>% 
                  dplyr::mutate_all(as.factor))

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
white = data %>% filter(raceethnicity=="NonHwhite") 
black = data %>% filter(raceethnicity=="NonHblack")
hispanic =data %>% filter(raceethnicity=="Hispanic")
#' ## whole sample white black and white hispanic 2 sample Wilcoxon test p value
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
a = map2(white %>% select(4:7), black %>% select(4:7), wilcox.test) %>% map(~ .$p.value) %>% as.data.frame() %>% `rownames<-`("white_black")
b = map2(white %>% select(4:7), hispanic %>% select(4:7), wilcox.test) %>% map(~ .$p.value) %>% as.data.frame() %>%  `rownames<-`("white_hispanic")

rbind(a,b)

#' ## whole sample white black and white hispanic 2 sample chisquare test p value
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
white_black = data %>% dplyr::filter(raceethnicity != "Hispanic")

a = map(white_black %>% select(2:3), chisq.test, white_black$raceethnicity) %>%
  map(~ .$p.value) %>%
  as.data.frame() %>% 
  `rownames<-`("white_black")

white_hispanic = data %>% dplyr::filter(raceethnicity != "NonHblack")

b = map(white_hispanic %>% select(2:3), chisq.test, white_hispanic$raceethnicity) %>%
  map(~ .$p.value) %>%
  as.data.frame() %>% 
  `rownames<-`("white_hispanic")
rbind(a,b)
# blood sample
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

white_blood = data_blood %>% filter(raceethnicity=="NonHwhite") 
black_blood = data_blood %>% filter(raceethnicity=="NonHblack")
hispanic_blood = data_blood %>% filter(raceethnicity=="Hispanic")
#' ## blood sample white black and white hispanic 2 sample Wilcoxon test p value
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
a = map2(white_blood %>% select(4:7), black_blood %>% select(4:7), wilcox.test) %>% map(~ .$p.value) %>% as.data.frame() %>% `rownames<-`("white_black")
b = map2(white_blood %>% select(4:7), hispanic_blood %>% select(4:7), wilcox.test) %>% map(~ .$p.value) %>% as.data.frame() %>%  `rownames<-`("white_hispanic")

rbind(a,b)

#' ## blood sample white black and white hispanic 2 sample chisquare test p value
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
white_black = data_blood %>% dplyr::filter(raceethnicity != "Hispanic")

a = map(white_black %>% select(2:3), chisq.test, white_black$raceethnicity) %>%
  map(~ .$p.value) %>%
  as.data.frame() %>% 
  `rownames<-`("white_black")

white_hispanic = data_blood %>% dplyr::filter(raceethnicity != "NonHblack")

b = map(white_hispanic %>% select(2:3), chisq.test, white_hispanic$raceethnicity) %>%
  map(~ .$p.value) %>%
  as.data.frame() %>% 
  `rownames<-`("white_hispanic")
rbind(a,b)


#' ##  freq distribution for countdiscrim (sometimes, often) in blood sample
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

countdiscrim = waves %>% select(AID, raceethnicity, totdiscrim1, matches("H5MN6[A-L]"), H5MN7) %>% 
  filter(AID %in% AID_blood) %>%
  mutate(countdiscrimwhy = case_when(totdiscrim1 == 0 ~ "no_discrim",
                                     totdiscrim1 == 1 & H5MN7 == 1 ~ "0",
                                     totdiscrim1 > 0 ~ select(., matches("H5MN6[A-L]")) %>% rowSums(na.rm = TRUE) %>% as.character()
                                    )) %>%
  mutate(totdiscrim1 = totdiscrim1 %>% as.character,
         countdiscrimwhy = factor(countdiscrimwhy,
                                  levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "no_discrim") %>% fct_rev))

#' ##  totdiscrim1
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

a = countdiscrim %>% select(raceethnicity, totdiscrim1) %>% table()
whole = countdiscrim %>% select(totdiscrim1) %>% table()

a %>% prop.table() 

rbind(a,whole)

#' ##  countdiscrimwhy
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

a = countdiscrim %>% select(raceethnicity, countdiscrimwhy) %>% table()

a %>% prop.table() %>% t
whole = countdiscrim %>% select(countdiscrimwhy) %>% table()
rbind(a,whole)

#' ###  count totdiscrim1 freq distribution
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
a = countdiscrim %>%
  filter(!is.na(totdiscrim1)) %>% 
  ggplot(aes(x = totdiscrim1))+
  geom_bar(stat = "count", position = "dodge", fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_text(aes(label=..count..),
            stat='count',
            size = 4,
            vjust = -0) +
  geom_text(aes(y = ((..count..)/sum(..count..)), 
                label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
            stat = "count", color = "darkblue", vjust = 1)+
  labs(title="whole blood")
a
w = countdiscrim %>%
  filter(!is.na(totdiscrim1)) %>% 
  filter(raceethnicity=="NonHwhite") %>% 
  ggplot(aes(x = totdiscrim1))+
  geom_bar(stat = "count", position = "dodge", fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_text(aes(label=..count..),
            stat='count',
            size = 4,
            vjust = 0) +
  geom_text(aes(y = ((..count..)/sum(..count..)), 
                label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
            stat = "count", color = "darkblue", vjust = 1)+
  labs(title="white")
w
b = countdiscrim %>%
  filter(!is.na(totdiscrim1)) %>% 
  filter(raceethnicity=="NonHblack") %>% 
  ggplot(aes(x = totdiscrim1))+
  geom_bar(stat = "count", position = "dodge", fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_text(aes(label=..count..),
            stat='count',
            size = 4,
            vjust = 0) +
  geom_text(aes(y = ((..count..)/sum(..count..)), 
                label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
            stat = "count", color = "darkblue", vjust = 1)+
  labs(title="black")
b
h = countdiscrim %>%
  filter(!is.na(totdiscrim1)) %>% 
  filter(raceethnicity=="Hispanic") %>% 
  ggplot(aes(x = totdiscrim1))+
  geom_bar(stat = "count", position = "dodge", fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_text(aes(label=..count..),
            stat='count',
            size = 4,
            vjust = 0) +
  geom_text(aes(y = ((..count..)/sum(..count..)), 
                label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
            stat = "count", color = "darkblue", vjust = 1)+
  labs(title="hispanic")
h
# ggpubr::ggarrange(a, w, b, h, ncol = 2, nrow = 2, labels = c("blood_sample", "white", "black", "hispanic"), common.legend = TRUE, legend="right") 




#' ###  countdiscrimwhy freq distribution
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
a = countdiscrim %>%
  filter(!is.na(countdiscrimwhy)) %>%
  filter(countdiscrimwhy!="no_discrim") %>% 
  ggplot(aes(x = countdiscrimwhy)) +
  geom_bar(stat = "count", position = "dodge", fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_text(aes(label=..count..),
            stat='count',
            size = 4,
            vjust = 0) +
  geom_text(aes(y = ((..count..)/sum(..count..)), 
                label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
            stat = "count", color = "darkblue", vjust = 1, size =3)+
  labs(title="whole blood")
a
w = countdiscrim %>%
  filter(!is.na(countdiscrimwhy)) %>%
  filter(countdiscrimwhy!="no_discrim") %>% 
  filter(raceethnicity=="NonHwhite") %>% 
  ggplot(aes(x = countdiscrimwhy)) +
  geom_bar(stat = "count", position = "dodge", fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_text(aes(label=..count..),
            stat='count',
            size = 4,
            vjust = 0) +
  geom_text(aes(y = ((..count..)/sum(..count..)), 
                label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
            stat = "count", color = "darkblue", vjust = 1, size =3)+
  labs(title="white")
w
b = countdiscrim %>%
  filter(!is.na(countdiscrimwhy)) %>%
  filter(countdiscrimwhy!="no_discrim") %>% 
  filter(raceethnicity=="NonHblack") %>% 
  ggplot(aes(x = countdiscrimwhy)) +
  geom_bar(stat = "count", position = "dodge", fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_text(aes(label=..count..),
            stat='count',
            size = 4,
            vjust = 0) +
  geom_text(aes(y = ((..count..)/sum(..count..)), 
                label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
            stat = "count", color = "darkblue", vjust = 1, size =3)+
  labs(title="black")
b
h = countdiscrim %>%
  filter(!is.na(countdiscrimwhy)) %>%
  filter(countdiscrimwhy!="no_discrim") %>% 
  filter(raceethnicity=="Hispanic") %>% 
  ggplot(aes(x = countdiscrimwhy)) +
  geom_bar(stat = "count", position = "dodge", fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_text(aes(label=..count..),
            stat='count',
            size = 4,
            vjust = 0) +
  geom_text(aes(y = ((..count..)/sum(..count..)), 
                label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
            stat = "count", color = "darkblue", vjust = 1, size =3)+
  labs(title="hispanic")
h
# ggpubr::ggarrange(a, w, b, h, ncol = 2, nrow = 2, labels = c("blood_sample", "white", "black", "hispanic"), common.legend = TRUE, legend="right") 

