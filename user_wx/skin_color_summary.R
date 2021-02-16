
#' ---
#' title: skin color and discrimination summary
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

waves <- readRDS("/home/share/preprocessed_two_batches/waves_17.11.2020.rds")

# compete t-tests for these variables by white versus DarkBlack, and also white versus LightMed. 
waves = waves %>% mutate(
  # raceethnicity = re %>%
  #   fct_recode(
  #     White = "1",
  #     # white nonhispanic
  #     Black = "2",
  #     # black nonhispanic
  #     NULL = "3",
  #     # asian nonhispanic
  #     NULL = "4",
  #     # other nonhispanic
  #     Hispanic = "5"
  #   ),
  color_byinterviewer = H3IR17 %>%
    as.character() %>% 
    as.factor %>% 
    fct_collapse(
      DarkBlack = c("1", "2"),
      LightMed = c("3", "4"),
      White = "5") %>%
    relevel(ref = "White")) 
# %>% 
#   filter(!is.na(raceethnicity))

## whole sample summary
#+ echo=F, eval=T, warning=FALSE, message=FALSE

data = waves %>%
  dplyr::select(color_byinterviewer, lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2)

# Hmisc::describe(data %>% 
#                   dplyr::mutate_all(as.factor))

#' ## blood  sample summary
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
data_blood = waves %>%
  filter(AID %in% AID_blood) %>% 
  select(color_byinterviewer, lowbirthweight, high_lowbirth, discrim1, discrim2, totdiscrim1, totdiscrim2)

Hmisc::describe(data_blood %>% 
                  dplyr::mutate_all(as.factor))

#' ## whole sample and blood sample plot
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
A = ggplot(data %>% filter(color_byinterviewer != "NA" & !is.na(lowbirthweight)), 
           aes(x = color_byinterviewer, fill = as.factor(lowbirthweight))) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

B = ggplot(data_blood %>% filter(color_byinterviewer != "NA" & !is.na(lowbirthweight)), 
           aes(x = color_byinterviewer, fill = as.factor(lowbirthweight))) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

ggpubr::ggarrange(A, B, ncol = 2, labels = c("full_sample", "blood_sample"), common.legend = TRUE, legend="right")


A = ggplot(data %>% filter(color_byinterviewer != "NA" & !is.na(high_lowbirth)),
           aes(x = color_byinterviewer, fill = as.factor(high_lowbirth))) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

B = ggplot(data_blood %>% filter(color_byinterviewer != "NA" & !is.na(high_lowbirth)),
           aes(x = color_byinterviewer, fill = as.factor(high_lowbirth))) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

ggpubr::ggarrange(A, B, ncol = 2, labels = c("full_sample", "blood_sample"), common.legend = TRUE, legend="right")


A = ggplot(data %>% filter(color_byinterviewer != "NA" & !is.na(discrim1)), 
           aes(x = discrim1, fill = color_byinterviewer)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

B = ggplot(data_blood %>% filter(color_byinterviewer != "NA" & !is.na(discrim1)),
           aes(x = discrim1, fill = color_byinterviewer)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

ggpubr::ggarrange(A, B, ncol = 2, labels = c("full_sample", "blood_sample"), common.legend = TRUE, legend="right")


A = ggplot(data %>% filter(color_byinterviewer != "NA" & !is.na(discrim2)),
           aes(x = discrim2, fill = color_byinterviewer)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

B = ggplot(data_blood %>% filter(color_byinterviewer != "NA" & !is.na(discrim2)),
           aes(x = discrim2, fill = color_byinterviewer)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

ggpubr::ggarrange(A, B, ncol = 2, labels = c("full_sample", "blood_sample"), common.legend = TRUE, legend="right")

A = ggplot(data %>% filter(color_byinterviewer != "NA" & !is.na(totdiscrim1)),
           aes(x = totdiscrim1, fill = color_byinterviewer)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

B = ggplot(data_blood %>% filter(color_byinterviewer != "NA" & !is.na(totdiscrim1)),
           aes(x = totdiscrim1, fill = color_byinterviewer)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

ggpubr::ggarrange(A, B, ncol = 2, labels = c("full_sample", "blood_sample"), common.legend = TRUE, legend="right")

A = ggplot(data %>% filter(color_byinterviewer != "NA" & !is.na(totdiscrim2)),
           aes(x = totdiscrim2, fill = color_byinterviewer)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

B = ggplot(data_blood %>% filter(color_byinterviewer != "NA" & !is.na(totdiscrim2)),
           aes(x = totdiscrim2, fill = color_byinterviewer)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE)

ggpubr::ggarrange(A, B, ncol = 2, labels = c("full_sample", "blood_sample"), common.legend = TRUE, legend="right")



# whole sample
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
White = data %>% filter(color_byinterviewer=="White") 
DarkBlack = data %>% filter(color_byinterviewer=="DarkBlack")
LightMed =data %>% filter(color_byinterviewer=="LightMed")
#' ## whole sample White DarkBlack and White LightMed 2 sample Wilcoxon test p value
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
a = map2(White %>% select(4:7), DarkBlack %>% select(4:7), wilcox.test) %>%
  map(~ .$p.value) %>% as.data.frame() %>% `rownames<-`("White_DarkBlack")
b = map2(White %>% select(4:7), LightMed %>% select(4:7), wilcox.test) %>% 
  map(~ .$p.value) %>% as.data.frame() %>%  `rownames<-`("White_LightMed")

rbind(a,b)

#' ## whole sample White DarkBlack and White LightMed 2 sample chisquare test p value
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
White_DarkBlack = data %>% dplyr::filter(color_byinterviewer != "LightMed")

a = map(White_DarkBlack %>% select(2:3), chisq.test, White_DarkBlack$color_byinterviewer) %>%
  map(~ .$p.value) %>%
  as.data.frame() %>% 
  `rownames<-`("White_DarkBlack")

White_LightMed = data %>% dplyr::filter(color_byinterviewer != "DarkBlack")

b = map(White_LightMed %>% select(2:3), chisq.test, White_LightMed$color_byinterviewer) %>%
  map(~ .$p.value) %>%
  as.data.frame() %>% 
  `rownames<-`("White_LightMed")
rbind(a,b)
# blood sample
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

White_blood = data_blood %>% filter(color_byinterviewer=="White") 
DarkBlack_blood = data_blood %>% filter(color_byinterviewer=="DarkBlack")
LightMed_blood = data_blood %>% filter(color_byinterviewer=="LightMed")
#' ## blood sample White DarkBlack and White LightMed 2 sample Wilcoxon test p value
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
a = map2(White_blood %>% select(4:7), DarkBlack_blood %>% select(4:7), wilcox.test) %>%
  map(~ .$p.value) %>% as.data.frame() %>% `rownames<-`("White_DarkBlack")
b = map2(White_blood %>% select(4:7), LightMed_blood %>% select(4:7), wilcox.test) %>% 
  map(~ .$p.value) %>% as.data.frame() %>%  `rownames<-`("White_LightMed")

rbind(a,b)

#' ## blood sample White DarkBlack and White LightMed 2 sample chisquare test p value
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
White_DarkBlack = data_blood %>% dplyr::filter(color_byinterviewer != "LightMed")

a = map(White_DarkBlack %>% select(2:3), chisq.test, White_DarkBlack$color_byinterviewer) %>%
  map(~ .$p.value) %>%
  as.data.frame() %>% 
  `rownames<-`("White_DarkBlack")

White_LightMed = data_blood %>% dplyr::filter(color_byinterviewer != "DarkBlack")

b = map(White_LightMed %>% select(2:3), chisq.test, White_LightMed$color_byinterviewer) %>%
  map(~ .$p.value) %>%
  as.data.frame() %>% 
  `rownames<-`("White_LightMed")
rbind(a,b)


#' ##  freq distribution of discrim1, discrim2, totdiscrim1, totdiscrim2, discrimwhy for blood sample
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

countdiscrim = waves %>% 
  select(AID, color_byinterviewer, discrim1, discrim2, totdiscrim1, totdiscrim2, matches("H5MN6[A-L]"), H5MN7) %>% 
  filter(AID %in% AID_blood) %>%
  mutate(countdiscrimwhy = case_when(totdiscrim1 == 0 ~ "no_discrim",
                                     totdiscrim1 == 1 & H5MN7 == 1 ~ "0",
                                     totdiscrim1 > 0 ~ select(., matches("H5MN6[A-L]")) %>% rowSums(na.rm = TRUE) %>% as.character()
  )) %>%
  mutate(totdiscrim1 = totdiscrim1 %>% as.character,
         countdiscrimwhy = factor(countdiscrimwhy,
                                  levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "no_discrim") %>% fct_rev))


#' ### discrim1 is the sum of 5 types of discrimination H5MN5[A-E]
#' * H5MN5[A-E], each range of 0:3, never(0), rarely(1), sometimes(2), often(3)
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

countdiscrim %>%
  filter(!is.na(discrim1)) %>% 
  ggplot(aes(x = discrim1))+
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

countdiscrim %>%
  filter(!is.na(discrim1)) %>% 
  filter(color_byinterviewer=="White") %>% 
  ggplot(aes(x = discrim1))+
  geom_bar(stat = "count", position = "dodge", fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_text(aes(label=..count..),
            stat='count',
            size = 4,
            vjust = 0) +
  geom_text(aes(y = ((..count..)/sum(..count..)), 
                label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
            stat = "count", color = "darkblue", vjust = 1)+
  labs(title="White")

countdiscrim %>%
  filter(!is.na(discrim1)) %>% 
  filter(color_byinterviewer=="DarkBlack") %>% 
  ggplot(aes(x = discrim1))+
  geom_bar(stat = "count", position = "dodge", fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_text(aes(label=..count..),
            stat='count',
            size = 4,
            vjust = 0) +
  geom_text(aes(y = ((..count..)/sum(..count..)), 
                label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
            stat = "count", color = "darkblue", vjust = 1)+
  labs(title="DarkBlack")

countdiscrim %>%
  filter(!is.na(discrim1)) %>% 
  filter(color_byinterviewer=="LightMed") %>% 
  ggplot(aes(x = discrim1))+
  geom_bar(stat = "count", position = "dodge", fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_text(aes(label=..count..),
            stat='count',
            size = 4,
            vjust = 0) +
  geom_text(aes(y = ((..count..)/sum(..count..)), 
                label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
            stat = "count", color = "darkblue", vjust = 1)+
  labs(title="LightMed")



#' ### discrim2 equals discrim1 if they feel the discrimination were due to race and ethnicity, otherwise 0
#' *  discrimination due to race and ethnicity: H5MN6A==1 (ancestry or national origin) or H5MN6D==1(race)
#' *  discrim2 ==0 : either they do not feel any discrimination or they do not think the discrimination were due to race and ethnicity

#+ echo=F, eval=T, warning=FALSE, message=FALSE 

countdiscrim %>%
  filter(!is.na(discrim2)) %>% 
  ggplot(aes(x = discrim2))+
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

countdiscrim %>%
  filter(!is.na(discrim2)) %>% 
  filter(color_byinterviewer=="White") %>% 
  ggplot(aes(x = discrim2))+
  geom_bar(stat = "count", position = "dodge", fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_text(aes(label=..count..),
            stat='count',
            size = 4,
            vjust = 0) +
  geom_text(aes(y = ((..count..)/sum(..count..)), 
                label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
            stat = "count", color = "darkblue", vjust = 1)+
  labs(title="White")

countdiscrim %>%
  filter(!is.na(discrim2)) %>% 
  filter(color_byinterviewer=="DarkBlack") %>% 
  ggplot(aes(x = discrim2))+
  geom_bar(stat = "count", position = "dodge", fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_text(aes(label=..count..),
            stat='count',
            size = 4,
            vjust = 0) +
  geom_text(aes(y = ((..count..)/sum(..count..)), 
                label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
            stat = "count", color = "darkblue", vjust = 1)+
  labs(title="DarkBlack")

countdiscrim %>%
  filter(!is.na(discrim2)) %>% 
  filter(color_byinterviewer=="LightMed") %>% 
  ggplot(aes(x = discrim2))+
  geom_bar(stat = "count", position = "dodge", fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_text(aes(label=..count..),
            stat='count',
            size = 4,
            vjust = 0) +
  geom_text(aes(y = ((..count..)/sum(..count..)), 
                label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
            stat = "count", color = "darkblue", vjust = 1)+
  labs(title="LightMed")




#' ###  totdiscrim1 is the sum of H5MN5[A-E] and H5MN7(binary)
#' * here H5MN5[A-E] are binarized: 0 for never and rarely, 1 for sometimes and often
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
countdiscrim %>%
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

countdiscrim %>%
  filter(!is.na(totdiscrim1)) %>% 
  filter(color_byinterviewer=="White") %>% 
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
  labs(title="White")

countdiscrim %>%
  filter(!is.na(totdiscrim1)) %>% 
  filter(color_byinterviewer=="DarkBlack") %>% 
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
  labs(title="DarkBlack")

countdiscrim %>%
  filter(!is.na(totdiscrim1)) %>% 
  filter(color_byinterviewer=="LightMed") %>% 
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
  labs(title="LightMed")



#' ### totdiscrim2 equals totdiscrim1 if they feel the discrimination were due to race and ethnicity, otherwise 0
#' *  discrimination due to race and ethnicity: H5MN6A==1 (ancestry or national origin) or H5MN6D==1(race)
#' *  discrim2 ==0 : either they do not feel any discrimination or they do not think the discrimination were due to race and ethnicity

#+ echo=F, eval=T, warning=FALSE, message=FALSE 


countdiscrim %>%
  filter(!is.na(totdiscrim2)) %>% 
  ggplot(aes(x = totdiscrim2))+
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

countdiscrim %>%
  filter(!is.na(totdiscrim2)) %>% 
  filter(color_byinterviewer=="White") %>% 
  ggplot(aes(x = totdiscrim2)) +
  geom_bar(stat = "count", position = "dodge", fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_text(aes(label=..count..),
            stat='count',
            size = 4,
            vjust = 0) +
  geom_text(aes(y = ((..count..)/sum(..count..)), 
                label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
            stat = "count", color = "darkblue", vjust = 1)+
  labs(title="White")

countdiscrim %>%
  filter(!is.na(totdiscrim2)) %>% 
  filter(color_byinterviewer=="DarkBlack") %>% 
  ggplot(aes(x = totdiscrim2))+
  geom_bar(stat = "count", position = "dodge", fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_text(aes(label=..count..),
            stat='count',
            size = 4,
            vjust = 0) +
  geom_text(aes(y = ((..count..)/sum(..count..)), 
                label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
            stat = "count", color = "darkblue", vjust = 1)+
  labs(title="DarkBlack")

countdiscrim %>%
  filter(!is.na(totdiscrim2)) %>% 
  filter(color_byinterviewer=="LightMed") %>% 
  ggplot(aes(x = totdiscrim2))+
  geom_bar(stat = "count", position = "dodge", fill = "lightblue") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_text(aes(label=..count..),
            stat='count',
            size = 4,
            vjust = 0) +
  geom_text(aes(y = ((..count..)/sum(..count..)), 
                label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
            stat = "count", color = "darkblue", vjust = 1)+
  labs(title="LightMed")


#' ###  countdiscrimwhy is created: if totdiscrim1>0, sum "H5MN6[A-L]" (reasons of discrimination)
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
countdiscrim %>%
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

countdiscrim %>%
  filter(!is.na(countdiscrimwhy)) %>%
  filter(countdiscrimwhy!="no_discrim") %>% 
  filter(color_byinterviewer=="White") %>% 
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
  labs(title="White")

countdiscrim %>%
  filter(!is.na(countdiscrimwhy)) %>%
  filter(countdiscrimwhy!="no_discrim") %>% 
  filter(color_byinterviewer=="DarkBlack") %>% 
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
  labs(title="DarkBlack")

countdiscrim %>%
  filter(!is.na(countdiscrimwhy)) %>%
  filter(countdiscrimwhy!="no_discrim") %>% 
  filter(color_byinterviewer=="LightMed") %>% 
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
  labs(title="LightMed")


