
#' ---
#' title: mediator summary
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
dat <- readRDS("/home/share/preprocessing/preprocessed_two_batches/all.batches.expression.set.tmm_waves_01.09.2021.rds")
AID_blood = dat@phenoData@data$AID

waves <- readRDS("/home/share/preprocessing/preprocessed_two_batches/waves_01.09.2021.rds")

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



data_blood = waves %>%
  filter(AID %in% AID_blood) %>% 
  select(raceethnicity, discrim1, discrim2, totdiscrim1, totdiscrim2)



#' ##  freq distribution of discrim1, discrim2, totdiscrim1, totdiscrim2, discrimwhy for blood sample
#+ echo=F, eval=T, warning=FALSE, message=FALSE 

countdiscrim = waves %>% select(AID, raceethnicity, discrim1, discrim2, totdiscrim1, totdiscrim2, matches("H5MN6[A-L]"), H5MN7) %>% 
  filter(AID %in% AID_blood) %>%
  mutate(countdiscrimwhy = case_when(totdiscrim1 == 0 ~ "no_discrim",
                                     totdiscrim1 == 1 & H5MN7 == 1 ~ "0",
                                     totdiscrim1 > 0 ~ select(., matches("H5MN6[A-L]")) %>% rowSums(na.rm = TRUE) %>% as.character()
  )) %>%
  mutate(totdiscrim1 = totdiscrim1 %>% as.character,
         countdiscrimwhy = factor(countdiscrimwhy,
                                  levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "no_discrim") %>% fct_rev))

colot_progress <- c("antiquewhite", "antiquewhite4", "antiquewhite3")
#' ### discrim1 is the sum of 5 types of discrimination H5MN5[A-E]
#' * H5MN5[A-E], each range of 0:3, never(0), rarely(1), sometimes(2), often(3)
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
B = ggplot(data_blood %>% filter(raceethnicity != "NA" & !is.na(discrim1)), aes(x = discrim1, fill = raceethnicity)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE) +
  scale_fill_manual(values = colot_progress, aesthetics = "fill")

B
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
  filter(raceethnicity=="NonHwhite") %>% 
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
  labs(title="white")

countdiscrim %>%
  filter(!is.na(discrim1)) %>% 
  filter(raceethnicity=="NonHblack") %>% 
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
  labs(title="black")

countdiscrim %>%
  filter(!is.na(discrim1)) %>% 
  filter(raceethnicity=="Hispanic") %>% 
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
  labs(title="hispanic")



#' ### discrim2 equals discrim1 if they feel the discrimination were due to race and ethnicity, otherwise 0
#' *  discrimination due to race and ethnicity: H5MN6A==1 (ancestry or national origin) or H5MN6D==1(race)
#' *  discrim2 ==0 : either they do not feel any discrimination or they do not think the discrimination were due to race and ethnicity

#+ echo=F, eval=T, warning=FALSE, message=FALSE 
B = ggplot(data_blood %>% filter(raceethnicity != "NA" & !is.na(discrim2)), aes(x = discrim2, fill = raceethnicity)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE) +
  scale_fill_manual(values = colot_progress, aesthetics = "fill")

B
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
  filter(raceethnicity=="NonHwhite") %>% 
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
  labs(title="white")

countdiscrim %>%
  filter(!is.na(discrim2)) %>% 
  filter(raceethnicity=="NonHblack") %>% 
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
  labs(title="black")

countdiscrim %>%
  filter(!is.na(discrim2)) %>% 
  filter(raceethnicity=="Hispanic") %>% 
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
  labs(title="hispanic")




#' ###  totdiscrim1 is the sum of H5MN5[A-E] and H5MN7(binary)
#' * here H5MN5[A-E] are binarized: 0 for never and rarely, 1 for sometimes and often
#+ echo=F, eval=T, warning=FALSE, message=FALSE 
B = ggplot(data_blood %>% filter(raceethnicity != "NA" & !is.na(totdiscrim1)), aes(x = totdiscrim1, fill = raceethnicity)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE) +
  scale_fill_manual(values = colot_progress, aesthetics = "fill")

B


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

countdiscrim %>%
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

countdiscrim %>%
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



#' ### totdiscrim2 equals totdiscrim1 if they feel the discrimination were due to race and ethnicity, otherwise 0
#' *  discrimination due to race and ethnicity: H5MN6A==1 (ancestry or national origin) or H5MN6D==1(race)
#' *  discrim2 ==0 : either they do not feel any discrimination or they do not think the discrimination were due to race and ethnicity

#+ echo=F, eval=T, warning=FALSE, message=FALSE 

B = ggplot(data_blood %>% filter(raceethnicity != "NA" & !is.na(totdiscrim2)), aes(x = totdiscrim2, fill = raceethnicity)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE) +
  scale_fill_manual(values = colot_progress, aesthetics = "fill")

B


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
  filter(raceethnicity=="NonHwhite") %>% 
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
  labs(title="white")

countdiscrim %>%
  filter(!is.na(totdiscrim2)) %>% 
  filter(raceethnicity=="NonHblack") %>% 
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
  labs(title="black")

countdiscrim %>%
  filter(!is.na(totdiscrim2)) %>% 
  filter(raceethnicity=="Hispanic") %>% 
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
  labs(title="hispanic")


#' ###  countdiscrimwhy is created: if totdiscrim1>0, sum "H5MN6[A-L]" (reasons of discrimination)
#+ echo=F, eval=T, warning=FALSE, message=FALSE 


B = ggplot(countdiscrim  %>% filter(raceethnicity != "NA" & !is.na(countdiscrimwhy)), aes(x = countdiscrimwhy, fill = raceethnicity)) +
  geom_bar(position = "stack", stat = "count", na.rm = TRUE) +
  scale_fill_manual(values = colot_progress, aesthetics = "fill")
  

B

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

countdiscrim %>%
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

countdiscrim %>%
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


