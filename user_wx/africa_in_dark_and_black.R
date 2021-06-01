#' ---
#' title: check proportions of Africa origin in black and dark skin color
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

#' #### procedure of creating africa origin:
#' 
#' * H3OD8 == 2 (candidates have chosen more than 2 origins but they chose Africa as the main origin) 
#' * H3OD8 == 997 & H3OD7[A-D] (candidates chose Africa origin only)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
waves <- readRDS("/home/share/preprocessed_two_batches/waves_07.12.2020.rds")
dat <- readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_17.11.2020.rds")
AID_blood = dat@phenoData@data$AID
waves_full <- readRDS("/home/share/preprocessed_two_batches/waves_full.rds")
# waves_full$H3OD8 %>% table
# waves_full$H3OD7A %>% table
# 
# a = waves_full %>%
#   select(H3OD8,  matches("H3OD7[A-D]")) %>%
#   filter(H3OD7A ==2 | H3OD7B == 2 | H3OD7C == 2 | H3OD7D == 2)
waves_full = waves_full %>% 
  mutate_at(.vars = vars(matches("H3OD7[A-D]")), .funs = list(africaonly = ~ . %in% c(2,0,997))) %>% 
  mutate(africaonly = ifelse(H3OD7A_africaonly==TRUE & H3OD7B_africaonly==TRUE & H3OD7C_africaonly==TRUE & H3OD7A_africaonly==TRUE, TRUE, FALSE),
         africa = case_when(H3OD8==2 ~ "Africa",
                            H3OD8==997 & (africaonly==TRUE)  ~ "Africa",
                            TRUE ~ "Non Africa"))

# b = data %>% select(H3OD8,  matches("H3OD7[A-D]"), africaonly, africa)

data = waves %>% 
  filter(AID %in% AID_blood) %>% 
  left_join(waves_full %>% select(AID, africa)) %>% 
  mutate( color_byinterviewer3 = H3IR17 %>%
            as.character() %>% 
            as.factor %>% 
            fct_collapse(
              DarkBlack = c("1", "2"),
              LightMed = c("3", "4"),
              White = "5") %>%
            relevel(ref = "White"),
          color_byinterviewer5 = H3IR17 %>%
            as.character() %>% 
            as.factor %>% 
            fct_recode(
              White = "5",
              Light = "4",
              Medium = "3",
              Dark = "2",
              Black = "1"
            ) %>% 
            relevel(ref = "White"))

#' cross table of being africa and skin color of 5 levles
#+ echo=F, eval=T, warning=FALSE, message=FALSE
descr::crosstab(data$africa, data$color_byinterviewer5, prop.r = T, prop.c = T, prop.t = FALSE, plot = F)

data2 = data %>% 
  mutate(color_darkblack = color_byinterviewer5 %>% 
           fct_collapse(
             Dark = "Dark",
             Black = "Black",
             NULL = c("White", "Light", "Medium")))
#' cross table of being africa and skin color of dark and black 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
descr::crosstab(data2$africa, data2$color_darkblack, prop.r = T, prop.c = T, prop.t = FALSE, plot = T)

#' chi square test shows there is no evidence of association of being Africa and skin color (being black or dark)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
chisq.test(data2$africa, data2$color_darkblack)
