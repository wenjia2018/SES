#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(here)
library(tidyverse)
library(rlang)
library(skimr)
library(gtools)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase)
library(nnet)
library(broom)
walk(dir(path = here("R"),full.names = TRUE), source)
load_data(reconciled = FALSE, remove_inflam = FALSE)
dt_color_snp <- readRDS("/home/share/dna_ancestry/dna/dt_color_snp.rds")
snp <- pData(dat) %>%
  dplyr::select(AID, 
                re,
                H3IR17) %>% 
  left_join(dt_color_snp) %>% 
  mutate_at(.vars = vars(starts_with("rs")),
            .funs = list(~ ifelse(.=="00", NA, .)))

dt <- snp %>%  
  mutate(
    color_byinterviewer3 = H3IR17 %>%
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
      relevel(ref = "White"),
    raceethnicity = re %>%
                  fct_recode(
                    white = "1",
                    # white nonhispanic
                    black = "2",
                    # black nonhispanic
                    NULL = "3",
                    # asian nonhispanic
                    NULL = "4",
                    # other nonhispanic
                    hispanic = "5"
                  ) %>% relevel(ref = "white")
  ) %>% 
  select(rs1110400, raceethnicity, color_byinterviewer3)

#' rs1110400 by race
#+ echo=F, eval=T, warning=FALSE, message=FALSE


descr::crosstab(dt$rs1110400,  dt$raceethnicity, prop.r = T, prop.c = T, prop.t = FALSE, plot = F)




#' rs1110400 by skin color
#+ echo=F, eval=T, warning=FALSE, message=FALSE

  descr::crosstab(dt$rs1110400,  dt$color_byinterviewer3, prop.r = T, prop.c = T, prop.t = FALSE, plot = F)

  