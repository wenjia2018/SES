library(tidyverse)
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase)
library(ggformula)
library(ggplot2)
library(corrplot)
library("ellipse")
library("RColorBrewer")

setwd("/home/share/projects/Mike/ses2020.7.17/ses/user_ms")
dat = readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_25.06.2020.rds")


phenSES <-pData(dat) %>% dplyr::select(edu_p,
                                income_pp1_log,
                                SEI_max_p_w12,
                                ses_composite_pp1,
                                edu_max,
                                income_hh_ff5,
                                SEI_ff5,
                                sss_5,
                                ses_sss_composite, 
                                ses_composite_ff5) %>%
                        mutate(edu_p = case_when(edu_p=="high" ~ 1,
                                                edu_p=="votec" ~ 2,
                                                edu_p=="college" ~ 3,
                                                edu_p=="post" ~ 4),
                               edu_max = case_when(edu_max=="high" ~ 1,
                                                   edu_max=="votec" ~ 2,
                                                   edu_max=="college" ~ 3,
                                                   edu_max=="post" ~ 4)) %>% 
  rename("Parental Education" = edu_p,
         "Parental Income"    = income_pp1_log,
         "Parental SEI"       = SEI_max_p_w12,
         "Parental SES Composite" =  ses_composite_pp1,
         "SES Composite 3" = ses_composite_ff5,
         "SES Composite" =  ses_sss_composite,
         "Education" = edu_max,
         "Income"    = income_hh_ff5,
         "SEI"       = SEI_ff5,
         "SSS" = sss_5
         )

phenSES <- phenSES %>% dplyr::select("Parental Education", "Parental Income"  ,
                        "Parental SEI", "Parental SES Composite" ,
                        "SES Composite 3" ,
                        "SES Composite" ,
                        "Education" ,
                        "Income"    ,
                        "SEI"       ,
                        "SSS")  

M<-cor(phenSES, use = "pairwise.complete.obs")

corrplot.mixed(M, lower = "number", upper = "circle", tl.pos = c("d", "lt", "n"), diag = c("n", "l", "u"), bg = "white", addgrid.col = "grey",
               lower.col = NULL, upper.col = NULL, plotCI = c("n", "square", "circle", "rect"), mar = c(0, 0, 0, 0))

corrplot::corrplot.mixed(M, lower="number", upper="circle", order = "original", sig.level = 0.05, insig="n", tl.pos="lt")

corrplot(M, method = "ellipse")