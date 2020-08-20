#' ---
#' title: Tables comparing the different samples for the Supplement
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#' Set global options 
#+ setup, warning=FALSE, message=FALSE
# knitr::opts_chunk$set(echo = FALSE)


library(tidyverse)
library(fastDummies)
library(foreign)
library(lubridate)
library(readxl)
library(haven)
library(Biobase)

# library(SASxport)

########################################################
# AUTHORS, JUSTIN CHUMBLEY, CECILIA POTENTE
########################################################

########################################################
# LOAD DATA
########################################################
# MULTIPLE AUTHORS OF THIS CODE: NEEDS MORE STYLE.
server = TRUE
if(server){
  data_input = "/home/share/data_input/" 
  data_output = "/home/share/preprocessed_two_batches" 
}else{
  data_input = "/Volumes/Share/data_input/" 
  data_output = "/Volumes/Share/preprocessed_two_batches/"
}

wave5 <- read.xport(str_c(data_input, "wave5.xpt")) %>% as_tibble
w5AID <- wave5$AID


setwd("/home/share/projects/bmi_lifecourse")
waves  <- readRDS("/home/share/preprocessed_two_batches/waves_25.06.2020.rds")  # path relative to 
datao  <- readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_25.06.2020.rds")
dataophen <- pData(datao) %>% as_tibble()
aidbatch1 <- dataophen %>% dplyr::select(AID)
waves_b1 <- waves %>% right_join(aidbatch1, by="AID") %>% transmute(batch1=1, AID=AID)
waves <- waves %>% left_join(waves_b1, by="AID")
waves <- waves %>% dplyr:: mutate(wav5tot= case_when(wave5==T ~ 1),
                                  wave5total=case_when(batch1==1 & wav5tot==1 ~ 0,
                                                       is.na(batch1) & wav5tot==1 ~ 1),
                                  non_missing_parentalSES = case_when(batch1==1 & wav5tot==1 & is.na(ses_composite_pp1)  ~ 0,
                                                                      batch1==1 & wav5tot==1 & !is.na(ses_composite_pp1) ~ 1))

waves <- waves %>% filter(wav5tot==1)

library(table1)
library(MatchIt) 

rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- waves[[name]]
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- t.test(y ~ waves$wave5total)$p.value
    } else {
      p <- chisq.test(table(y, droplevels(waves$wave5total)))$p.value
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

label(waves$sex_interv)    <- "Sex"
label(waves$re_2)    <- "Black Non-Hispanic"
label(waves$re_3)    <- "Asian Non-Hispanic"
label(waves$re_4)    <- "Other Non-Hispanic"
label(waves$re_5)    <- "Hispanic"
label(waves$sex_interv_m)    <- "Male"
label(waves$BirthY)    <- "Birth Year"
label(waves$edu_max_votec)    <- "Vocational Education"
label(waves$edu_max_post)    <- "Higher than college"
label(waves$edu_max_college)    <- "College"
label(waves$ses_sss_composite)    <- "ses4 (edu,income,SEI,sss)"
label(waves$ses_composite_ff5)    <- "ses3 (edu,income,SEI)"
label(waves$SEI_ff5)    <- "social economic index SEI"
label(waves$sss_5)    <- "SSS"
label(waves$income_hh_ff5)    <- "Income"
label(waves$work_collar_ff5)    <- "Work Collar"
label(waves$work_collar_rf_f12)    <- "Work collar father"
label(waves$work_collar_rm_f12)    <- "Work collar mother"
label(waves$ses_composite_pp1)    <- "ses (edu,income,SEI)"
label(waves$edu_p)    <- "parent education"
label(waves$SEI_max_p_w12)    <- "Parent SEI"
label(waves$income_pp1_log)    <- "Parent Income"

waves$wave5total<- factor(waves$wave5total, levels=0:1, labels=c("Wave 5 Batch 1", "Wave 5 No Batch 1"))
table1::table1(~ sex_interv_m + re_2+re_3+ re_4+ re_5+BirthY+ses_sss_composite+ses_composite_ff5+sss_5+SEI_ff5+edu_max+income_hh_ff5+work_collar_ff5+work_collar_rf_f12+work_collar_rm_f12+ses_composite_pp1+edu_p+SEI_max_p_w12+income_pp1_log | wave5total, data = waves, overall=T)

table1::table1(~ sex_interv_m + re_2+re_3+ re_4+ re_5+BirthY+ses_sss_composite+ses_composite_ff5+sss_5+SEI_ff5+edu_max+income_hh_ff5+work_collar_ff5+work_collar_rf_f12+work_collar_rm_f12+ses_composite_pp1+edu_p+SEI_max_p_w12+income_pp1_log | wave5total, data = waves, droplevel=F, render=rndr, render.strat=rndr.strat, render.default=T, overall=F)

waves$non_missing_parentalSES     <- factor(waves$non_missing_parentalSES, levels=0:1, labels=c("Missing Parental SES", "No Missing Parental SES"))
table1::table1(~ sex_interv_m + re_2+re_3+ re_4+ re_5+BirthY+ses_sss_composite+ses_composite_ff5+sss_5+SEI_ff5+edu_max+income_hh_ff5+work_collar_ff5+work_collar_rf_f12+work_collar_rm_f12+ses_composite_pp1+edu_p+SEI_max_p_w12+income_pp1_log | non_missing_parentalSES, data = waves, droplevel=F, render=rndr, render.strat=rndr.strat, render.default=T, overall=F)