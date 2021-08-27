
#' ---
#' title: model comparision metrics 
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE



#' ## 1. compare two models with transformed outcome (same oucome but with log transform)
#' * model comparison metrics is different even with the same data after transformation
#+ echo=T, eval=T, warning=FALSE, message=FALSE
library(tidyverse)
library(BayesFactor)
fit = lm(carb ~ . , data = mtcars)
fit_log = lm(log(carb) ~ ., data = mtcars)
BIC(fit)
BIC(fit_log)



# set.seed(12313231)
# x = rnorm(1000, 2, 3)
# m = 1.5 * x + rnorm(1000)
# y = 4 * x + rnorm(1000)
# 
# fit1 = lm(y ~ x + m)
# BIC(fit1)
# 
# fit1 %>% summary
# 
# fit2 = lm(m ~ x + y)
# BIC(fit2)
# 
# fit2 %>% summary


#' ## 2. compare two models, the true data generating process is x -> m -> y
#' * model comparison metrics is favoring the alternative data generating process over the true data generating process
#+ echo=T, eval=T, warning=FALSE, message=FALSE

set.seed(12313222)
x = rnorm(1000, 10, 3)
m = 4 + 1.5 * x + rnorm(1000)
y = 2.5 + 3 * x + 0.6 * m + rnorm(1000)

dt = tibble(x, m , y)

fit1 = lm(y ~ x + m)
mod1 = lmBF(y ~ x + m, data = dt)
# linearReg.R2stat(N=1000, p=2, R2=0.9928, rscale = 0.353553390593274)[['bf']]
fit1 %>% summary

#' ### the true data generating process is x -> m -> y
#+ echo=T, eval=T, warning=FALSE, message=FALSE
BIC(fit1)
AIC(fit1)


fit2 = lm(m ~ x + y)
mod2 = lmBF(m ~ x + y, data = dt)
fit2 %>% summary
# linearReg.R2stat(N=1000, p=2, R2=0.965, rscale = 0.353553390593274)[['bf']]

#' ### the alternative data generating process is x -> y ->  m
#+ echo=T, eval=T, warning=FALSE, message=FALSE
BIC(fit2)
AIC(fit2)

