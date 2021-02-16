
#' ---
#' title: discrimation variable distribution fit
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE


#' ## totdiscrim1_gamma, totdiscrim2_gamma, countdiscrimwhy_gamma
#+ echo=F, eval=T, warning=FALSE, message=FALSE
set.seed(123)
library(here)
library(tidyverse)
library(fitdistrplus)

walk(dir(path = here("R"),full.names = TRUE), source)
fit_m4 = partial(fit_m4, n_perm = 1000) # specify n_perm
load_data(reconciled = FALSE, remove_inflam = FALSE)
df<- pData(dat) %>%
  mutate(
    totdiscrim1_category = case_when(
      totdiscrim1 %in% c(0, 1, 2, 3) ~ totdiscrim1,
      totdiscrim1 %in% c(4, 5, 6) ~ 4),
    countdiscrimwhy = ifelse(countdiscrimwhy=="no_discrim", NA, countdiscrimwhy) %>% as.numeric(),
    totdiscrim2_binary = ifelse(totdiscrim2 > 0, 1, 0),
    discrim2_binary = ifelse(discrim2 > 0, 1, 0),
    # add 1 for all in oder to fit a gamma distribution in glm (exponential has to be fitted by gamma first then
    # set the dispersion=1)
    totdiscrim1_gamma = totdiscrim1 + 0.001, #exponential fit glm for mediation, add 0.001 to avoid 0
    totdiscrim2_gamma = totdiscrim2 + 0.001, #exponential fit glm for mediation, add 0.001 to avoid 0
    countdiscrimwhy_gamma = countdiscrimwhy + 0.001#exponential fit glm for mediation, add 0.001 to avoid 0
    
  ) %>% 
  dplyr::select(totdiscrim1_gamma, totdiscrim2_gamma, countdiscrimwhy_gamma)





#' ### totdiscrim1 
#+ echo=F, eval=T, warning=FALSE, message=FALSE

y = df$totdiscrim1_gamma
keep = y %>% complete.cases()
y = y[keep]

fit.gamma <- fitdist(y, distr = "gamma", method = "mle")
summary(fit.gamma)
plot(fit.gamma)

# 
# fit.pois <- fitdist(y-0.001, distr = "pois", method = "mle")
# summary(fit.pois)
# plot(fit.pois)
# 
# fit.geom <- fitdist(y-0.001, distr = "geom", method = "mle")
# summary(fit.geom)
# plot(fit.geom)


#' ### totdiscrim2 
#+ echo=F, eval=T, warning=FALSE, message=FALSE

y = df$totdiscrim2_gamma
keep = y %>% complete.cases()
y = y[keep]

fit.gamma <- fitdist(y, distr = "gamma", method = "mle")
summary(fit.gamma)
plot(fit.gamma)

# fit.pois <- fitdist(y-0.001, distr = "pois", method = "mle")
# summary(fit.pois)
# plot(fit.pois)
# 
# fit.geom <- fitdist(y-0.001, distr = "geom", method = "mle")
# summary(fit.geom)
# plot(fit.geom)

#' ### countdiscrimwhy 
#+ echo=F, eval=T, warning=FALSE, message=FALSE

y = df$countdiscrimwhy_gamma
keep = y %>% complete.cases()
y = y[keep]

fit.gamma <- fitdist(y, distr = "gamma", method = "mle")
summary(fit.gamma)
plot(fit.gamma)

# fit.pois <- fitdist(y-0.001, distr = "pois", method = "mle")
# summary(fit.pois)
# plot(fit.pois)
# 
# fit.geom <- fitdist(y-0.001, distr = "geom", method = "mle")
# summary(fit.geom)
# plot(fit.geom)
