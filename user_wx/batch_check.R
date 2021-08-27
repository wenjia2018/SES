
#' ---
#' title: batch check
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE, include=FALSE, comment=NA

#' ### 3 batches
#' *  with plate information
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA

dat_batch1 <- readRDS("/home/share/preprocessing/preprocessed/dat_ref.rds")
dat_batch12 <- readRDS("/home/share/preprocessing/preprocessed_two_batches/dt_batches1_2_steve_waves_17.11.2020.rds")
dat_batch123 <- readRDS("/home/share/preprocessing/preprocessed_two_batches/all.batches.expression.set.rle_waves_05.08.2021.rds")


#' *  batch1
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA

dat_batch1$Plate %>% table()

#' *  batch1 and batch2
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
dat_batch12$Plate %>% table()

#' *  batch1, batch2 and batch3
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
dat_batch123$Plate %>% table()
