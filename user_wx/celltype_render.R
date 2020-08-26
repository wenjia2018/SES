
#' ---
#' title: celltype composition
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE


library(tidyverse)
library(rlang)
library(compositions)
# library(easyCODA)
example3 = readRDS("/home/share/scratch/xu/example3_celltype.rds")
b=example3$out[[9]]$result$b*10

barplot(b, args.legend = list(x = "topright", cex = 0.7),  main = "Cell type distribution (proportion)")

barplot(b %>% sort(),
        las=2,
                cex.names=0.4,
                # col = heat.colors(21),
                # col = rainbow(21),
                # col =RColorBrewer::brewer.pal(1:length(b),"Purples"),
                beside = TRUE,
                horiz = TRUE,
                legend.text = NULL,
        main = "Cell type distribution (magnitude)"
        )
# barplot(b*10,  
#         names.arg = names(b),
#         las=2,
#         cex.names=0.4,
#         # col = heat.colors(21),
#         # col = rainbow(21),
#         # col =RColorBrewer::brewer.pal(1:length(b),"Purples"),
#         beside = TRUE,
#         horiz = TRUE,
#         legend.text = NULL,
#         # xaxt = "n",#remove x label
#         # args.legend = list(x = "topright", cex = 0.35),
#         )


#'  rmarkdown::render("/home/xu/ses-1/user_wx/celltype_render.R")
