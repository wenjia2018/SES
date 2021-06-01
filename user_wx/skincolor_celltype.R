
#' ---
#' title: celltype plotting for skincolor
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE


#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=8
library(tidyverse)
library(here)
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase)
library(enrichplot)
library(gridExtra)
library(grid)
library(ggpubr)
library(dbr) # my package
library(stringi)
walk(dir(path = here("R"),full.names = TRUE), source)
source("/home/xu/ses-1/user_wx/extract_v2.R")
load_data(reconciled = FALSE, remove_inflam = FALSE)
source("/home/xu/ses-1/user_wx/plotting_utils.R")
control = "ancestryPC_ses"

 example3 = readRDS("/home/xu/ses-1/user_wx/skincolor_celltype_eqtl005_aging_composite_ancestry_20.05.2021.rds")
 example3 = example3$out[[1]]$result$example0
 plot_celltype = function(example_data, focal) {
   b = example_data %>% filter(treatment %>% str_detect(focal)) %>% hoist(out,b = list("result","b")) %>% pull(b) %>% `[[`(1)
   N = length(b)
   b = b*10

   data = data.frame(celltype=names(b %>% sort()), value = b %>% sort(), row.names=NULL)
   data$celltype = data$celltype %>% 
     str_replace_all("\\.", " ") %>% 
     str_to_title() %>% 
     str_squish()
   
   data$celltype <- factor(data$celltype, # Factor levels in decreasing order
                           levels = data$celltype[order(data$value, decreasing = FALSE)])
   
   ggplot(data = data, aes(x = celltype, y = value, fill=celltype))+
     geom_bar(stat="identity")+
     scale_fill_grey()+
     ylim(0, 0.5)+
     geom_hline(yintercept = 1/N, linetype="dotted",color = "red") +
     coord_flip() +
     labs(caption = paste("Figure S3. Cell composition as a function of skin color", focal,"\n Add Health (n = 1655)"), face="bold", x = "Cell Types", y = "Magnitude") +
     theme(
       axis.ticks.y = element_blank(),
       legend.position = "none",
       panel.background = element_rect(fill = "transparent"), # bg of the panel
       plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
       panel.grid.major = element_blank(), # get rid of major grid
       panel.grid.minor = element_blank(), # get rid of minor grid
       legend.background = element_rect(fill = "transparent"), # get rid of legend bg
       legend.box.background = element_rect(fill = "transparent"), # get rid of
       # LABLES APPEARANCE
       plot.caption = element_text(hjust = 0.5, size=14, face= "bold", colour= "black", margin = margin(t = 20, unit = "pt")),
       axis.title.x = element_text(size=12, face="bold", colour = "black"),
       axis.title.y = element_text(size=12, face="bold",colour = "black"),
       axis.text.x = element_text(size=12, face="bold", colour = "black"), 
       axis.text.y = element_text(size=12, face="bold", colour = "black") 
     ) 
 }
 plot_celltype(example3, focal = "DarkBlack")
 plot_celltype(example3, focal = "LightMed")
 