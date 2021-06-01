
#' ---
#' title: aging composite ancestry control with fdr correction
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
# skincolor_eqtl005_aging_composite_ancestry_11.05.2021.rds m7 m8 without genowide mediation
example_skincolor3 = readRDS("~/ses-1/user_wx/skincolor_eqtl005_aging_composite_ancestry_11.05.2021.rds")

p_eqtl = 0.05

#' ### 1 omnibus fdr correction with aging composite ancestry controls
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=10, fig.height=8
a = omnibus_plot(example_skincolor3, p_eqtl)
a



#' ### 2 aging composite ancestry control pca with fdr corrections
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10
c = pca_plotting(example_skincolor3, p_eqtl, "fdr")
c
#' ### figure1 two panels composite ancestry control with fdr corrections
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10
pdf("./user_wx/figure.pdf", width = 12, height = 7,onefile=FALSE)
ggarrange(a, c, ncol = 2, labels = c("A", "B"), common.legend = TRUE, legend="right")
dev.off()
#' ### 3 aging composite ancestry control pca mediation with fdr corrections (no significant mediation)
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

# d = p_eqtl %>% map_df(outm7med, mediators, control,example) 

# d %>% filter(p_med<0.05)

#' ## 4 results from race stratum
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

example_skincolor3_black = readRDS("~/ses-1/user_wx/skincolor3_NonHblack_strata_eqtl005_aging_composite_ancestry_12.05.2021.rds")
example_skincolor3_white = readRDS("~/ses-1/user_wx/skincolor3_NonHwhite_strata_eqtl005_aging_composite_ancestry_12.05.2021.rds")
example_skincolor3_hispanic = readRDS("~/ses-1/user_wx/skincolor3_Hispanic_strata_eqtl005_aging_composite_ancestry_12.05.2021.rds")

#' ### omnibus fdr correction with aging composite ancestry controls black race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

sba = omnibus_plot(example_skincolor3_black, p_eqtl)

sba


#' ### omnibus fdr correction with aging composite ancestry controls white race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

swa = omnibus_plot(example_skincolor3_white, p_eqtl)

swa

#' ### omnibus fdr correction with aging composite ancestry controls Hispanic race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

sha = omnibus_plot(example_skincolor3_hispanic, p_eqtl)

sha

#' ### aging composite ancestry control pca with fdr correction in black race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

sbc = pca_plotting(example_skincolor3_black, p_eqtl, "fdr")
sbc


#' ### aging composite ancestry control pca with fdr corrections in white race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

swc = pca_plotting(example_skincolor3_white, p_eqtl, "fdr")
swc


#' ### aging composite ancestry control pca with fdr corrections in hispanic race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

shc = pca_plotting(example_skincolor3_hispanic, p_eqtl, "fdr")
shc

ggarrange(sbc, swc, ncol = 2, labels = c("A", "B"), common.legend = TRUE, legend="right")

