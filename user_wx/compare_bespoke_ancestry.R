
#' ---
#' title: compare different ancestry control and correction methods
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
walk(dir(path = here("R"),full.names = TRUE), source)
source("/home/xu/ses-1/user_wx/extract_v2.R")

source("/home/xu/ses-1/user_wx/plotting_utils.R")

example_skincolor3 = readRDS("~/ses-1/user_wx/skincolor_eqtl005_aging_composite_ancestry_11.05.2021.rds")
example_skincolor3_bespoke = readRDS("~/ses-1/user_wx/color3_bespoke_15.03.2021.rds")

p_eqtl = 0.05

#' ### omnibus fdr correction with different ancestry controls
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=10, fig.height=8
a = omnibus_plot(example_skincolor3, p_eqtl)

b = omnibus_plot(example_skincolor3_bespoke, p_eqtl)


# ggarrange(a, b+ rremove("y.text") + rremove("ylab")+ rremove("y.ticks"), ncol = 2, labels = c("A: bespoke ancestry aging_composite", "B: bespoke ancestry each signature"), common.legend = TRUE, legend="right")

ggarrange(a, b, ncol = 2, labels = c("A:aging_composite ancestry", "B: bespoke ancestry each signature"),
          vjust = 1,font.label = list(size = 10), common.legend = TRUE, legend="right")

#' ### aging composite ancestry control pca with different corrections
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10
c = pca_plotting(example_skincolor3, p_eqtl, "fdr")
d = pca_plotting(example_skincolor3, p_eqtl, "bonferroni")

e = ggarrange(c, d, ncol = 2, labels = c("A: fdr correction", "B: bonferroni correction"), font.label = list(size = 10), common.legend = TRUE, legend="right")
title <- expression(atop(bold("Figure:"), scriptstyle(" aging_composite ancestry")))
annotate_figure(e, bottom=text_grob(title))

#' ### bespoke ancestry control pca with different corrections
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10
f = pca_plotting(example_skincolor3_bespoke, p_eqtl, "fdr")
g = pca_plotting(example_skincolor3_bespoke, p_eqtl, "bonferroni")

h = ggarrange(f, g, ncol = 2, labels = c("A: fdr correction", "B: bonferroni correction"), font.label = list(size = 10), common.legend = TRUE, legend="right")
title <- expression(atop(bold("Figure:"), scriptstyle("each signature ancestry ")))
annotate_figure(h, bottom=text_grob(title))

#' ### aging composite ancestry control pca and signature bespoke ancestry control pca with fdr correction
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

i = ggarrange(c, f, ncol = 2, labels = c("A: aging composite ancestry", "B: bespoke ancestry"), font.label = list(size = 10), common.legend = TRUE, legend="right")
title <- expression(atop(bold("Figure:"), scriptstyle("fdr correction")))
annotate_figure(i, bottom=text_grob(title))


#' ### aging composite ancestry control pca and signature bespoke ancestry control pca with bonferroni correction
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

j = ggarrange(d, g, ncol = 2, labels = c("A: aging composite ancestry", "B: bespoke ancestry"), font.label = list(size = 10), common.legend = TRUE, legend="right")
title <- expression(atop(bold("Figure:"), scriptstyle("bonferroni correction")))
annotate_figure(j, bottom=text_grob(title))



#' ## results from race stratum
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

example_skincolor3_black = readRDS("~/ses-1/user_wx/skincolor3_NonHblack_strata_eqtl005_aging_composite_ancestry_12.05.2021.rds")
example_skincolor3_white = readRDS("~/ses-1/user_wx/skincolor3_NonHwhite_strata_eqtl005_aging_composite_ancestry_12.05.2021.rds")
example_skincolor3_hispanic = readRDS("~/ses-1/user_wx/skincolor3_Hispanic_strata_eqtl005_aging_composite_ancestry_12.05.2021.rds")


example_skincolor3_bespoke_black = readRDS("/home/xu/ses-1/user_wx/color3_bespoke_NonHblack_strata_16.03.2021.rds")
example_skincolor3_bespoke_white =  readRDS("/home/xu/ses-1/user_wx/color3_bespoke_NonHwhite_strata_16.03.2021.rds")
example_skincolor3_bespoke_hispanic = readRDS("/home/xu/ses-1/user_wx/color3_bespoke_Hispanic_strata_16.03.2021.rds")

#' ### omnibus fdr correction with different ancestry controls black race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

sba = omnibus_plot(example_skincolor3_black, p_eqtl)

sbb = omnibus_plot(example_skincolor3_bespoke_black, p_eqtl)

ggarrange(sba, sbb, ncol = 2, labels = c("A:aging_composite ancestry", "B: bespoke ancestry each signature"),
          vjust = 1,font.label = list(size = 10), common.legend = TRUE, legend="right")


#' ### omnibus fdr correction with different ancestry controls white race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

swa = omnibus_plot(example_skincolor3_white, p_eqtl)

swb = omnibus_plot(example_skincolor3_bespoke_white, p_eqtl)

ggarrange(swa, swb, ncol = 2, labels = c("A:aging_composite ancestry", "B: bespoke ancestry each signature"),
          vjust = 1,font.label = list(size = 10), common.legend = TRUE, legend="right")



#' ### omnibus fdr correction with different ancestry controls Hispanic race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

sha = omnibus_plot(example_skincolor3_hispanic, p_eqtl)

shb = omnibus_plot(example_skincolor3_bespoke_hispanic, p_eqtl)

ggarrange(sha, shb, ncol = 2, labels = c("A:aging_composite ancestry", "B: bespoke ancestry each signature"),
          vjust = 1,font.label = list(size = 10), common.legend = TRUE, legend="right")


#' ### aging composite ancestry control pca with different corrections in black race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

sbc = pca_plotting(example_skincolor3_black, p_eqtl, "fdr")
sbd = pca_plotting(example_skincolor3_black, p_eqtl, "bonferroni")

sbe = ggarrange(sbc, sbd, ncol = 2, labels = c("A: fdr correction", "B: bonferroni correction"), font.label = list(size = 10), common.legend = TRUE, legend="right")
title <- expression(atop(bold("Figure:"), scriptstyle("aging composite ancestry")))
annotate_figure(sbe, bottom=text_grob(title))


#' ### bespoke ancestry control pca with different corrections in black race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

sbf = pca_plotting(example_skincolor3_bespoke_black, p_eqtl, "fdr")
sbg = pca_plotting(example_skincolor3_bespoke_black, p_eqtl, "bonferroni")

sbh = ggarrange(sbf, sbg, ncol = 2, labels = c("A: fdr correction", "B: bonferroni correction"), font.label = list(size = 10), common.legend = TRUE, legend="right")
title <- expression(atop(bold("Figure:"), scriptstyle("each signature ancestry")))
annotate_figure(sbh, bottom=text_grob(title))


#' ### aging composite ancestry control pca and signature bespoke ancestry control pca with fdr correction in black race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

sbi = ggarrange(sbc, sbf, ncol = 2, labels = c("A: aging composite ancestry", "B: bespoke ancestry"), font.label = list(size = 10), common.legend = TRUE, legend="right")
title <- expression(atop(bold("Figure:"), scriptstyle("fdr correction")))
annotate_figure(sbi, bottom=text_grob(title))


#' ### aging composite ancestry control pca and signature bespoke ancestry control pca with bonferroni correction in black race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

sbj = ggarrange(sbd, sbg, ncol = 2, labels = c("A: aging composite ancestry", "B: bespoke ancestry"), font.label = list(size = 10), common.legend = TRUE, legend="right")
title <- expression(atop(bold("Figure:"), scriptstyle("bonferroni correction")))
annotate_figure(sbj, bottom=text_grob(title))




#' ### aging composite ancestry control pca with different corrections in white race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

swc = pca_plotting(example_skincolor3_white, p_eqtl, "fdr")
swd = pca_plotting(example_skincolor3_white, p_eqtl, "bonferroni")

swe = ggarrange(swc, swd, ncol = 2, labels = c("A: fdr correction", "B: bonferroni correction"), font.label = list(size = 10), common.legend = TRUE, legend="right")
title <- expression(atop(bold("Figure:"), scriptstyle(" aging_composite ancestrye")))
annotate_figure(swe, bottom=text_grob(title))


#' ### bespoke ancestry control pca with different corrections in white race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

swf = pca_plotting(example_skincolor3_bespoke_white, p_eqtl, "fdr")
swg = pca_plotting(example_skincolor3_bespoke_white, p_eqtl, "bonferroni")

swh = ggarrange(swf, swg, ncol = 2, labels = c("A: fdr correction", "B: bonferroni correction"), font.label = list(size = 10), common.legend = TRUE, legend="right")
title <- expression(atop(bold("Figure:"), scriptstyle("each signature ancestry")))
annotate_figure(swh, bottom=text_grob(title))


#' ### aging composite ancestry control pca and signature bespoke ancestry control pca with fdr correction in white race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

swi = ggarrange(swc, swf, ncol = 2, labels = c("A: aging composite ancestry", "B: bespoke ancestry"), font.label = list(size = 10), common.legend = TRUE, legend="right")
title <- expression(atop(bold("Figure:"), scriptstyle("fdr correction")))
annotate_figure(swi, bottom=text_grob(title))


#' ### aging composite ancestry control pca and signature bespoke ancestry control pca with bonferroni correction in white race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

swj = ggarrange(swd, swg, ncol = 2, labels = c("A: aging composite ancestry", "B: bespoke ancestry"), font.label = list(size = 10), common.legend = TRUE, legend="right")
title <- expression(atop(bold("Figure:"), scriptstyle("bonferroni correction")))
annotate_figure(swj, bottom=text_grob(title))




#' ### aging composite ancestry control pca with different corrections in hispanic race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

shc = pca_plotting(example_skincolor3_hispanic, p_eqtl, "fdr")
shd = pca_plotting(example_skincolor3_hispanic, p_eqtl, "bonferroni")

she = ggarrange(shc, shd, ncol = 2, labels = c("A: fdr correction", "B: bonferroni correction"), font.label = list(size = 10), common.legend = TRUE, legend="right")
title <- expression(atop(bold("Figure:"), scriptstyle(" aging_composite ancestrye")))
annotate_figure(she, bottom=text_grob(title))


#' ### bespoke ancestry control pca with different corrections in hispanic race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

shf = pca_plotting(example_skincolor3_bespoke_hispanic, p_eqtl, "fdr")
shg = pca_plotting(example_skincolor3_bespoke_hispanic, p_eqtl, "bonferroni")

shh = ggarrange(shf, shg, ncol = 2, labels = c("A: fdr correction", "B: bonferroni correction"), font.label = list(size = 10), common.legend = TRUE, legend="right")
title <- expression(atop(bold("Figure:"), scriptstyle("each signature ancestry")))
annotate_figure(shh, bottom=text_grob(title))


#' ### aging composite ancestry control pca and signature bespoke ancestry control pca with fdr correction in hispanic race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

shi = ggarrange(shc, shf, ncol = 2, labels = c("A: aging composite ancestry", "B: bespoke ancestry"), font.label = list(size = 10), common.legend = TRUE, legend="right")
title <- expression(atop(bold("Figure:"), scriptstyle("fdr correction")))
annotate_figure(shi, bottom=text_grob(title))


#' ### aging composite ancestry control pca and signature bespoke ancestry control pca with bonferroni correction in hispanic race straturm
#+ echo=F, eval=T, warning=FALSE, message=FALSE, fig.width=12, fig.height=10

shj = ggarrange(shd, shg, ncol = 2, labels = c("A: aging composite ancestry", "B: bespoke ancestry"), font.label = list(size = 10), common.legend = TRUE, legend="right")
title <- expression(atop(bold("Figure:"), scriptstyle("bonferroni correction")))
annotate_figure(shj, bottom=text_grob(title))

