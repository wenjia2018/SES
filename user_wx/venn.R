#' ---
#' title: DE results for ses
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
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase)
library(ggformula)
library(here)
library(dbr)
library(venn)
walk(dir(path = here("R"),full.names = TRUE), source)
############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################

load_data(reconciled = FALSE, remove_inflam = FALSE)
define_treatments_and_controls()
recode_variables_in_dat()
# print(abbreviations)
funcs = str_subset(abbreviations$shorthand, "^m")
funcs = funcs %>% str_subset("m[6-8]")
# explicitly assign ncomp as the smallest number of table signatures gene numbers

ncomp = signatures$outcome_set[table1]%>% map_dfc(length) %>% unlist() %>%  min
if(0){
  example1 =
    args %>%
    filter(treatment %in%c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5", "sss_5",
                           # "ses_composite_pp1", "edu_p", "SEI_max_p_w12", "income_pp1_log"
    ),
    gene_set_name == "whole_genome_and_tfbm",
    names(controls) == "all") %>% slice_(2) %>% 
    mutate(out = pmap(., safely(model_fit), funcs),
           controls = names(controls))
}


#' ## venn without removing disease signatures
#+ echo=F, eval=T, warning=FALSE, message=FALSE

example1 = readRDS("/home/xu/ses-1/user_wx/DE_with1k.rds")


# ERRORS?
error = example1 %>%
  hoist(out, "error") %>%
  mutate(error = map(error, as.character)) %>%
  unnest(error)

# RESULTS
temp = example1 %>%
  hoist(out, "result")


result = temp %>%
  pluck("result") %>%
  set_names(temp$treatment)

DE = result %>% 
  map(pluck("ttT")) %>% 
  map(~ filter(.x,adj.P.Val<=0.05) %>% pull(gene))


names(DE) = c("Education", "Income", "Occupation", "SES Composite", "Subjective Social Status")
names(DE) = c("E", "I", "O", "S4", "Sss")


#' ## after removing all signatures in table 1

#+ echo=F, eval=T, warning=FALSE, message=FALSE

sig = Reduce(union, signatures$outcome_set[table1])

DE_removesig = map(DE, ~ setdiff(.x, sig))
names(DE_removesig) = c("Education", "Income", "Occupation", "SES Composite", "Subjective Social Status")
names(DE_removesig) = c("E", "I", "O", "S4", "Sss")
par(mfrow=c(1,2))

par(mar = rep(.5,4))
venn::venn(DE, ilabels = FALSE, zcolor = "style",snames = " ",
           ellipse = FALSE, opacity = 0.15, ilcs = 1.2,
           box = FALSE,
           sncs = 0.3)

venn::venn(DE_removesig, ilabels = TRUE, zcolor = "style",
           ellipse = FALSE, opacity = 0.15, ilcs = 1.2,
           box = FALSE,
           sncs = 0.3)

legend("topleft", legend=c("A", "B"),
       col=c("red", "blue"), cex=0.8)
mtext("Figure 2: Venn diagram",
      side = 1, line = -3, 
      outer = TRUE, cex = 1)
mtext("Left(whole genome), Right(removing disease signature)",
      side=1, line=-2, cex=1, outer=TRUE)  


vp1=venn.diagram(DE, filename  = NULL,
             fill =1:5,
             alpha = 0.1,
             main.cex=2,
             main.just = c(1,1),
             cex = 2, cat.cex= 1.7, cat.col=c(1:5), margin=0.3,
             cat.just=list(c(0.6,-0.1) , c(0.7,-0.5) , c(1,0.6) , c(0.2,1.2) , c(0.6 ,-2.2))
             )
# cat.just up - down +, left +, righ -
# margin increase +

grid.draw(vp1)

vp2=venn.diagram(DE_removesig, filename  = NULL,
                fill =1:5,
                alpha = 0.1,
                main.cex=2,
                main.just = c(1,1),
                cex = 2, cat.cex=1.7, cat.col=c(1:5), margin=0.3,
                cat.just=list(c(0.6,-0.1) , c(0.7,-0.5) , c(1,0.6) , c(0.2,1.2) , c(0.6 ,-2.2))
)
# cat.just up - down +, left +, righ -
# margin increase +

grid.draw(vp2)

grid.arrange(gTree(children=vp1), gTree(children=vp2), ncol=2)

#'`rmarkdown::render("/home/xu/ses-1/user_wx/venn.R")`
