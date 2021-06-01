#' ---
#' title: Supplementary Material
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE, include=FALSE
#+ echo=FALSE

set.seed(123)
library(here)
library(tidyverse)
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase)
library(enrichplot)
library(stringi)
library(dbr) # my package
#'## Descriptives

walk(dir(path = here("R"),full.names = TRUE), source)
fit_m4 = partial(fit_m4, n_perm = 1000) # specify n_perm

# WHICH EXAMPLES TO RUN? 
example4 <- example3 <- example2 <- example1 <- example0 <- FALSE
example4 <- TRUE
load_data(reconciled = FALSE, remove_inflam = TRUE)
des <- phenoData(dat)@data %>% as_tibble()

library(table1)
table1::label(des$ses_sss_composite) <- "SES Composite"
table1::label(des$income_hh_ff5) <- "Income"
table1::label(des$SEI_ff5) <- "Occupation"
table1::label(des$edu_max) <- "Education"
table1::label(des$sss_5) <- "Subjective Social Status"
table1::label(des$sex_interv) <- "Sex"
table1::label(des$re) <- "Race and Ethnicity"
table1::label(des$age_w5) <- "Age at Wave V"
table1::label(des$BirthY) <- "Year of birth"
table1::label(des$W5REGION) <- "Region"
table1::label(des$kit_biow5) <- "Kit"
table1::label(des$tube_biow5) <- "Tube"
table1::label(des$travel_biow5) <- "Travelling"
table1::label(des$stress_perceived) <- "Stress Perceived"
table1::label(des$bills) <- "Problems Paying Bills"
table1::label(des$currentsmoke) <- "Current Smoking"
table1::label(des$w5bmi) <- "BMI Wave V"
table1::label(des$insurance_lack) <- "Lack of Health Insurance"
table1::label(des$lowbirthweight) <- "Low Birthweight"


table1::table1(~ses_sss_composite+income_hh_ff5+SEI_ff5+sss_5+edu_max+sex_interv+re+age_w5+BirthY+
                 W5REGION+kit_biow5+tube_biow5+travel_biow5+stress_perceived+bills+currentsmoke+
                 w5bmi+insurance_lack+lowbirthweight | sex_interv, data = des)


example0 = readRDS("/home/share/projects/aging/example0_mediation_all_5_tot.rds")

ex2 <- example0 %>%
  hoist(out, p = list("result", "m8_fdr", 1, "other", "m")) %>% 
  filter(controls=="all") %>% 
  dplyr::select(treatment, gene_set_name, p) %>% 
  dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  mutate(treatment= case_when(treatment =="work_collar_ff5" ~ "Occupation Work Collar",
                              treatment =="edu_max" ~ "Education" ,
                              treatment =="income_hh_ff5" ~ "Income"     ,
                              treatment =="SEI_ff5" ~ "Occupation"      ,
                              treatment =="ses_sss_composite" ~ "SES Composite"  ,
                              treatment =="sss_5" ~ "Subjective Social Status",
                              treatment =="ses_composite_ff5"  ~ "SES Composite 3")) %>% 
  mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
         gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
         gene_set_name= case_when(gene_set_name =="Aging Up Cl4" ~ "Lysosome Metabolism" ,
                                  gene_set_name =="Aging Up Cl1" ~ "Innate/Adaptive Immunity",
                                  gene_set_name =="Aging Up Cl3" ~ "Fatty Acid Metabolism" ,
                                  gene_set_name =="Aging Up Cl2" ~ "Actine Regulation",
                                  gene_set_name =="Aging Down Cl3" ~ "Immune Genes",
                                  gene_set_name =="Aging Down Cl2" ~ "Ribosome"  ,
                                  gene_set_name =="Aging Down Cl1b" ~ "Mitochondrial",
                                  gene_set_name =="Aging Down Cl1a" ~ "RNA Metabolism",
                                  gene_set_name =="Aging Down Cl1c" ~ "DNA Replication",
                                  gene_set_name =="Aging" ~ "Aging All Genes"))


ex2hist <-ex2 %>% unnest(p) %>% dplyr:: filter(gene_set_name== c("Aging All Genes"))

#'## Histograms
#
#'### Aging All Genes Histogram

ggplot(ex2hist, aes(x = logFC, fill = as.factor(treatment))) +
  geom_density() +
  facet_wrap( ~ treatment, scales = "free")  +
  theme(legend.position = "none",
        panel.margin.x = unit(1, "lines"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color = "grey70"),
        axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
        axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
        plot.title = element_text(size = 20, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
        plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
        strip.text = element_text(size = 7),
        text = element_text(family = "Georgia"))+
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")


#'### Lysosome Metabolism Genes Histogram

ex2hist <-ex2 %>% unnest(p) %>% dplyr:: filter(gene_set_name== c("Lysosome Metabolism"))

ggplot(ex2hist, aes(x = logFC, fill = as.factor(treatment))) +
  geom_density() +
  facet_wrap( ~ treatment, scales = "free")  +
  theme(legend.position = "none",
        panel.margin.x = unit(1, "lines"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color = "grey70"),
        axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
        axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
        plot.title = element_text(size = 20, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
        plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
        strip.text = element_text(size = 7),
        text = element_text(family = "Georgia"))+
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")

#'### Innate/Adaptive Immunity Metabolism Genes Histogram

ex2hist <-ex2 %>% unnest(p) %>% dplyr:: filter(gene_set_name== c("Innate/Adaptive Immunity"))

ggplot(ex2hist, aes(x = logFC, fill = as.factor(treatment))) +
  geom_density() +
  facet_wrap( ~ treatment, scales = "free")  +
  theme(legend.position = "none",
        panel.margin.x = unit(1, "lines"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color = "grey70"),
        axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
        axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
        plot.title = element_text(size = 20, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
        plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
        strip.text = element_text(size = 7),
        text = element_text(family = "Georgia"))+
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")

#'### Fatty Acid Metabolism Genes Histogram

ex2hist <-ex2 %>% unnest(p) %>% dplyr:: filter(gene_set_name== c("Fatty Acid Metabolism"))

ggplot(ex2hist, aes(x = logFC, fill = as.factor(treatment))) +
  geom_density() +
  facet_wrap( ~ treatment, scales = "free")  +
  theme(legend.position = "none",
        panel.margin.x = unit(1, "lines"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color = "grey70"),
        axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
        axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
        plot.title = element_text(size = 20, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
        plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
        strip.text = element_text(size = 7),
        text = element_text(family = "Georgia"))+
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")

#'### Actine Regulation Metabolism Genes Histogram

ex2hist <-ex2 %>% unnest(p) %>% dplyr:: filter(gene_set_name== c("Actine Regulation"))

ggplot(ex2hist, aes(x = logFC, fill = as.factor(treatment))) +
  geom_density() +
  facet_wrap( ~ treatment, scales = "free")  +
  theme(legend.position = "none",
        panel.margin.x = unit(1, "lines"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color = "grey70"),
        axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
        axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
        plot.title = element_text(size = 20, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
        plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
        strip.text = element_text(size = 7),
        text = element_text(family = "Georgia"))+
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")

#'### Immune Genes Histogram

ex2hist <-ex2 %>% unnest(p) %>% dplyr:: filter(gene_set_name== c("Immune Genes"))

ggplot(ex2hist, aes(x = logFC, fill = as.factor(treatment))) +
  geom_density() +
  facet_wrap( ~ treatment, scales = "free")  +
  theme(legend.position = "none",
        panel.margin.x = unit(1, "lines"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color = "grey70"),
        axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
        axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
        plot.title = element_text(size = 20, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
        plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
        strip.text = element_text(size = 7),
        text = element_text(family = "Georgia"))+
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")

#'### Ribosome Genes Histogram

ex2hist <-ex2 %>% unnest(p) %>% dplyr:: filter(gene_set_name== c("Ribosome"))

ggplot(ex2hist, aes(x = logFC, fill = as.factor(treatment))) +
  geom_density() +
  facet_wrap( ~ treatment, scales = "free")  +
  theme(legend.position = "none",
        panel.margin.x = unit(1, "lines"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color = "grey70"),
        axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
        axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
        plot.title = element_text(size = 20, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
        plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
        strip.text = element_text(size = 7),
        text = element_text(family = "Georgia"))+
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")

#'### Mitochondrial Genes Histogram

ex2hist <-ex2 %>% unnest(p) %>% dplyr:: filter(gene_set_name== c("Mitochondrial"))

ggplot(ex2hist, aes(x = logFC, fill = as.factor(treatment))) +
  geom_density() +
  facet_wrap( ~ treatment, scales = "free")  +
  theme(legend.position = "none",
        panel.margin.x = unit(1, "lines"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color = "grey70"),
        axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
        axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
        plot.title = element_text(size = 20, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
        plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
        strip.text = element_text(size = 7),
        text = element_text(family = "Georgia"))+
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")

#'### RNA Metabolism Genes Histogram

ex2hist <-ex2 %>% unnest(p) %>% dplyr:: filter(gene_set_name== c("RNA Metabolism"))

ggplot(ex2hist, aes(x = logFC, fill = as.factor(treatment))) +
  geom_density() +
  facet_wrap( ~ treatment, scales = "free")  +
  theme(legend.position = "none",
        panel.margin.x = unit(1, "lines"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color = "grey70"),
        axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
        axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
        plot.title = element_text(size = 20, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
        plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
        strip.text = element_text(size = 7),
        text = element_text(family = "Georgia"))+
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")

#'### DNA Replication Genes Histogram

ex2hist <-ex2 %>% unnest(p) %>% dplyr:: filter(gene_set_name== c("DNA Replication"))

ggplot(ex2hist, aes(x = logFC, fill = as.factor(treatment))) +
  geom_density() +
  facet_wrap( ~ treatment, scales = "free")  +
  theme(legend.position = "none",
        panel.margin.x = unit(1, "lines"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color = "grey70"),
        axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
        axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
        plot.title = element_text(size = 20, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
        plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
        strip.text = element_text(size = 7),
        text = element_text(family = "Georgia"))+
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")

set.seed(123)
library(here)
library(tidyverse)
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase)
library(enrichplot)
library(dbr) # my package
library(rlist)
walk(dir(path = here("R"),full.names = TRUE), source)
load_data(reconciled = FALSE, remove_inflam = TRUE)

############################################################
# after obtained analysis results containing mm8_fdr
# pathway presenting
# please specify the example object and the signature you
# want to look at
############################################################

#in case you want to recompute all pathways (TRUE)
recomputeGSEA<-F

#in case you want to keep all pathways (TRUE), otherwise we keep just
#the pathways which have at least one aging gene among their leading genes
keepallpathways<-F

# function gsea_webgestalt calculate functional pathway using preferred database:
# pathway_KEGG and pathway_Reactome for example
# you could change the argument in the function 
gsea_webgestalt = function(treatment, ttT, file_output, enrichMethod = "GSEA", enrichDatabase="pathway_Reactome"){
  rankFile = str_c(file_output,"/", treatment,".rnk")
  
  ttT %>%
    dplyr::select(gene, logFC) %>% 
    write.table(file = rankFile, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  
  enrichResult = WebGestaltR::WebGestaltR(interestGeneFile = rankFile,
                                          interestGeneType = "genesymbol",
                                          enrichMethod = enrichMethod,
                                          organism = "hsapiens",
                                          enrichDatabase = enrichDatabase,
                                          fdrThr = 1,
                                          minNum=5,
                                          perNum = 1000,
                                          outputDirectory = file_output)
}

if(recomputeGSEA==TRUE){ 
  #gene_signature = "aging_mRNA"
  # create a folder in your current directory to save all the intermediate files and report generated at last from
  # webgestalt, you could name your new folder differently
  tempfolder = "temp_webgestalt"
  dir.create(tempfolder)
  example =readRDS(file="/home/share/projects/aging/tt_whole_genome_5.rds")
  
  args_whole_genome = 
    example %>% 
    hoist(out, ttT = list("result", "ttT")) %>% 
    dplyr::select(treatment, ttT) %>% 
    mutate(file_output = str_c(getwd(), "/", tempfolder))
  
  complete_tables_wholegenome = args_whole_genome %>%
    pmap(gsea_webgestalt) %>% 
    map(~ .x %>% filter(FDR < 0.05)) %>%
    set_names(args_whole_genome$treatment)
  complete_tables_wholegenome %>% saveRDS(file="/home/share/projects/aging/complete_tables_whole_genome_5.rds")}else{
    complete_tables_wholegenome= readRDS(file="/home/share/projects/aging/complete_tables_whole_genome_5.rds")
  }



if(keepallpathways==TRUE){ 
  # extract sets of pathway for each treatment
  gsea_genesetnames2 = 
    complete_tables_wholegenome %>%
    map(~ .$geneSet)
  
  # union of all the pathway
  complete_tables2 =
    complete_tables_wholegenome %>%
    reduce(rbind) %>%
    select(geneSet, description) %>% 
    unique()
  
  # calculate the complement of the pathway for each treatment
  #gsea_genesetnames_complement = map(gsea_genesetnames, ~setdiff(complete_tables$geneSet, .x)) # the complement of these reactome terms
  
  gsea_genesetnames_complement2 = map(gsea_genesetnames2, ~setdiff(complete_tables2$geneSet, .x)) # the complement of these reactome terms
  
  # the universe of reactome terms considered here is "all reactome sets deemed significantly related to at least one ses predictor"
  # there are then 2^5 intersections in the venn diagram (i.e. intersections over the 5 ses predictors, with each being the complement or not)
  
  universe= 
    list(gsea_genesetnames_complement2,
         gsea_genesetnames2) %>% 
    transpose() 
  
  get_venn_cell = 
    function(P){ 
      
      # Get the intersection for each collection of "literals" of ABCDE (a "literal"
      # of A is either A or the complement of A, i.e. A') for each row in matrix
      # of indexes P. These fill the 2^D cells of the Ven diagram.
      map(1:dim(P)[1], 
          # pluck the corresponding subset 
          ~ map2(names(universe),
                 P[.x, ], 
                 ~ pluck(universe, .x, .y + 1)
          ) 
      ) %>% 
        # interestion over all
        map(reduce, intersect)
    }
  
  crossing(!!!set_names(rerun(length(universe), 0:1), 
                        names(universe))) %>% 
    mutate(geneSet = get_venn_cell(.)) %>% 
    unnest(geneSet) %>% 
    left_join(complete_tables2, by = "geneSet")  %>% 
    knitr::kable()
  
  
  
#  venn::venn(gsea_genesetnames2, ilabels = TRUE, 
 #            zcolor = "style", ellipse = FALSE,
  #           opacity = 0.15)
}else{
  #now selecting only the pathways which have at least one leading gene in the aging signature
  # extract sets of pathway for each treatment
  
  complete_tables_wholegenome2 = map(1:length(complete_tables_wholegenome), ~ complete_tables_wholegenome[[.x]]$userId %>%  set_names(complete_tables_wholegenome[[.x]]$geneSet) %>% map(str_split, ";")) %>% flatten()
  #complete_tables_wholegenome %>% rlist::list.any( %in% signatures$outcome_set$aging_mRNA)
  
  bla<-list.flatten(complete_tables_wholegenome2, use.names = TRUE)
  list_of_treat <- list()
  for (i in 1:183){
    list_of_treat[[i]] <- sum(bla[[i]]%in%signatures$outcome_set$aging_mRNA)
    names(list_of_treat[[i]]) =(names(bla[i]))
  }
  physioleading<-cbind(list_of_treat) %>% unlist() %>% as_tibble()
  
  physioleading <-physioleading %>% mutate(physio=names(bla))
  overlapping <-physioleading %>% dplyr::filter(value>0) %>% dplyr::select(physio) %>% c()
  vect <-overlapping$physio
  library(purrr)
  complete_tables_wholegenome3<-map(complete_tables_wholegenome, ~filter(.x, geneSet %in% vect))
  
  
  gsea_genesetnames2 = 
    complete_tables_wholegenome3 %>%
    map(~ .$geneSet)
  
  # union of all the pathway
  complete_tables2 =
    complete_tables_wholegenome3 %>%
    reduce(rbind) %>%
    select(geneSet, description) %>% 
    unique()
  
  
  # calculate the complement of the pathway for each treatment
  #gsea_genesetnames_complement = map(gsea_genesetnames, ~setdiff(complete_tables$geneSet, .x)) # the complement of these reactome terms
  
  gsea_genesetnames_complement2 = map(gsea_genesetnames2, ~setdiff(complete_tables2$geneSet, .x)) # the complement of these reactome terms
  
  # the universe of reactome terms considered here is "all reactome sets deemed significantly related to at least one ses predictor"
  # there are then 2^5 intersections in the venn diagram (i.e. intersections over the 5 ses predictors, with each being the complement or not)
  
  universe= 
    list(gsea_genesetnames_complement2,
         gsea_genesetnames2) %>% 
    transpose() 
  
  get_venn_cell = 
    function(P){ 
      
      # Get the intersection for each collection of "literals" of ABCDE (a "literal"
      # of A is either A or the complement of A, i.e. A') for each row in matrix
      # of indexes P. These fill the 2^D cells of the Ven diagram.
      map(1:dim(P)[1], 
          # pluck the corresponding subset 
          ~ map2(names(universe),
                 P[.x, ], 
                 ~ pluck(universe, .x, .y + 1)
          ) 
      ) %>% 
        # interestion over all
        map(reduce, intersect)
    }
  
  crossing(!!!set_names(rerun(length(universe), 0:1), 
                        names(universe))) %>% 
    mutate(geneSet = get_venn_cell(.)) %>% 
    unnest(geneSet) %>% 
    left_join(complete_tables2, by = "geneSet")  %>% 
    rename(Education=edu_max,
           Income=income_hh_ff5,
           Occupation= SEI_ff5,
           "SES Composite" = ses_sss_composite,
           "Subjective Social Status" = sss_5) %>%  
    knitr::kable()
  
  names(gsea_genesetnames2)[1] <- "Education"
  names(gsea_genesetnames2)[2] <- "Income"
  names(gsea_genesetnames2)[3] <- "Occupation"
  names(gsea_genesetnames2)[4] <- "SES Composite"
  names(gsea_genesetnames2)[5] <- "Subjective Social Status"
  
  
 # venn::venn(gsea_genesetnames2, ilabels = TRUE, 
#             zcolor = "style", ellipse = T,
  #           opacity = 0.15)
}

#'# Detail on Physiological Pathways
crossing(!!!set_names(rerun(length(universe), 0:1), 
                      names(universe))) %>% 
  mutate(geneSet = get_venn_cell(.)) %>% 
  unnest(geneSet) %>% 
  left_join(complete_tables2, by = "geneSet")  %>% 
  rename(Education=edu_max,
         Income=income_hh_ff5,
         Occupation= SEI_ff5,
         "SES Composite" = ses_sss_composite,
         "Subjective Social Status" = sss_5) %>%  
  knitr::kable()


set.seed(123)
library(here)
library(tidyverse)
library(rlang)
library(skimr)
library(furrr)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase)
library(enrichplot)
library(stringi)
library(dbr) # my package
library(openxlsx)

walk(dir(path = here("R"),full.names = TRUE), source)
fit_m4 = partial(fit_m4, n_perm = 1000) # specify n_perm

# WHICH EXAMPLES TO RUN? 
example4 <- example3 <- example2 <- example1 <- example0 <- FALSE
example4 <- TRUE


load_data(reconciled = FALSE, remove_inflam = TRUE)
define_treatments_and_controls()
recode_variables_in_dat_bespoke()
print(abbreviations)
funcs = str_subset(abbreviations$shorthand, "^m") 

#'## Gene set rotation test

library(openxlsx)

#'### SES Composite

de<-de_and_tfbm(treatment = "ses_sss_composite", controls = controls$all, de_only=F)
de$gene_set_test %>%  
  knitr::kable()

#'### Income

de<-de_and_tfbm(treatment = "income_hh_ff5", controls = controls$all, de_only=F)
de$gene_set_test %>%  
  knitr::kable()

#'### Occupation

de<-de_and_tfbm(treatment = "SEI_ff5", controls = controls$all, de_only=F)
de$gene_set_test %>%  
  knitr::kable()

#'### SSS

de<-de_and_tfbm(treatment = "sss_5", controls = controls$all, de_only=F)
de$gene_set_test %>%  
  knitr::kable()

#'### Education

de<-de_and_tfbm(treatment = "edu_max", controls = controls$all, de_only=F)
de$gene_set_test %>%  
  knitr::kable()
