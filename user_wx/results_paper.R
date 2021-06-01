#' ---
#' title: Results of Association between Aging Signature (and its clusters) and SES Indicators
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#' code_folding : hide
#+ warning=FALSE, message=FALSE, include=FALSE
recomputeexamplefiles<-F

if(recomputeexamplefiles==TRUE){ 
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

walk(dir(path = here("R"),full.names = TRUE), source)
fit_m4 = partial(fit_m4, n_perm = 1000) # specify n_perm

# WHICH EXAMPLES TO RUN? 
example4 <- example3 <- example2 <- example1 <- example0 <- FALSE
example4 <- TRUE

############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################

load_data(reconciled = FALSE, remove_inflam = TRUE)
define_treatments_and_controls()
recode_variables_in_dat_bespoke()
print(abbreviations)
funcs = str_subset(abbreviations$shorthand, "^m") 
# explicitly assign ncomp as the smallest number of table signatures gene numbers

ncomp = signatures$outcome_set[table1]%>% map_dfc(length) %>% unlist() %>%  min
ncomp = 5
fit_pca_util = partial(fit_pca_util, ncomp = 5) # specify n_perm

example =
  args %>%
  filter(gene_set_name == "whole_genome_and_tfbm",
         names(controls) == "all") %>%
  mutate(out = pmap(., safely(model_fit), funcs),
         controls = names(controls))
example %>% saveRDS(file="/home/share/projects/aging/tt_whole_genome_5.rds")
funcs="m8"
example0 =
  args %>%
  filter(is.element(gene_set_name, table1),
         names(controls) == "all") %>% 
  mutate(out = pmap(., safely(model_fit), funcs),
         controls = names(controls))

saveRDS(example0, "/home/share/projects/aging/example0_mediation_all_5_tot.rds")
}


# figure 1
# 3 panels

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
library(ggpubr)
library(here)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(tidyverse)
library(stringi)
walk(dir(path = here("R"),full.names = TRUE), source)

example0 = readRDS("/home/share/projects/aging/example0_mediation_all_5_tot.rds")

# example0 = readRDS("/home/share/scratch/xu/example0_w5bmi.rds")
# example0 = readRDS("/home/share/scratch/xu/example0_w5bmi_removeinflam_6pc_trial2.rds")
# panel A
ex0 <- example0 %>%
  hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
  filter(controls=="all") %>% 
  dplyr::select(treatment, gene_set_name, p) %>% 
  dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5"))

ex_0 <- ex0 %>%
  mutate(pval=case_when(p<0.0001 ~ 0.0001,
                        p<0.001 ~ 0.001,
                        p<0.01 ~ 0.01,
                        p<0.05 ~ 0.05,
                        p>0.05 ~ 100),
         pval2=case_when(p<0.0001 ~ 100000,
                         p<0.001 ~ 25000,
                         p<0.01 ~ 15000,
                         p<0.05 ~ 10000,
                         p>0.05 ~ 0.0000001),
         treatment= case_when(treatment =="edu_max" ~ "Education" ,
                              treatment =="income_hh_ff5" ~ "Income"     ,
                              treatment =="SEI_ff5" ~ "Occupation"      ,
                              treatment =="ses_sss_composite" ~ "SES Composite"  ,
                              treatment =="sss_5" ~ "Subjective Social Status",
                              treatment =="ses_composite_ff5"  ~ "SES Composite 3"))  %>%
  mutate(treatment = factor(treatment, levels = c("SES Composite 3", "SES Composite", "Education", "Income",
                                                  "Occupation", "Subjective Social Status","Occupation Work Collar"))) %>%
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


panelA = ggplot(ex_0, aes(treatment, gene_set_name, size = pval2, alpha = 0.4)) +
  geom_point(fill = "red", color="navy") +
  theme_bw() +
  labs(
    # title = "Figure 1. Associations between Indicators of Socioeconomic Status 
    #         and mRNA-Based Disease Signatures, Add Health 
    #         (p-values reported, FDR-corrected for whole genome)",
    y = "Aging mRNA Clusters",
    x = "SES Indicators") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_size_continuous(name = "P-value",
                        range = c(0, 14),
                        limits = c(0.0000001, 100000), breaks = c(0.0000001, 10000, 15000, 25000, 100000),
                        labels = c("n.s.", "p<0.05", "p<0.01", "p<0.001", "p<0.0001"))+
  scale_alpha(guide = 'none')

title <- expression(atop(bold("Figure 1:"),
                         scriptstyle("panel A:Associations between Indicators of Socioeconomic Status and mRNA-Based Disease Signatures, Add Health (p-values reported, FDR-corrected for whole genome)")))

# title <- expression(atop(bold("Figure 1:"),
# scriptstyle("panel A:Associations between Indicators of Socioeconomic Status and mRNA-Based Disease Signatures, Add Health (p-values reported, FDR-corrected for whole genome) \n panel B: Significant PCs \n panel C: pathways")))

annotate_figure(panelA,bottom = text_grob(title)
                # bottom = text_grob("Figure 1: \n panel A: Associations between Indicators of Socioeconomic Status and mRNA-Based Disease Signatures, Add Health (p-values reported, FDR-corrected for whole genome \n panel B: Significant PCs \n panel C: pathways")
                #           fig.lab = "Figure 1:
                #           panel A - Associations between Indicators of Socioeconomic Status
                #      and mRNA-Based Disease Signatures, Add Health
                # (p-values reported, FDR-corrected for whole genome)",
                #           fig.lab.face = "bold",
                #           fig.lab.pos = "bottom.right"
)

#'## Figure 1
#'
panelA


example = readRDS("/home/share/projects/aging/example0_mediation_all_5_tot.rds")

temp =
  example %>%
  hoist(out, ttT = list("result", "m8_fdr", 1, "other", "m")) %>%
  dplyr::filter(ttT != "NULL") %>%
  mutate(gene_sig = map(ttT, ~ dplyr::filter(., adj.P.Val<0.05) %>% pull(gene))) %>%
  dplyr::filter(controls == "all") %>%
  dplyr::filter(gene_sig %>% map_dfc(~ length(.))>0) %>%
  unnest_longer(gene_sig)

temp2<-temp %>%
  hoist(out, p_omnibus = list("result", "m8_fdr", 1, "p")) %>%
  hoist(out, logFC = list("result", "m8_fdr", 1, "detail", "logFC")) %>%
  hoist(out, med_single=list("result", "m8_fdr", 1, "mediation_single")) %>%
  unnest_longer(med_single) %>%
  hoist(med_single, med_single_p = list("result", "p")) %>%
  filter(!is.na(med_single_p)) %>%
  hoist(med_single, med_single_prop = list("result", "other", "med_prop")) %>%
  hoist(med_single, med_single_beta = list("result", "other", "med_beta")) %>%
  dplyr::select(-controls, -ttT, -out, -med_single, -gene_sig) %>%
  mutate(med_single_p= med_single_p %>% str_remove("<") %>% as.numeric()) %>% 
  distinct() %>% 
  #filter(med_single_p<0.05) %>%
  mutate(mediator=stri_extract(med_single_id, regex='[^_]*'),
         #genename=stri_extract(med_single_id, regex='[_^]*'),
         pvalueadj=p.adjust(med_single_p, method="BH")) %>% 
  group_by(treatment, gene_set_name,mediator) %>% 
  summarise(min=min(pvalueadj)) %>% 
  dplyr::filter(min<0.05) %>% 
  kableExtra::kable(caption="Mediation by single gene") %>%
  kableExtra::kable_styling()

temp<-temp %>%
  hoist(out, p_omnibus = list("result", "m8_fdr", 1, "p")) %>%
  hoist(out, logFC = list("result", "m8_fdr", 1, "detail", "logFC")) %>%
  hoist(out, med_single=list("result", "m8_fdr", 1, "mediation_single")) %>%
  unnest_longer(med_single) %>%
  hoist(med_single, med_single_p = list("result", "p")) %>%
  filter(!is.na(med_single_p)) %>%
  hoist(med_single, med_single_prop = list("result", "other", "med_prop")) %>%
  hoist(med_single, med_single_beta = list("result", "other", "med_beta")) %>%
  dplyr::select(-controls, -ttT, -out, -med_single) %>%
  mutate(med_single_p= med_single_p %>% str_remove("<") %>% as.numeric()) %>% 
  distinct() %>% 
  #filter(med_single_p<0.05) %>%
  mutate(mediator=stri_extract(med_single_id, regex='[^_]*'),
         pvalueadj=p.adjust(med_single_p, method="BH")) %>% 
  group_by(treatment, gene_set_name,mediator) %>% 
  summarise(min=min(pvalueadj)) %>% 
  dplyr::filter(min<0.05) 

temp_bmi <- temp %>% dplyr::filter(mediator=="w5bmi")
temp_smoke <- temp %>% dplyr::filter(mediator=="currentsmoke")
temp_bills <- temp %>% dplyr::filter(mediator=="bills")

#bmi mediator
exB <- temp_bmi %>% 
  mutate(pval2=case_when(min<0.0001 ~ 100000,
                         min<0.001 ~ 25000,
                         min<0.01 ~ 15000,
                         min<0.05 ~ 10000,
                         min>0.05 ~ 0.0000001)
  ) %>% 
  dplyr::select(treatment, gene_set_name, mediator, min, pval2) #%>%  dplyr::filter(p < threshold_med) 

exB_data_with = exB %>% 
  mutate(treatment= case_when(treatment =="work_collar_ff5" ~ "Occupation Work Collar",
                              treatment =="edu_max" ~ "Education" ,
                              treatment =="income_hh_ff5" ~ "Income"     ,
                              treatment =="SEI_ff5" ~ "Occupation"      ,
                              treatment =="ses_sss_composite" ~ "SES Composite"  ,
                              treatment =="sss_5" ~ "Subjective Social Status",
                              treatment =="ses_composite_ff5"  ~ "SES Composite 3"))  %>%
  mutate(treatment = factor(treatment, levels = c("SES Composite 3", "SES Composite", "Education", "Income",
                                                  "Occupation", "Subjective Social Status","Occupation Work Collar"
  ))) %>%
  mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
         gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
         gene_set_name= case_when(gene_set_name =="Aging Up Cl4" ~ "Lysosome Metabolism" ,
                                  gene_set_name =="Aging Up Cl1" ~ "Innate/Adaptive Immunity",
                                  gene_set_name =="Aging Down Cl3" ~ "Immune Genes"      ,
                                  gene_set_name =="Aging Down Cl2" ~ "Ribosome"  ,
                                  gene_set_name =="Aging Down Cl1b" ~ "Mitochondrial",
                                  gene_set_name =="Aging" ~ "Aging All Genes"))


#exB_data_with = exB_data_with %>% mutate(gene_set_name = factor(gene_set_name, levels = c("Aging", "Aging Up Cl4", "Aging Up Cl3","Aging Up Cl2",
#                                                                                "Aging Up Cl1", "Aging Down Cl3",
#                                                                                "Aging Down Cl2", "Aging Down Cl1c", "Aging Down Cl1b",
#                                                                                "Aging Down Cl1a")) %>% fct_rev)

exB_data_with = exB_data_with %>%  mutate(treatment = factor(treatment, levels = c("SES Composite",
                                                                         "Education",
                                                                         "Income",
                                                                         "Occupation",
                                                                         "Subjective Social Status")))
#current smoke
exB2 <- temp_smoke %>% 
  mutate(pval2=case_when(min<0.0001 ~ 100000,
                         min<0.001 ~ 25000,
                         min<0.01 ~ 15000,
                         min<0.05 ~ 10000,
                         min>0.05 ~ 0.0000001)
  ) %>% 
  dplyr::select(treatment, gene_set_name, mediator, min, pval2) #%>%  dplyr::filter(p < threshold_med) 

exB_data_with2 = exB2 %>% 
  mutate(treatment= case_when(treatment =="work_collar_ff5" ~ "Occupation Work Collar",
                              treatment =="edu_max" ~ "Education" ,
                              treatment =="income_hh_ff5" ~ "Income"     ,
                              treatment =="SEI_ff5" ~ "Occupation"      ,
                              treatment =="ses_sss_composite" ~ "SES Composite"  ,
                              treatment =="sss_5" ~ "Subjective Social Status",
                              treatment =="ses_composite_ff5"  ~ "SES Composite 3"))  %>%
  mutate(treatment = factor(treatment, levels = c("SES Composite 3", "SES Composite", "Education", "Income",
                                                  "Occupation", "Subjective Social Status","Occupation Work Collar"
  ))) %>%
  mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
         gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
         gene_set_name= case_when(gene_set_name =="Aging Up Cl4" ~ "Lysosome Metabolism" ,
                                  gene_set_name =="Aging Up Cl1" ~ "Innate/Adaptive Immunity",
                                  gene_set_name =="Aging Down Cl3" ~ "Immune Genes"      ,
                                  gene_set_name =="Aging Down Cl2" ~ "Ribosome"  ,
                                  gene_set_name =="Aging Down Cl1b" ~ "Mitochondrial",
                                  gene_set_name =="Aging" ~ "Aging All Genes"))

#bills
#current smoke
exB3 <- temp_bills %>% 
  mutate(pval2=case_when(min<0.0001 ~ 100000,
                         min<0.001 ~ 25000,
                         min<0.01 ~ 15000,
                         min<0.05 ~ 10000,
                         min>0.05 ~ 0.0000001)
  ) %>% 
  dplyr::select(treatment, gene_set_name, mediator, min, pval2) #%>%  dplyr::filter(p < threshold_med) 

exB_data_with3 = exB3 %>% 
  mutate(treatment= case_when(treatment =="work_collar_ff5" ~ "Occupation Work Collar",
                              treatment =="edu_max" ~ "Education" ,
                              treatment =="income_hh_ff5" ~ "Income"     ,
                              treatment =="SEI_ff5" ~ "Occupation"      ,
                              treatment =="ses_sss_composite" ~ "SES Composite"  ,
                              treatment =="sss_5" ~ "Subjective Social Status",
                              treatment =="ses_composite_ff5"  ~ "SES Composite 3"))  %>%
  mutate(treatment = factor(treatment, levels = c("SES Composite 3", "SES Composite", "Education", "Income",
                                                  "Occupation", "Subjective Social Status","Occupation Work Collar"
  ))) %>%
  mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
         gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
         gene_set_name= case_when(gene_set_name =="Aging Up Cl4" ~ "Lysosome Metabolism" ,
                                  gene_set_name =="Aging Up Cl1" ~ "Innate/Adaptive Immunity",
                                  gene_set_name =="Aging Down Cl3" ~ "Immune Genes"      ,
                                  gene_set_name =="Aging Down Cl2" ~ "Ribosome"  ,
                                  gene_set_name =="Aging Down Cl1b" ~ "Mitochondrial",
                                  gene_set_name =="Aging" ~ "Aging All Genes"))

exB_data = full_join(exB_data_with %>% dplyr::select(1,2,3,4,5),
                     exB_data_with2 %>% dplyr::select(1,2,3,4,5), 
                     by.x = treatment, by.y = "gene_set_name") %>% 
  full_join(exB_data_with3 %>% dplyr::select(1,2,3,4,5), 
            by.x = treatment, by.y = "gene_set_name")

# color order from bottom to top
axiscolor = c("grey30", "grey30", "grey30", "darkseagreen1","darkseagreen4" ,"grey30","red", "darkred")

panelB <-ggplot(exB_data, aes(treatment, gene_set_name, size = pval2, fill = mediator, colour = mediator)) +
  geom_point(stroke = 1.5, shape = 21, alpha = 0.5, position = position_dodge(0.5)) +
  scale_color_manual(values = c("deeppink4","goldenrod3","darkblue")) + 
  theme_bw() +
  labs(
    # title = "Figure 1. Associations between Indicators of Socioeconomic Status 
    #         and mRNA-Based Disease Signatures, Add Health 
    #         (p-values reported, FDR-corrected for whole genome)",
    y = "Aging mRNA Clusters",
    x = "SES Indicators") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        text = element_text(size=10, face = "bold")) +
  scale_size_continuous(name = "P-value",
                        range = c(0, 16),
                        limits = c(10000, 100000), breaks = c(10000, 15000, 25000, 100000),
                        labels = c("p<0.05", "p<0.01", "p<0.001", "p<0.0001"))+
  scale_alpha(guide = 'none')+
  guides(#shape = guide_legend(override.aes = list(size = 10)),
    fill = guide_legend(override.aes = list(size = 10)))+
  scale_fill_manual(values = c("deeppink4","goldenrod3","darkblue")) 

#'## Figure 2
panelB


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


#'## Venn Diagram

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
  
  
  
  venn::venn(gsea_genesetnames2, ilabels = TRUE, 
             zcolor = "style", ellipse = FALSE,
             opacity = 0.15)
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
  
  
  venn::venn(gsea_genesetnames2, ilabels = TRUE, 
             zcolor = "style", ellipse = T,
             opacity = 0.15)
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

#' Direct effect!

example = readRDS("/home/share/projects/aging/example0_mediation_all_5_tot.rds")

temp =
  example %>%
  hoist(out, ttT = list("result", "m8_fdr", 1, "other", "m")) %>%
  dplyr::filter(ttT != "NULL") %>%
  mutate(gene_sig = map(ttT, ~ dplyr::filter(., adj.P.Val<0.05) %>% pull(gene))) %>%
  dplyr::filter(controls == "all") %>%
  dplyr::filter(gene_sig %>% map_dfc(~ length(.))>0) %>%
  unnest_longer(gene_sig)

temp2<-temp %>%
  hoist(out, p_omnibus = list("result", "m8_fdr", 1, "p")) %>%
  hoist(out, logFC = list("result", "m8_fdr", 1, "detail", "logFC")) %>%
  hoist(out, med_single=list("result", "m8_fdr", 1, "mediation_single")) %>%
  unnest_longer(med_single) %>%
  hoist(med_single, med_single_p = list("result", "p")) %>%
  dplyr::filter(!is.na(med_single_p)) %>%
  hoist(med_single, med_single_prop = list("result", "other", "med_prop")) %>%
  hoist(med_single, med_single_beta = list("result", "other", "med_beta")) %>%
  hoist(med_single, med_single_ade_p = list("result", "other", "med_ADE_p")) %>%
  dplyr::select(-controls, -ttT, -out, -med_single, -gene_sig) %>%
  mutate(med_single_p= med_single_p %>% str_remove("<") %>% as.numeric()) %>% 
  distinct() %>% 
  #filter(med_single_p<0.05) %>%
  mutate(mediator=stri_extract(med_single_id, regex='[^_]*'),
         #genename=stri_extract(med_single_id, regex='[_^]*'),
         pvalueadj=p.adjust(med_single_p, method="BH")) %>% 
  group_by(treatment, gene_set_name,mediator) %>% 
  summarise(min=min(pvalueadj)) %>% 
  dplyr::filter(min<0.05) %>% 
  kableExtra::kable(caption="Mediation by single gene") %>%
  kableExtra::kable_styling()

temp3<-temp %>%
  hoist(out, p_omnibus = list("result", "m8_fdr", 1, "p")) %>%
  hoist(out, logFC = list("result", "m8_fdr", 1, "detail", "logFC")) %>%
  hoist(out, med_single=list("result", "m8_fdr", 1, "mediation_single")) %>%
  unnest_longer(med_single) %>%
  hoist(med_single, med_single_p = list("result", "p")) %>%
  dplyr::filter(!is.na(med_single_p)) %>%
  hoist(med_single, med_single_prop = list("result", "other", "med_prop")) %>%
  hoist(med_single, med_single_beta = list("result", "other", "med_beta")) %>%
  hoist(med_single, med_single_ade_p = list("result", "other", "med_ADE_p")) %>%
  dplyr::select(-controls, -ttT, -out, -med_single) %>%
  mutate(med_single_p= med_single_p %>% str_remove("<") %>% as.numeric()) %>%
  mutate(med_single_ade_p= med_single_ade_p %>% str_remove("<") %>% as.numeric()) %>% 
  distinct() %>% 
  #filter(med_single_p<0.05) %>%
  mutate(mediator=stri_extract(med_single_id, regex='[^_]*'),
         gene1=stri_extract_last(med_single_id, regex='([^_]+$)'),
         pvalueadj_ade=p.adjust(med_single_ade_p, method="BH"))%>%
  dplyr::filter(pvalueadj_ade<0.05) %>% 
  group_by(treatment, gene_set_name,mediator) %>% 
  summarise(min=min(pvalueadj_ade)) %>% 
  dplyr::filter(min<0.05) 

temp4<-temp %>%
  hoist(out, p_omnibus = list("result", "m8_fdr", 1, "p")) %>%
  hoist(out, logFC = list("result", "m8_fdr", 1, "detail", "logFC")) %>%
  hoist(out, med_single=list("result", "m8_fdr", 1, "mediation_single")) %>%
  unnest_longer(med_single) %>%
  hoist(med_single, med_single_p = list("result", "p")) %>%
  dplyr::filter(!is.na(med_single_p)) %>%
  hoist(med_single, med_single_prop = list("result", "other", "med_prop")) %>%
  hoist(med_single, med_single_beta = list("result", "other", "med_beta")) %>%
  hoist(med_single, med_single_ade_p = list("result", "other", "med_ADE_p")) %>%
  dplyr::select(-controls, -ttT, -out, -med_single) %>%
  mutate(med_single_p= med_single_p %>% str_remove("<") %>% as.numeric()) %>%
  mutate(med_single_ade_p= med_single_ade_p %>% str_remove("<") %>% as.numeric()) %>% 
  distinct() %>% 
  #filter(med_single_p<0.05) %>%
  mutate(mediator=stri_extract(med_single_id, regex='[^_]*'),
         gene1=stri_extract_last(med_single_id, regex='([^_]+$)'),
         pvalueadj=p.adjust(med_single_p, method="BH"))%>%
  filter(pvalueadj<0.05) %>% 
  group_by(treatment, gene_set_name,mediator, gene1) %>% 
  summarise(min=min(pvalueadj)) %>% 
  distinct(treatment, gene_set_name,mediator, min) %>% 
  dplyr::filter(min<0.05) 


temp_bmi <- temp3 %>% dplyr::filter(mediator=="w5bmi")
temp_smoke <- temp3 %>% dplyr::filter(mediator=="currentsmoke")
temp_bills <- temp3 %>% dplyr::filter(mediator=="bills")

#bmi mediator
exB <- temp_bmi %>% 
  mutate(pval2=case_when(min<0.0001 ~ 100000,
                         min<0.001 ~ 25000,
                         min<0.01 ~ 15000,
                         min<0.05 ~ 10000,
                         min>0.05 ~ 0.0000001)
  ) %>% 
  dplyr::select(treatment, gene_set_name, mediator, min, pval2) #%>%  dplyr::filter(p < threshold_med) 

exB_data_with = exB %>% 
  mutate(treatment= case_when(treatment =="work_collar_ff5" ~ "Occupation Work Collar",
                              treatment =="edu_max" ~ "Education" ,
                              treatment =="income_hh_ff5" ~ "Income"     ,
                              treatment =="SEI_ff5" ~ "Occupation"      ,
                              treatment =="ses_sss_composite" ~ "SES Composite"  ,
                              treatment =="sss_5" ~ "Subjective Social Status",
                              treatment =="ses_composite_ff5"  ~ "SES Composite 3"))  %>%
  mutate(treatment = factor(treatment, levels = c("SES Composite 3", "SES Composite", "Education", "Income",
                                                  "Occupation", "Subjective Social Status","Occupation Work Collar"
  ))) %>%
  mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
         gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
         gene_set_name= case_when(gene_set_name =="Aging Up Cl4" ~ "Lysosome Metabolism" ,
                                  gene_set_name =="Aging Up Cl1" ~ "Innate/Adaptive Immunity",
                                  gene_set_name =="Aging Down Cl3" ~ "Immune Genes"      ,
                                  gene_set_name =="Aging Down Cl2" ~ "Ribosome"  ,
                                  gene_set_name =="Aging Down Cl1b" ~ "Mitochondrial",
                                  gene_set_name =="Aging" ~ "Aging All Genes"))


exB_data_with = exB_data_with %>%  mutate(treatment = factor(treatment, levels = c("SES Composite",
                                                                                   "Education",
                                                                                   "Income",
                                                                                   "Occupation",
                                                                                   "Subjective Social Status")))
#current smoke
exB2 <- temp_smoke %>% 
  mutate(pval2=case_when(min<0.0001 ~ 100000,
                         min<0.001 ~ 25000,
                         min<0.01 ~ 15000,
                         min<0.05 ~ 10000,
                         min>0.05 ~ 0.0000001)
  ) %>% 
  dplyr::select(treatment, gene_set_name, mediator, min, pval2) #%>%  dplyr::filter(p < threshold_med) 

exB_data_with2 = exB2 %>% 
  mutate(treatment= case_when(treatment =="work_collar_ff5" ~ "Occupation Work Collar",
                              treatment =="edu_max" ~ "Education" ,
                              treatment =="income_hh_ff5" ~ "Income"     ,
                              treatment =="SEI_ff5" ~ "Occupation"      ,
                              treatment =="ses_sss_composite" ~ "SES Composite"  ,
                              treatment =="sss_5" ~ "Subjective Social Status",
                              treatment =="ses_composite_ff5"  ~ "SES Composite 3"))  %>%
  mutate(treatment = factor(treatment, levels = c("SES Composite 3", "SES Composite", "Education", "Income",
                                                  "Occupation", "Subjective Social Status","Occupation Work Collar"
  ))) %>%
  mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
         gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
         gene_set_name= case_when(gene_set_name =="Aging Up Cl4" ~ "Lysosome Metabolism" ,
                                  gene_set_name =="Aging Up Cl1" ~ "Innate/Adaptive Immunity",
                                  gene_set_name =="Aging Down Cl3" ~ "Immune Genes"      ,
                                  gene_set_name =="Aging Down Cl2" ~ "Ribosome"  ,
                                  gene_set_name =="Aging Down Cl1b" ~ "Mitochondrial",
                                  gene_set_name =="Aging" ~ "Aging All Genes"))

#bills
#current smoke
exB3 <- temp_bills %>% 
  mutate(pval2=case_when(min<0.0001 ~ 100000,
                         min<0.001 ~ 25000,
                         min<0.01 ~ 15000,
                         min<0.05 ~ 10000,
                         min>0.05 ~ 0.0000001)
  ) %>% 
  dplyr::select(treatment, gene_set_name, mediator, min, pval2) #%>%  dplyr::filter(p < threshold_med) 

exB_data_with3 = exB3 %>% 
  mutate(treatment= case_when(treatment =="work_collar_ff5" ~ "Occupation Work Collar",
                              treatment =="edu_max" ~ "Education" ,
                              treatment =="income_hh_ff5" ~ "Income"     ,
                              treatment =="SEI_ff5" ~ "Occupation"      ,
                              treatment =="ses_sss_composite" ~ "SES Composite"  ,
                              treatment =="sss_5" ~ "Subjective Social Status",
                              treatment =="ses_composite_ff5"  ~ "SES Composite 3"))  %>%
  mutate(treatment = factor(treatment, levels = c("SES Composite 3", "SES Composite", "Education", "Income",
                                                  "Occupation", "Subjective Social Status","Occupation Work Collar"
  ))) %>%
  mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
         gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
         gene_set_name= case_when(gene_set_name =="Aging Up Cl4" ~ "Lysosome Metabolism" ,
                                  gene_set_name =="Aging Up Cl1" ~ "Innate/Adaptive Immunity",
                                  gene_set_name =="Aging Down Cl3" ~ "Immune Genes"      ,
                                  gene_set_name =="Aging Down Cl2" ~ "Ribosome"  ,
                                  gene_set_name =="Aging Down Cl1b" ~ "Mitochondrial",
                                  gene_set_name =="Aging" ~ "Aging All Genes"))

exB_data = full_join(exB_data_with %>% dplyr::select(1,2,3,4,5),
                     exB_data_with2 %>% dplyr::select(1,2,3,4,5), 
                     by.x = treatment, by.y = "gene_set_name") %>% 
  full_join(exB_data_with3 %>% dplyr::select(1,2,3,4,5), 
            by.x = treatment, by.y = "gene_set_name")

# color order from bottom to top
axiscolor = c("grey30", "grey30", "grey30", "darkseagreen1","darkseagreen4" ,"grey30","red", "darkred")

panelC <-ggplot(exB_data, aes(treatment, gene_set_name, size = pval2, fill = mediator, colour = mediator)) +
  geom_point(stroke = 1.5, shape = 21, alpha = 0.5, position = position_dodge(0.5)) +
  scale_color_manual(values = c("deeppink4","goldenrod3","darkblue")) + 
  theme_bw() +
  labs(
    # title = "Figure 1. Associations between Indicators of Socioeconomic Status 
    #         and mRNA-Based Disease Signatures, Add Health 
    #         (p-values reported, FDR-corrected for whole genome)",
    y = "mRNA Signatures",
    x = "SES Indicators") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        text = element_text(size=10, face = "bold")) +
  scale_size_continuous(name = "P-value",
                        range = c(0, 16),
                        limits = c(10000, 100000), breaks = c(10000, 15000, 25000, 100000),
                        labels = c("p<0.05", "p<0.01", "p<0.001", "p<0.0001"))+
  scale_alpha(guide = 'none')+
  guides(#shape = guide_legend(override.aes = list(size = 10)),
    fill = guide_legend(override.aes = list(size = 10)))+
  scale_fill_manual(values = c("deeppink4","goldenrod3","darkblue")) 

#'## Figure 3 Direct Effect
panelC


#'# Kegg Pathways
#'recomputeexamplefiles<-F


if(recomputeexamplefiles==TRUE){ 
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
  
  walk(dir(path = here("R"),full.names = TRUE), source)
  fit_m4 = partial(fit_m4, n_perm = 1000) # specify n_perm
  
  # WHICH EXAMPLES TO RUN? 
  example4 <- example3 <- example2 <- example1 <- example0 <- FALSE
  example4 <- TRUE
  
  ############################################################
  # LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
  ############################################################
  
  load_data(reconciled = FALSE, remove_inflam = TRUE)
  define_treatments_and_controls()
  recode_variables_in_dat_bespoke()
  print(abbreviations)
  funcs = str_subset(abbreviations$shorthand, "^m") 
  # explicitly assign ncomp as the smallest number of table signatures gene numbers
  
  ncomp = signatures$outcome_set[table1]%>% map_dfc(length) %>% unlist() %>%  min
  ncomp = 5
  fit_pca_util = partial(fit_pca_util, ncomp = 5) # specify n_perm
  
  example =
    args %>%
    filter(gene_set_name == "whole_genome_and_tfbm",
           names(controls) == "all") %>%
    mutate(out = pmap(., safely(model_fit), funcs),
           controls = names(controls))
  example %>% saveRDS(file="/home/share/projects/aging/tt_whole_genome_5.rds")
  funcs="m8"
  example0 =
    args %>%
    filter(is.element(gene_set_name, table1),
           names(controls) == "all") %>% 
    mutate(out = pmap(., safely(model_fit), funcs),
           controls = names(controls))
  
  saveRDS(example0, "/home/share/projects/aging/example0_mediation_all_5_tot.rds")
}

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
library(ggpubr)
library(here)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(knitr)

# figure 1
# 3 panels

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
library(ggpubr)
library(here)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(tidyverse)
library(stringi)
walk(dir(path = here("R"),full.names = TRUE), source)

example0 = readRDS("/home/share/projects/aging/example0_mediation_all_5_tot.rds")

# example0 = readRDS("/home/share/scratch/xu/example0_w5bmi.rds")
# example0 = readRDS("/home/share/scratch/xu/example0_w5bmi_removeinflam_6pc_trial2.rds")
# panel A
ex0 <- example0 %>%
  hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
  filter(controls=="all") %>% 
  dplyr::select(treatment, gene_set_name, p) %>% 
  dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5"))

ex_0 <- ex0 %>%
  mutate(pval=case_when(p<0.0001 ~ 0.0001,
                        p<0.001 ~ 0.001,
                        p<0.01 ~ 0.01,
                        p<0.05 ~ 0.05,
                        p>0.05 ~ 100),
         pval2=case_when(p<0.0001 ~ 100000,
                         p<0.001 ~ 25000,
                         p<0.01 ~ 15000,
                         p<0.05 ~ 10000,
                         p>0.05 ~ 0.0000001),
         treatment= case_when(treatment =="edu_max" ~ "Education" ,
                              treatment =="income_hh_ff5" ~ "Income"     ,
                              treatment =="SEI_ff5" ~ "Occupation"      ,
                              treatment =="ses_sss_composite" ~ "SES Composite"  ,
                              treatment =="sss_5" ~ "Subjective Social Status",
                              treatment =="ses_composite_ff5"  ~ "SES Composite 3"))  %>%
  mutate(treatment = factor(treatment, levels = c("SES Composite 3", "SES Composite", "Education", "Income",
                                                  "Occupation", "Subjective Social Status","Occupation Work Collar"))) %>%
  mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
         gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title())


panelA = ggplot(ex_0, aes(treatment, gene_set_name, size = pval2, alpha = 0.4)) +
  geom_point(fill = "red", color="navy") +
  theme_bw() +
  labs(
    # title = "Figure 1. Associations between Indicators of Socioeconomic Status 
    #         and mRNA-Based Disease Signatures, Add Health 
    #         (p-values reported, FDR-corrected for whole genome)",
    y = "mRNA Signatures",
    x = "SES Indicators") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_size_continuous(name = "P-value",
                        range = c(0, 14),
                        limits = c(0.0000001, 100000), breaks = c(0.0000001, 10000, 15000, 25000, 100000),
                        labels = c("n.s.", "p<0.05", "p<0.01", "p<0.001", "p<0.0001"))+
  scale_alpha(guide = 'none')

title <- expression(atop(bold("Figure 1:"),
                         scriptstyle("panel A:Associations between Indicators of Socioeconomic Status and mRNA-Based Disease Signatures, Add Health (p-values reported, FDR-corrected for whole genome)")))

# title <- expression(atop(bold("Figure 1:"),
# scriptstyle("panel A:Associations between Indicators of Socioeconomic Status and mRNA-Based Disease Signatures, Add Health (p-values reported, FDR-corrected for whole genome) \n panel B: Significant PCs \n panel C: pathways")))

annotate_figure(panelA,bottom = text_grob(title)
                # bottom = text_grob("Figure 1: \n panel A: Associations between Indicators of Socioeconomic Status and mRNA-Based Disease Signatures, Add Health (p-values reported, FDR-corrected for whole genome \n panel B: Significant PCs \n panel C: pathways")
                #           fig.lab = "Figure 1:
                #           panel A - Associations between Indicators of Socioeconomic Status
                #      and mRNA-Based Disease Signatures, Add Health
                # (p-values reported, FDR-corrected for whole genome)",
                #           fig.lab.face = "bold",
                #           fig.lab.pos = "bottom.right"
)



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
keepallpathways<-T

# function gsea_webgestalt calculate functional pathway using preferred database:
# pathway_KEGG and pathway_Reactome for example
# you could change the argument in the function 
gsea_webgestalt = function(treatment, ttT, file_output, enrichMethod = "GSEA", enrichDatabase="pathway_KEGG"){
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
  complete_tables_wholegenome %>% saveRDS(file="/home/share/projects/aging/complete_tables_whole_genome_5_kegg.rds")}else{
    complete_tables_wholegenome= readRDS(file="/home/share/projects/aging/complete_tables_whole_genome_5_kegg.rds")
  }


#'## Venn Diagram

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
  
  
  
  venn::venn(gsea_genesetnames2, ilabels = TRUE, 
             zcolor = "style", ellipse = FALSE,
             opacity = 0.15)
}else{
  #now selecting only the pathways which have at least one leading gene in the aging signature
  # extract sets of pathway for each treatment
  
  complete_tables_wholegenome2 = map(1:length(complete_tables_wholegenome), ~ complete_tables_wholegenome[[.x]]$userId %>%  set_names(complete_tables_wholegenome[[.x]]$geneSet) %>% map(str_split, ";")) %>% flatten()
  #complete_tables_wholegenome %>% rlist::list.any( %in% signatures$outcome_set$aging_mRNA)
  
  bla<-list.flatten(complete_tables_wholegenome2, use.names = TRUE)
  list_of_treat <- list()
  for (i in 1:13){
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
  
  
  venn::venn(gsea_genesetnames2, ilabels = TRUE, 
             zcolor = "style", ellipse = T,
             opacity = 0.15)
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

############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################

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


