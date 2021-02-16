#' ---
#' title: Results of Association between Aging Signature (and its clusters) and SES Indicators
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE

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

############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################
example = readRDS("/home/share/projects/aging/example0_aging_clusters_pc5.rds")
# create a folder in your current directory to save all the intermediate files and report generated at last from
# webgestalt, you could name your new folder differently
tempfolder2 = "temp_webgestalt"
dir.create(tempfolder2)
# create input for the gsea_webgestalt function
# ttT is the results extracted from m8_fdr
args = 
  example %>% 
  hoist(out, ttT = list("result", "m8_fdr", 1, "other", "m")) %>% 
  dplyr::filter(gene_set_name == "aging_mRNA") %>% # for a specific signature 
  dplyr::select(treatment, ttT) %>% 
  mutate(file_output = str_c(getwd(), "/", tempfolder2)) %>% 
  slice(3, 5, 9, 10) # for test
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

# get pathway for FRD < 0.05
# be noted that the results are slightly different every time you run it, as it is based on permutations
# you could change perNum to enlarge the number of permutation
complete_tables = args %>%
  pmap(gsea_webgestalt) %>% 
  map(~ .x %>% filter(FDR < 0.05)) %>%
  set_names(args$treatment)

#BiocManager::install("reactome.db")
library("reactome.db")
library("rlist")
xx <- as.list(reactomePATHID2EXTID)

treatment = c(
  "ses_sss_composite", "income_hh_ff5", "SEI_ff5",  "sss_5"
)


list_of_vecs <- list()
list_of_treat <- list()
list_long <- list()
list_long_ses_sss_composite <- list()
list_long_SEI_ff5 <- list()
list_long_sss_5 <- list()

for (v in treatment){
  list_of_treat[[v]] <- complete_tables[[v]]$geneSet  
}
for (c in treatment){
  list_of_vecs =list() 
  for (i in list_of_treat[[c]]){
    list_of_vecs[[i]]  <-( list.filter(xx[[i]]))
  }
  list_long[[c]]<- (list_of_vecs)
  
}

# convert our whole gene list from hugo name to entrez
results =  readRDS("/home/share/preprocessed_two_batches/entrezgeneid.rds")

a = list_long %>% map_depth(2, function(x) plyr::mapvalues(x, results$entrezgene_id, results$hgnc_symbol, warn_missing = F))
