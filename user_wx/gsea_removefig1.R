# gsea removing all signatures in figure1

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
library(WebGestaltR)
library(dbr) # my package
library(MendelianRandomization)
walk(dir(path = here("R"),full.names = TRUE), source)
fit_m4 = partial(fit_m4, n_perm = 1000) # specify n_perm


############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################

load_data(reconciled = FALSE, remove_inflam = FALSE)
define_treatments_and_controls()
recode_variables_in_dat()

sig = Reduce(union, signatures$outcome_set[table1[1:11]])

# ses4 as an example
if(0){
  ses4 <- readxl::read_excel("./user_wx/DE_ses_fulllist_gsea.xlsx", sheet = 1)
  ses4 = ses4 %>%
    filter(!(gene %in%sig)) %>% 
    dplyr::select(1,2)
  write.table(ses4, file ="./user_wx/ses4_removefig1.rnk", quote = FALSE,
              col.names = FALSE, row.names = FALSE,
              sep = "\t")
  
  rankFile <- "./user_wx/ses4_removefig1.rnk"
  
  outputDirectory <- getwd()
  # enrichResult_kegg <- WebGestaltR(enrichMethod="GSEA", organism="hsapiens",
  #                                  enrichDatabase="pathway_KEGG", interestGeneFile=rankFile,
  #                                  interestGeneType="genesymbol", sigMethod="top", topThr=10, minNum=5,
  #                                  outputDirectory=outputDirectory)
  # 
  # 
  # enrichResult_kegg %>% kableExtra::kable() %>% kableExtra::kable_styling()
  set.seed(1234123493)
  enrichResult_reactome <- WebGestaltR(enrichMethod="GSEA", organism="hsapiens",
                                       enrichDatabase="pathway_Reactome", interestGeneFile=rankFile,
                                       interestGeneType="genesymbol",
                                       # sigMethod="top", topThr=10,
                                       fdrThr = 1,
                                       minNum=5,
                                       outputDirectory=outputDirectory)
  
  
  enrichResult_reactome %>% kableExtra::kable() %>% kableExtra::kable_styling() 
}


# 
outputDirectory <- getwd()
gsea_webgestalt = function(sheet_name){
  # ses <- readxl::read_excel("./user_wx/DE_ses_fulllist_gsea.xlsx", sheet = sheet_name)
  # ses = ses %>%
  #   filter(!(gene %in%sig)) %>% 
  #   dplyr::select(1,2)
  # 
  # write.table(ses, file =str_c("./user_wx/ses_", sheet_name,"_removefig1.rnk"),
  #             quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  
  rankFile <- str_c("./user_wx/ses_", sheet_name,"_removefig1.rnk")
  set.seed(1234554321)
  enrichResult_reactome <- WebGestaltR(enrichMethod="GSEA", organism="hsapiens",
                                       enrichDatabase="pathway_KEGG", interestGeneFile=rankFile,
                                       interestGeneType="genesymbol",
                                       # sigMethod="top", topThr=10,
                                       fdrThr = 1,
                                       minNum=5,
                                       perNum = 1000,
                                       outputDirectory=outputDirectory)
  
}
full_gsea = c("ses4", "income", "edu", "SEI", "sss") %>%
  set_names() %>% 
  map(gsea_webgestalt)
# full_gsea %>% readRDS("./user_wx/ses_gsea_webgestaltR_removefig1A.rds")

complete_tables = full_gsea %>% map(~ .x %>% filter(FDR<0.05))

# 
# enrichResult_reactome %>% kableExtra::kable() %>% kableExtra::kable_styling()
gsea_genesetnames = list(
  edu = complete_tables$edu$geneSet,
  income = complete_tables$income$geneSet,
  SEI = complete_tables$SEI$geneSet,
  ses4 = complete_tables$ses4$geneSet,
  sss = complete_tables$sss$geneSet)


complete_tables <- complete_tables %>% reduce(rbind) %>% select(geneSet, description) %>% unique()
gsea_genesetnames_complement <- map(gsea_genesetnames, ~setdiff(complete_tables$geneSet, .x)) # the complement of these reactome terms

# the universe of reactome terms considered here is "all reactome sets deemed significantly related to at least one ses predictor"
# there are then 2^5 intersections in the venn diagram (i.e. intersections over the 5 ses predictors, with each being the complement or not)

universe= 
  list(gsea_genesetnames_complement,
       gsea_genesetnames) %>% 
  transpose()

crossing(edu =0:1, 
         income =0:1, 
         SEI =0:1,
         ses4 =0:1,
         sss =0:1) %>%
  mutate(x = 
           pmap(., 
                function(edu, income, SEI, ses4, sss) 
                  list(universe$edu[edu + 1], 
                       universe$income[income + 1], 
                       universe$SEI[SEI + 1],
                       universe$ses4[ses4 + 1], 
                       universe$sss[sss + 1]))  %>%
           map(flatten)  %>% 
           map(reduce, intersect)) %>% 
  unnest(x) %>% 
  left_join(complete_tables, by = c("x" = "geneSet"))  %>% 
  knitr::kable()



venn::venn(gsea_genesetnames, ilabels = TRUE, 
           zcolor = "style", ellipse = FALSE,
           opacity = 0.15)
venn::venn(gsea_genesetnames, ilabels = FALSE, zcolor = "style",snames = " ",
           ellipse = FALSE, opacity = 0.15, ilcs = 1.2,
           box = FALSE,
           sncs = 0.1)
