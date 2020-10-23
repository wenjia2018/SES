# gsea

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


ses4 <- readxl::read_excel("./user_wx_RESTORED/DE_removetable1signatures_ses4income_fulllist_Brandt.xlsx", sheet = 1)
ses4 = ses4 %>% dplyr::select(1,2)
write.table(ses4, file ="./user_wx/ses4.rnk", quote = FALSE,
            col.names = FALSE, row.names = FALSE,
            sep = "\t")

rankFile <- "./user_wx/ses4.rnk"

outputDirectory <- getwd()
enrichResult_kegg <- WebGestaltR(enrichMethod="GSEA", organism="hsapiens",
                            enrichDatabase="pathway_KEGG", interestGeneFile=rankFile,
                            interestGeneType="genesymbol", sigMethod="top", topThr=10, minNum=5,
                            outputDirectory=outputDirectory)


enrichResult_kegg %>% kableExtra::kable() %>% kableExtra::kable_styling()
# http://www.webgestalt.org/results/1603186800/#
enrichResult_reactome <- WebGestaltR(enrichMethod="GSEA", organism="hsapiens",
                                 enrichDatabase="pathway_Reactome", interestGeneFile=rankFile,
                                 interestGeneType="genesymbol", sigMethod="top", topThr=10, minNum=5,
                                 outputDirectory=outputDirectory)


enrichResult_reactome %>% kableExtra::kable() %>% kableExtra::kable_styling()
# 
# http://www.webgestalt.org/results/1603187038/#




outputDirectory <- getwd()
gsea_webgestalt = function(sheet_name){
  ses <- readxl::read_excel("./user_wx/DE_ses_fulllist_gsea.xlsx", sheet = sheet_name)
  ses = ses %>%
    # filter(!(gene %in%sig)) %>%
    dplyr::select(1,2)

  write.table(ses, file =str_c("./user_wx/ses_", sheet_name,".rnk"),
              quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  
  rankFile <- str_c("./user_wx/ses_", sheet_name,".rnk")

  enrichResult_reactome <- WebGestaltR(enrichMethod="GSEA", organism="hsapiens",
                                       enrichDatabase="pathway_Reactome", interestGeneFile=rankFile,
                                       interestGeneType="genesymbol",
                                       # sigMethod="top", topThr=10,
                                       fdrThr = 1,
                                       minNum=5,
                                       perNum = 1000,
                                       outputDirectory=outputDirectory)
  
}
a = c("ses4", "income", "edu", "SEI", "sss") %>%
  set_names() %>% 
  map(gsea_webgestalt)

a %>% saveRDS("./user_wx/ses_gsea_webgestaltR_removefig1A_10kperm.rds")



