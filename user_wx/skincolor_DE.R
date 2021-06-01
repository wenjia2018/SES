

#+ echo=F, eval=T, warning=FALSE, message=FALSE

p_eqtl <- c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
control = "ancestryPC_ses"

library(tidyverse)
# functions to extract data
source("/home/xu/ses-1/user_wx/extract_v2.R")

#' ### DE genes for skincolor5 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
color5 <- readRDS("~/ses-1/user_wx/skincolor5_bespoke_DE_01.03.2021.rds")
ttT = p_eqtl %>% map(outttT, control, color5)
temp = ttT %>%
  bind_rows() %>%
  filter(gene_sig %>% map_dfc(~ length(.))>0) %>% 
  mutate(p_eqtl = p_eqtl %>% format(scientific =T)) %>%
  select(-ttT, -controls)
black = temp %>%
  filter(treatment %>% str_detect("Black"))

black$gene_sig

dark = temp %>%
  filter(treatment %>% str_detect("Dark"))

dark$gene_sig
#' ### DE genes for skincolor3 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
color3 <- readRDS("~/ses-1/user_wx/skincolor3_bespoke_DE_01.03.2021.rds")
ttT = p_eqtl %>% map(outttT, control, color3)
temp = ttT %>%
  bind_rows() %>%
  filter(gene_sig %>% map_dfc(~ length(.))>0) %>% 
  mutate(p_eqtl = p_eqtl %>% format(scientific =T)) %>%
  select(-ttT, -controls)
ttT[[1]]$ttT[[1]]
temp$gene_sig[[1]]  %>% openxlsx::write.xlsx("./user_wx/color3_bespoke_DE_gene.xlsx")

file_output = str_c(getwd(),"/temp_webgestalt/")
rankFile = str_c(file_output, "color3_darkblack",".rnk")
ttT[[1]]$ttT[[1]] %>%
  dplyr::select(gene, logFC) %>% 
  write.table(file = rankFile, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

enrichResult = WebGestaltR::WebGestaltR(interestGeneFile = rankFile,
                                        interestGeneType = "genesymbol",
                                        enrichMethod = "GSEA",
                                        enrichDatabase="pathway_Reactome",
                                        organism = "hsapiens",
                                        fdrThr = 1,
                                        minNum=5,
                                        perNum = 1000,
                                        outputDirectory = file_output)
