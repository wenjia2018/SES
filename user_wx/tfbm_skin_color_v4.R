#' ---
#' title: DE and tfbm for skin color
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
library(readxl)
library(parsnip)
library(workflows)
library(Biobase)
library(enrichplot)
library(gridExtra)
library(grid)
library(ggpubr)
library(dbr) # my package
walk(dir(path = here("R"),full.names = TRUE), source)
load_data(reconciled = FALSE, remove_inflam = FALSE)
source("/home/xu/ses-1/user_wx/extract_v2.R")
# example_race = readRDS("~/ses-1/user_wx/race_bespoke_15.03.2021.rds")
example = readRDS("/home/xu/ses-1/user_wx/color3_bespoke_28.03.2021.rds")

example_skincolor3 = example %>% filter(table1 %>% str_detect("aging|whole_genome"))
control = "ancestryPC_ses"
p_eqtl = 0.05
tempfolder = "temp_webgestalt"
#' ## DE and TFBM genowide for skin color
#+ echo=F, eval=T, warning=FALSE, message=FALSE



#' ### DE gene for Darkblack 
#+ echo=F, eval=T, warning=FALSE, message=FALSE

temp = p_eqtl %>%
  outm10_whole_genome(., control, data = example_skincolor3) %>% 
  filter(ttT!="NULL") %>% 
  mutate(DE_gene = map(ttT, ~ filter(., adj.P.Val< 0.05) %>% pull(gene))) %>%
  filter(map_lgl(DE_gene, ~ length(.)>0))
DE_gene = list2(!!str_c(temp$treatment, temp$p_eqtl) := temp %>% pull(DE_gene)) %>% 
  unlist(recursive = F)

venn::venn(list(DE_gene = DE_gene[[1]],
                Aging_Down = signatures$outcome_set$aging_down_mRNA,
                Aging_Up = signatures$outcome_set$aging_up_mRNA
                # Aging_Cluster_Complement = signatures$outcome_set$aging_cluster_complement_mRNA
),
ilabels = FALSE, zcolor = "style",snames = " ",
ellipse = FALSE, opacity = 0.15, ilcs = 1.2,
box = FALSE,
sncs = 0.8)


# darkblack_DE_list = temp$ttT[[1]] %>% dplyr::select(1,2)
# write.table(darkblack_DE_list, file ="./user_wx/darkblack_DE_list.rnk", quote = FALSE,
#             col.names = FALSE, row.names = FALSE,
#             sep = "\t")
# 
# rankFile <- "./user_wx/darkblack_DE_list.rnk"
# 
# outputDirectory <- getwd()
# enrichResult_kegg <- WebGestaltR::WebGestaltR(enrichMethod="GSEA", organism="hsapiens",
#                                  enrichDatabase="pathway_KEGG", interestGeneFile=rankFile,
#                                  interestGeneType="genesymbol", sigMethod="top", topThr=10, minNum=5,
#                                  outputDirectory=outputDirectory)
# 
# enrichResult_kegg %>% saveRDS("./user_wx/darkblack_kegg.rds")
# 
# enrichResult_kegg %>% kableExtra::kable() %>% kableExtra::kable_styling()
# # http://www.webgestalt.org/results/1603186800/#
# enrichResult_reactome <- WebGestaltR::WebGestaltR(enrichMethod="GSEA", organism="hsapiens",
#                                      enrichDatabase="pathway_Reactome", interestGeneFile=rankFile,
#                                      interestGeneType="genesymbol", sigMethod="top", topThr=10, minNum=5,
#                                      outputDirectory=outputDirectory)
# 
# enrichResult_reactome %>% saveRDS("./user_wx/darkblack_reactome.rds")
# enrichResult_reactome %>% kableExtra::kable() %>% kableExtra::kable_styling()
# 
# 
if(reproduciable <- FALSE) {
  temp = p_eqtl %>%
    outm10_whole_genome(., control, data = example_skincolor3)
  rankFile = str_c("./temp_webgestalt/", temp$treatment,".rnk")

   DE_ttT = p_eqtl %>%
    outm10_whole_genome(., control, data = example_skincolor3) %>%
    select(-tfbm_all, -tfbm_immue)

    args =
    DE_ttT %>%
    select(treatment, ttT) %>%
    mutate(file_output = str_c(getwd(), "/", tempfolder))

  complete_tables = args %>%
    pmap(gsea_webgestalt) %>%
    set_names(args$treatment)

  complete_tables %>% saveRDS("./user_wx/skincolor_gsea_0.05_perm10000.rds")
}
# 
# 
# complete_tables = readRDS("/home/xu/ses-1/user_wx/skincolor_gsea_0.05.rds")
#' ### skincolor darkblack gsea pathway
#+ echo=F, eval=T, warning=FALSE, message=FALSE

complete_tables = readRDS("/home/xu/ses-1/user_wx/fgsea_skincolor.rds")
darkblack = 
  complete_tables %>% 
  filter(treatment == "color_byinterviewer3_DarkBlack") %>% 
  pluck("gsea")
  
darkblack[[1]] %>% 
  filter(padj < 0.05) %>%
  select(-leadingEdge) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### skincolor LigthMed gsea pathway
#+ echo=F, eval=T, warning=FALSE, message=FALSE
lightmed = 
  complete_tables %>% 
  filter(treatment == "color_byinterviewer3_LightMed") %>% 
  pluck("gsea")

lightmed[[1]] %>% 
  filter(padj < 0.05) %>%
  select(-leadingEdge) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### tfbm_all (650 TF, corrected for MC using FDR)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

tfbm_all = p_eqtl %>%
  outm10_whole_genome(., control, data = example_skincolor3) %>% 
  filter(ttT!="NULL") %>% 
  select(-tfbm_immue,-ttT) %>% 
  unnest(tfbm_all) %>% 
  mutate_at(.vars = vars(c("p_under", "p_over", "m_uni.p.value", "m_cov.p.value")),
            .funs = list( ~ p.adjust(., method = "fdr"))) %>%
  filter(p_under<0.05|p_over<0.05|m_uni.p.value<0.05|m_cov.p.value<0.05)
tfbm_all %>% kableExtra::kable() %>% kableExtra::kable_styling()

#' ### tfbm_immue (14 TF, corrected for MC using FDR)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
tfbm_immue = p_eqtl %>%
  outm10_whole_genome(., control, data = example_skincolor3) %>% 
  filter(ttT!="NULL") %>% 
  select(-tfbm_all,-ttT) %>% 
  unnest(tfbm_immue) %>% 
  mutate_at(.vars = vars(c("p_under", "p_over", "m_uni.p.value", "m_cov.p.value")),
            .funs = list( ~ p.adjust(., method = "fdr"))) %>%
  filter(p_under<0.05|p_over<0.05|m_uni.p.value<0.05|m_cov.p.value<0.05)

tfbm_immue %>% kableExtra::kable() %>% kableExtra::kable_styling()

# DE  %>% openxlsx::write.xlsx("./user_wx/color3_bespoke_DE_gene.xlsx")


#' ## TFBM for each signature for skin color
#+ echo=F, eval=T, warning=FALSE, message=FALSE


#' ### tfbm_all (650 TF, 14 signatures, corrected for MC using FDR)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
tfbm_all = outm10(p_eqtl, control, data = example_skincolor3) %>% 
  select(p_eqtl, treatment, gene_set_name, tfbm) %>% 
  unnest_wider(tfbm) %>% 
  filter(tfbm_all!= "NULL") %>% 
  unnest(tfbm_all) %>% 
  mutate(m_uni.p.adj = p.adjust(m_uni.p.value,method = "fdr")) %>%
  select(-p_under, -p_over, -m_cov.p.value , - tfbm_immue)

tfbm_all %>% 
  filter(m_uni.p.adj<0.05) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()
#' ### tfbm_immue (14 TF, 14 signatures, corrected for MC)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
tfbm_immue = outm10(p_eqtl, control, data = example_skincolor3) %>% 
  select(p_eqtl, treatment, gene_set_name, tfbm) %>% 
  unnest_wider(tfbm) %>% 
  filter(tfbm_immue!= "NULL") %>% 
  unnest(tfbm_immue) %>% 
  mutate(m_uni.p.adj = p.adjust(m_uni.p.value,method = "fdr")) %>%
  select(-p_under, -p_over, -m_cov.p.value , - tfbm_all)


tfbm_immue %>% 
  filter(m_uni.p.adj<0.05) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ## TFBM enrichment analysis for prior signatures
#+ echo=F, eval=T, warning=FALSE, message=FALSE
example = readRDS("/home/xu/ses-1/user_wx/color3_bespoke_29.03.2021.rds")
example_skincolor3 = example %>% filter(table1 %>% str_detect("aging|whole_genome"))


#' ### tfbm_all (14 signatures, corrected for MC using FDR)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
tfbm_all = outm11(p_eqtl, control, data = example_skincolor3) %>% 
  select(p_eqtl, treatment, gene_set_name, tfbm) %>% 
  unnest_wider(tfbm) %>% 
  filter(tfbm_all!= "NULL") %>% 
  select(-treatment, - tfbm_immue) %>% 
  unnest(tfbm_all) %>% 
  select(-m_uni.p.value, -m_cov.p.value) %>% 
  unique() %>% 
  mutate(
    p_under_adj = p_under %>% p.adjust (method = "fdr"),
    p_over_adj = p_over %>% p.adjust(method = "fdr")
  ) %>% 
  filter(p_under_adj<0.05|p_over_adj<0.05)

tfbm_all %>% 
  unique() %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ### tfbm_immue (14 signatures, corrected for MC using FDR)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
tfbm_immue = outm11(p_eqtl, control, data = example_skincolor3) %>% 
  select(p_eqtl, treatment, gene_set_name, tfbm) %>% 
  unnest_wider(tfbm) %>% 
  filter(tfbm_all!= "NULL") %>% 
  select(-treatment, - tfbm_all) %>% 
  unnest(tfbm_immue) %>% 
  select(-m_uni.p.value, -m_cov.p.value) %>% 
  unique() %>% 
  mutate(
    p_under_adj = p_under %>% p.adjust (method = "fdr"),
    p_over_adj = p_over %>% p.adjust(method = "fdr")
  ) %>% 
  filter(p_under_adj<0.05|p_over_adj<0.05)


tfbm_immue %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

