#' ---
#' title: DE and tfbm for skin color in race strata nonhispanic black
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
load_data(reconciled = FALSE, remove_inflam = FALSE)
source("/home/xu/ses-1/user_wx/extract_v2.R")
# example_race = readRDS("~/ses-1/user_wx/race_bespoke_15.03.2021.rds")
example = readRDS("/home/xu/ses-1/user_wx/color3_bespoke_NonHblack_strata_28.03.2021.rds")
example_skincolor3 = example %>% filter(table1 %>% str_detect("aging|whole_genome"))
control = "ancestryPC_ses"
p_eqtl = 0.01
tempfolder = "temp_webgestalt"
#' ### DE and TFBM genowide for skin color
#+ echo=F, eval=T, warning=FALSE, message=FALSE
# DE_gene = p_eqtl %>%
#   set_names() %>% 
#   map(~ outm10_whole_genome(., control, data = example_skincolor3) %>% 
#         mutate(DE_gene = map(ttT, ~ filter(., adj.P.Val< 0.05) %>% pull(gene))) %>%
#         filter(map_lgl(DE_gene, ~ length(.)>0)) %>%
#         pull(DE_gene)) %>% 
#   unlist(recursive = F)



#' #### DE gene for Darkblack 
#+ echo=F, eval=T, warning=FALSE, message=FALSE

temp = p_eqtl %>%
  outm10_whole_genome(., control, data = example_skincolor3) %>% 
  filter(ttT!="NULL") %>% 
  mutate(DE_gene = map(ttT, ~ filter(., adj.P.Val< 0.05) %>% pull(gene))) %>%
  filter(map_lgl(DE_gene, ~ length(.)>0))
DE_gene = list2(!!str_c(temp$treatment, temp$p_eqtl) := temp %>% pull(DE_gene)) %>% 
  unlist(recursive = F)

venn::venn(list(DE_gene = DE_gene$color_byinterviewer3_DarkBlack0.01,
                Aging_Down = signatures$outcome_set$aging_down_mRNA,
                Aging_Up = signatures$outcome_set$aging_up_mRNA,
                Aging_Cluster_Complement = signatures$outcome_set$aging_cluster_complement_mRNA),
           ilabels = FALSE, zcolor = "style",snames = " ",
           ellipse = FALSE, opacity = 0.15, ilcs = 1.2,
           box = FALSE,
           sncs = 0.8)

DE_ttT = p_eqtl %>%
  outm10_whole_genome(., control, data = example_skincolor3) %>% 
  select(-tfbm_all, -tfbm_immue) %>% 
  filter(ttT!="NULL") 



if(reproduciable <- FALSE) {
  rankFile = str_c("./temp_webgestalt/", DE_ttT$treatment,"_in_NonhispanicBlack_stratum.rnk")
  
  args = 
    DE_ttT %>%
    select(treatment, ttT) %>% 
    mutate(file_output = str_c(getwd(), "/", tempfolder))
  
  complete_tables = args %>%
    pmap(gsea_webgestalt) %>%
    set_names(args$treatment)
  
  complete_tables %>% saveRDS("./user_wx/skincolor_nonhispanicblack_strata_gsea.rds")
}


complete_tables = readRDS("/home/xu/ses-1/user_wx/skincolor_nonhispanicblack_strata_gsea.rds")
#' ##### skincolor darkblack gsea pathway
#+ echo=F, eval=T, warning=FALSE, message=FALSE
complete_tables$color_byinterviewer3_DarkBlack %>% 
  filter(FDR < 0.05) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()


#' #### tfbm_all (corrected for MC using FDR)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

tfbm_all = p_eqtl %>%
  outm10_whole_genome(., control, data = example_skincolor3) %>% 
  filter(ttT!="NULL") %>% 
  mutate_at(.vars = vars(c("tfbm_all", "tfbm_immue")),
            .funs = list(~ map(., ~ mutate_if(., is.numeric, p.adjust, method = "fdr")))) %>%
  mutate(`tellis_p_under < 0.05` = tfbm_all %>% map(~ filter(.x, p_under < 0.05) %>% pull("tfbm")),
         `tellis_p_over < 0.05` = tfbm_all %>% map(~ filter(.x, p_over < 0.05) %>% pull("tfbm")),
         `regression p uni < 0.05` = tfbm_all %>% map(~ filter(.x, m_uni.p.value< 0.05) %>% pull("tfbm")),
         `regression p cov < 0.05` = tfbm_all %>% map(~ filter(.x, m_cov.p.value< 0.05) %>% pull("tfbm"))) %>% 
  select(-tfbm_all, - tfbm_immue, -ttT)

tfbm_all %>% kableExtra::kable() %>% kableExtra::kable_styling()

#' #### tfbm_immue (not corrected for MC)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
tfbm_immue = p_eqtl %>%
  outm10_whole_genome(., control, data = example_skincolor3) %>%
  filter(ttT!="NULL") %>% 
  mutate(`tellis_p_under < 0.05` = tfbm_immue %>% map(~ filter(.x, p_under < 0.05) %>% pull("tfbm")),
         `tellis_p_over < 0.05` = tfbm_immue %>% map(~ filter(.x, p_over < 0.05) %>% pull("tfbm")),
         `regression p uni < 0.05` = tfbm_immue %>% map(~ filter(.x, m_uni.p.value< 0.05) %>% pull("tfbm")),
         `regression p cov < 0.05` = tfbm_immue %>% map(~ filter(.x, m_cov.p.value< 0.05) %>% pull("tfbm"))) %>% 
  select(-tfbm_all, - tfbm_immue, -ttT)

tfbm_immue %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()

# DE  %>% openxlsx::write.xlsx("./user_wx/color3_bespoke_DE_gene.xlsx")


#' ### TFBM for each signature for skin color
#+ echo=F, eval=T, warning=FALSE, message=FALSE

temp = p_eqtl %>%
  outm10(., control, data = example_skincolor3) %>% 
  filter(ttT!= "NULL") %>% 
  mutate(DE_gene = map(ttT, ~ filter(., adj.P.Val< 0.05) %>% pull(gene))) %>%
  filter(map_lgl(DE_gene, ~ length(.)>0))
# DE_gene = list2(!!str_c(temp$treatment, temp$p_eqtl) := temp %>% pull(DE_gene)) %>% 
#   unlist(recursive = F)
# 
# venn::venn(list(DE_gene = DE_gene$color_byinterviewer3_DarkBlack0.01,
#                 Aging_Down = signatures$outcome_set$aging_down_mRNA,
#                 Aging_Up = signatures$outcome_set$aging_up_mRNA,
#                 Aging_Cluster_Complement = signatures$outcome_set$aging_cluster_complement_mRNA),
#            ilabels = FALSE, zcolor = "style",snames = " ",
#            ellipse = FALSE, opacity = 0.15, ilcs = 1.2,
#            box = FALSE,
#            sncs = 0.8)

#' #### tfbm_all (corrected for MC using FDR)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
tfbm_all = outm10(p_eqtl, control, data = example_skincolor3) %>% 
  select(p_eqtl, treatment, gene_set_name, tfbm) %>% 
  unnest_wider(tfbm) %>% 
  filter(tfbm_all!= "NULL") %>% 
  # mutate_at(.vars = vars(c("tfbm_all", "tfbm_immue")),
  #           .funs = list(adj.p = ~ .x)) %>% 
  # mutate_at(.vars = vars(c("tfbm_all_adj.p", "tfbm_immue_adj.p")),
  #           .funs = list(~ map(., ~ mutate_if(., is.numeric, p.adjust, method = "fdr")))) %>%
  # mutate_at(.vars = vars(c("tfbm_all_adj.p", "tfbm_immue_adj.p")),
  #           .funs = list(~ map(., ~ mutate(., ind = pmin(p_under, p_over, m_uni.p.value, m_cov.p.value, na.rm=T))))) %>% 
  # mutate_at(.vars = vars(c("tfbm_all_adj.p", "tfbm_immue_adj.p")),
  #           .funs = list(~ map(., ~ filter(., ind <0.05)))) %>% 
  # filter(map_lgl(tfbm_all_adj.p, ~ dim(.)[1] != 0))
  
  mutate_at(.vars = vars(c("tfbm_all", "tfbm_immue")),
            .funs = list(~ map(., ~ mutate_if(., is.numeric, p.adjust, method = "fdr")))) %>%
  mutate(
    # `tellis_p_under < 0.05` = tfbm_all %>% map(~ filter(.x, p_under < 0.05) %>% pull("tfbm")),
    #     `tellis_p_over < 0.05` = tfbm_all %>% map(~ filter(.x, p_over < 0.05) %>% pull("tfbm")),
    `regression p uni < 0.05` = tfbm_all %>% map(~ filter(.x, m_uni.p.value< 0.05) %>% pull("tfbm")),
    # `regression p cov < 0.05` = tfbm_all %>% map(~ filter(.x, m_cov.p.value< 0.05) %>% pull("tfbm"))
  ) %>% 
  select(-tfbm_all, - tfbm_immue)

tfbm_all %>%
  filter(map_lgl(`regression p uni < 0.05`, ~ length(.)>0)) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' #### tfbm_immue (not corrected for MC)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
tfbm_immue = outm10(p_eqtl, control, data = example_skincolor3) %>% 
  select(p_eqtl, treatment, gene_set_name, tfbm) %>% 
  unnest_wider(tfbm) %>% 
  filter(tfbm_all!= "NULL") %>% 
  # mutate_at(.vars = vars(c("tfbm_all", "tfbm_immue")),
  #           .funs = list(adj.p = ~ .x)) %>% 
  # mutate_at(.vars = vars(c("tfbm_all_adj.p", "tfbm_immue_adj.p")),
  #           .funs = list(~ map(., ~ mutate_if(., is.numeric, p.adjust, method = "fdr")))) %>%
  # mutate_at(.vars = vars(c("tfbm_all_adj.p", "tfbm_immue_adj.p")),
  #           .funs = list(~ map(., ~ mutate(., ind = pmin(p_under, p_over, m_uni.p.value, m_cov.p.value, na.rm=T))))) %>% 
  # mutate_at(.vars = vars(c("tfbm_all_adj.p", "tfbm_immue_adj.p")),
  #           .funs = list(~ map(., ~ filter(., ind <0.05)))) %>% 
  # filter(map_lgl(tfbm_all_adj.p, ~ dim(.)[1] != 0))
  # 
  # mutate_at(.vars = vars(c("tfbm_all", "tfbm_immue")),
#           .funs = list(~ map(., ~ mutate_if(., is.numeric, p.adjust, method = "fdr")))) %>%
mutate(
  # `tellis_p_under < 0.05` = tfbm_all %>% map(~ filter(.x, p_under < 0.05) %>% pull("tfbm")),
  #     `tellis_p_over < 0.05` = tfbm_all %>% map(~ filter(.x, p_over < 0.05) %>% pull("tfbm")),
  `regression p uni < 0.05` = tfbm_immue %>% map(~ filter(.x, m_uni.p.value< 0.05) %>% pull("tfbm")),
  # `regression p cov < 0.05` = tfbm_all %>% map(~ filter(.x, m_cov.p.value< 0.05) %>% pull("tfbm"))
) %>% 
  select(-tfbm_all, - tfbm_immue)

tfbm_immue %>%
  filter(map_lgl(`regression p uni < 0.05`, ~ length(.)>0)) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()



