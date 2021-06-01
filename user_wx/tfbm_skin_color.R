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
example = readRDS("/home/xu/ses-1/user_wx/color3_bespoke_28.03.2021.rds")
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
                Aging_Up = signatures$outcome_set$aging_up_mRNA
                # Aging_Cluster_Complement = signatures$outcome_set$aging_cluster_complement_mRNA
                ),
           ilabels = FALSE, zcolor = "style",snames = " ",
           ellipse = FALSE, opacity = 0.15, ilcs = 1.2,
           box = FALSE,
           sncs = 0.8)

#' #### Geonowide DE gene for Darkblack logFC plotting (663 genes)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

p_eqtl %>% 
  outttT_sig(control, data = example_skincolor3) %>% 
  unnest(ttT) %>% 
  ggplot(aes(x = logFC, fill = as.factor(treatment))) +
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



#' #### Darkblack GOResults biological pathway from GOrilla(Gene Ontology enRIchment anaLysis and visuaLizAtion tool)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

# DE = p_eqtl %>%
# set_names() %>%
# map(~ outm10_whole_genome(., control, data = example_skincolor3) %>%
# mutate(DE_gene = map(ttT, ~ filter(., adj.P.Val<0.05) %>% pull(gene))) %>% filter(map_lgl(DE_gene, ~ length(.)>0)) %>% pull(DE_gene)) %>%
# unlist(recursive = F)
# View(DE)
# DE  %>% openxlsx::write.xlsx("./user_wx/color3_bespoke_DE_gene.xlsx")



knitr::include_graphics("/home/xu/ses-1/user_wx/GO_db.png")

GO_db = read_excel("/home/xu/ses-1/user_wx/GO_db.xlsx")

GO_db %>% 
  select(1:5) %>% 
  filter(`FDR q-value` < 0.05) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()

#' #### LightMed GOResults biological pathway from GOrilla(Gene Ontology enRIchment anaLysis and visuaLizAtion tool)
#+ echo=F, eval=T, warning=FALSE, message=FALSE

knitr::include_graphics("/home/xu/ses-1/user_wx/GO_lm.png")

GO_lm = read_excel("/home/xu/ses-1/user_wx/GO_lm.xlsx")

GO_lm %>% 
  select(1:5) %>% 
  filter(`FDR q-value` < 0.05) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_styling()

#' #### intersection of GOResults biological pathways for darkblack and lightmed
#+ echo=F, eval=T, warning=FALSE, message=FALSE

venn::venn(list(DarkBlack = GO_db$`GO Term`,
                LightMed = GO_lm$`GO Term`),
ilabels = FALSE, zcolor = "style",snames = " ",
ellipse = FALSE, opacity = 0.15, ilcs = 1.2,
box = FALSE,
sncs = 0.8)


DE_ttT = p_eqtl %>%
  outm10_whole_genome(., control, data = example_skincolor3) %>% 
  select(-tfbm_all, -tfbm_immue)



if(reproduciable <- FALSE) {
  rankFile = str_c("./temp_webgestalt/", temp$treatment,".rnk")
  
  args = 
    DE_ttT %>%
    select(treatment, ttT) %>% 
    mutate(file_output = str_c(getwd(), "/", tempfolder))
  
  complete_tables = args %>%
    pmap(gsea_webgestalt) %>%
    set_names(args$treatment)
  
  complete_tables %>% saveRDS("./user_wx/skincolor_gsea.rds")
}


complete_tables = readRDS("/home/xu/ses-1/user_wx/skincolor_gsea.rds")
#' ##### skincolor darkblack gsea pathway
#+ echo=F, eval=T, warning=FALSE, message=FALSE
complete_tables$color_byinterviewer3_DarkBlack %>% 
  filter(FDR < 0.05) %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

#' ##### skincolor LigthMed gsea pathway
#+ echo=F, eval=T, warning=FALSE, message=FALSE
complete_tables$color_byinterviewer3_LightMed %>% 
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

tfbm_immue %>% kableExtra::kable() %>% kableExtra::kable_styling()

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

#' ### TFBM for skin color (loose method, use predefiend signature genes as subset of genes for tfbm)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
example = readRDS("/home/xu/ses-1/user_wx/color3_bespoke_29.03.2021.rds")
example_skincolor3 = example %>% filter(table1 %>% str_detect("aging|whole_genome"))


#' #### tfbm_all (corrected for MC using FDR)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
tfbm_all = outm11(p_eqtl, control, data = example_skincolor3) %>% 
  select(p_eqtl, treatment, gene_set_name, tfbm) %>% 
  unnest_wider(tfbm) %>% 
  filter(tfbm_all!= "NULL") %>%
  mutate_at(.vars = vars(c("tfbm_all", "tfbm_immue")),
            .funs = list(~ map(., ~ mutate_if(., is.numeric, p.adjust, method = "fdr")))) %>%
  mutate(
    `tellis_p_under < 0.05` = tfbm_all %>% map(~ filter(.x, p_under < 0.05) %>% pull("tfbm")),
    `tellis_p_over < 0.05` = tfbm_all %>% map(~ filter(.x, p_over < 0.05) %>% pull("tfbm")),
    `regression p uni < 0.05` = tfbm_all %>% map(~ filter(.x, m_uni.p.value< 0.05) %>% pull("tfbm")),
    `regression p cov < 0.05` = tfbm_all %>% map(~ filter(.x, m_cov.p.value< 0.05) %>% pull("tfbm"))
  ) %>% 
  select(-tfbm_all, - tfbm_immue)

tfbm_all %>% 
  # filter(map_lgl(`regression p uni < 0.05`, ~ length(.)>0)) %>% 
  filter(map_lgl(`regression p uni < 0.05`, ~ length(.)>0) | map_lgl(`regression p cov < 0.05`, ~ length(.)>0) |
           map_lgl(`tellis_p_under < 0.05`, ~ length(.)>0) | map_lgl(`tellis_p_over < 0.05`, ~ length(.)>0)) %>% 
  mutate(treatment = treatment %>% str_remove("color_byinterviewer3_"))
# %>% 
#   kableExtra::kable() %>%
#   kableExtra::kable_styling()
#' #### tfbm_immue (not corrected for MC)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
tfbm_immue = outm11(p_eqtl, control, data = example_skincolor3) %>% 
  select(p_eqtl, treatment, gene_set_name, tfbm) %>% 
  unnest_wider(tfbm) %>% 
  filter(tfbm_all!= "NULL") %>% 
  mutate(
    `tellis_p_under < 0.05` = tfbm_immue %>% map(~ filter(.x, p_under < 0.05) %>% pull("tfbm")),
    `tellis_p_over < 0.05` = tfbm_immue %>% map(~ filter(.x, p_over < 0.05) %>% pull("tfbm")),
    `regression p uni < 0.05` = tfbm_immue %>% map(~ filter(.x, m_uni.p.value< 0.05) %>% pull("tfbm")),
    `regression p cov < 0.05` = tfbm_immue %>% map(~ filter(.x, m_cov.p.value< 0.05) %>% pull("tfbm"))
  ) %>% 
  select(-tfbm_all, - tfbm_immue)


tfbm_immue %>% 
  filter(map_lgl(`regression p uni < 0.05`, ~ length(.)>0) | map_lgl(`regression p cov < 0.05`, ~ length(.)>0) |
           map_lgl(`tellis_p_under < 0.05`, ~ length(.)>0) | map_lgl(`tellis_p_over < 0.05`, ~ length(.)>0)) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()
