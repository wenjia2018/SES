
#' ---
#' title: Omnibus test results for ses using RLE normalization
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->


#' ###  regression of gene on SES indicator with basic controls
#' * gene ~ ses + controls
#' * about 12,000 genes(differs in each regression)
#' * p values are genowide FDR corrected
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
library(tidyverse)
example_RLE_DE <- readRDS("~/ses-1/user_wx/example_RLE_DE_v2.rds")
data = example_RLE_DE %>% 
  hoist(out, ttT = list("result", "ttT")) %>% 
  filter(map_lgl(.$ttT, ~dim(.)[1]!=0)) %>% 
  mutate(gene_sig = map(ttT, ~ filter(., adj.P.Val<0.05) %>% pull(gene)) %>% set_names(treatment)) %>% 
  # mutate(min_p = map(ttT, ~ dplyr::slice_min(.,adj.P.Val) %>% pull(adj.P.Val)) %>% set_names(treatment)) %>% 
  select(-out, -gene_set_name) 

#' ###  intersection of DE genes for each SES indicator
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
venn::venn(data$gene_sig, zcolor = "style",
           ellipse = FALSE, opacity = 0.15, ilcs = 1,
           box = FALSE,
           sncs = 1)



UpSetR::upset(fromList(data$gene_sig), order.by = "freq")


outcome_set = readRDS("/home/share/preprocessing/preprocessed_two_batches/allpossiblegene.rds")

outcome_set = outcome_set$outcome_set[table1]



temp = NULL
for(i in 1: (data$treatment %>% length)) {
    temp[[i]] = outcome_set %>% map(~ .x %in% data$gene_sig[[i]] %>% sum) %>% bind_rows
}
temp = temp %>% bind_rows() %>% t %>% `colnames<-`(data$treatment)


outcome_set_noinflame = outcome_set %>% map(~ setdiff(.x, outcome_set$inflam1k_mRNA))

temp_noinflame = NULL
for(i in 1: (data$treatment %>% length)) {
  temp_noinflame[[i]] = outcome_set_noinflame %>% map(~ .x %in% data$gene_sig[[i]] %>% sum) %>% bind_rows
}
temp_noinflame = temp_noinflame %>% bind_rows() %>% t %>% `colnames<-`(data$treatment)

#' ###  Number DE genes in each predefined signature for each SES indicator
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA


temp %>% kableExtra::kable() %>% kableExtra::kable_styling()

#' ### Number DE genes in each predefined signature(excluding inflammation 1k) for each SES indicator
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA



temp_noinflame %>% kableExtra::kable() %>% kableExtra::kable_styling()

#' ### Analog to fig1 panleA in pnas paper 
#' * BUT here the size of bubbles means the number of significant genes in this signature
#' * without cell type controls
#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA


temp %>% 
  as.data.frame() %>% 
  rownames_to_column("signature") %>% 
  pivot_longer(cols = (2:6), names_to = "treatment", values_to ="No.SigGene") %>% 
  mutate(No.SigGene = ifelse(No.SigGene==0, NA, No.SigGene),
         `1KIgene` = "with 1KI") %>% 
  bind_rows(temp_noinflame %>% 
              as.data.frame() %>% 
              rownames_to_column("signature") %>% 
              pivot_longer(cols = (2:6), names_to = "treatment", values_to ="No.SigGene") %>% 
              mutate(No.SigGene = ifelse(No.SigGene==0, NA, No.SigGene),
                     `1KIgene` = "without 1KI")) %>% 
  ggplot(mapping = aes(x = treatment, size = No.SigGene, y = signature, fill = `1KIgene`)) + 
  geom_point( shape = 21, alpha = 0.4) +
  scale_fill_manual(values = c("darkblue", "goldenrod3")) +
  scale_color_manual(values = c("darkblue", "goldenrod3")) +
  scale_size_continuous(breaks = c(10, 20, 50, 100)) 

