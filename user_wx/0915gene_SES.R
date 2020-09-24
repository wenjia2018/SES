# estimate pearson correlations (list wise) for each gene and SES.
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
library(here)
library(dbr)
library(venn)
walk(dir(path = here("R"),full.names = TRUE), source)
############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################

load_data(reconciled = FALSE, remove_inflam = FALSE)
define_treatments_and_controls()
recode_variables_in_dat()
treatment = "ses_sss_composite"
  
  # Specify whole-genome regression of rna on design
  y <- dat %>% Biobase::exprs() %>% t %>% as.data.frame() %>%  rownames_to_column(var = "AID") %>%
    mutate(AID = substring(AID, 2, 100) %>% as.character())
  X <- dat %>% phenoData() %>% .@data %>% dplyr::select(AID, all_of(treatment))
  gene = y %>% dplyr::select(-1) %>% as.list()

  pearson_test <-  gene %>%
    map(~ cor.test(.x, X[,2]))
  
  a = pearson_test %>%
    map(~ tibble::tibble(cor = purrr::pluck(.x, "estimate", .default = NA), 
                        p = purrr::pluck(.x, "p.value", .default = NA))) %>%
    enframe(name = "gene") %>%
    unnest("value") %>% 
    mutate(adj.p = p %>% p.adjust(method = "fdr"))

  a %>%
    arrange(adj.p, decreasing = T) %>%
    mutate_at(.vars = vars(c("p","adj.p")),
              .funs = list(. %>% format(digits = 3, scientific =T))) %>% 
    openxlsx::write.xlsx("./user_wx/gene_ses4_pearsoncor.xlsx")

  
  ses4_unique_down <- read_excel("user_wx/DE_removetable1signatures_ses4income_fulllist.xlsx",
                                                              sheet = "ses4_unique_down",
                                 col_names = "gene")
  income_unique_down <- read_excel("user_wx/DE_removetable1signatures_ses4income_fulllist.xlsx", 
                                      sheet = "income_unique_down", col_names = "gene")
  ses4_income_intersection_down <- read_excel("user_wx/DE_removetable1signatures_ses4income_fulllist.xlsx", 
                                                  sheet = "ses4_income_intersection_down", 
                                                 col_names = "gene")
  

  ses4_unique_down_cor = ses4_unique_down %>% left_join(a)
  ses4_unique_down_cor %>% openxlsx::write.xlsx("./user_wx/gene_ses4uniquedown_pearsoncor.xlsx")

  
  # Divide SES composite into quintiles.
# Report the mean abundance level of each gene at each quintile level. 
 
   X_quintiles = X %>% ntile(5) %>% as.factor()
   data = tibble(group = X_quintiles, y %>% dplyr::select(-1)) 
   gene_mean_bygroup = data %>%
     group_by(group) %>%
     summarise_all( ~ mean(., na.rm=TRUE))
   
   temp = t(gene_mean_bygroup) %>% as.data.frame() %>% rownames_to_column( var = "gene")
   
   plot.ts(gene_mean_bygroup[1:5,2:11])
   # get row and colnames in order
   # temp <-data.table:: transpose(gene_mean_bygroup)
   # colnames(temp) <- rownames(gene_mean_bygroup)
   # rownames(temp) <- colnames(gene_mean_bygroup)
   # temp = rownames_to_column(temp, var = "gene")
   temp[1,1] = ""
   
   
   temp %>% openxlsx::write.xlsx("./user_wx/gene_mean_ses4group.xlsx")
   ses4_unique_down_groupmean = ses4_unique_down %>% left_join(temp)
   
   ses4_unique_down_groupmean %>% openxlsx::write.xlsx("./user_wx/ses4uniquedown_mean_ses4group.xlsx")
   
   