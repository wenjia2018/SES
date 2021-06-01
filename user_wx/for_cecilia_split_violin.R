example = readRDS("~/ses-1/user_wx/skincolor_eqtl005_mediation_singel_gene_aging_11.05.2021.rds")
example = readRDS("~/ses-1/user_wx/skincolor3_wholegenome_aging_composite_ancestry_14.05.2021.rds")
library(stringi)
example = example$out[[1]]$result$example0

temp =
  example %>%
  hoist(out, ttT = list("result", "m8_fdr", 1, "mediation_single")) %>% 
  unnest_longer(ttT) %>% 
  unnest_longer(ttT) %>% 
  hoist(ttT, med_single_prop = list("result", "other", "med_prop")) %>%
  hoist(ttT, med_single_beta = list("result", "other", "med_ACME")) %>%
  hoist(ttT, med_single_ade = list("result", "other", "med_ADE")) %>% 
  dplyr::select(-controls, -ttT, -out) %>% 
  dplyr::mutate(mediator=stri_extract(ttT_id, regex='[^_]*'),
                gene1=stri_extract_last(ttT_id, regex='([^_]+$)')) %>% 
  mutate(treatment= case_when(treatment =="color_byinterviewer3_DarkBlack" ~ "Dark Black",
                              treatment =="color_byinterviewer3_LightMed" ~ "Light Medium")) %>% 
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
                                  gene_set_name =="Aging Cluster Complement" ~ "Aging Complement",
                                  gene_set_name =="Aging" ~ "Aging All Genes",
                                  TRUE ~ gene_set_name))


temp2<-temp %>%
  mutate(med_single_prop= as.numeric(med_single_prop),
         med_single_beta= as.numeric( med_single_beta),
         med_single_ade =as.numeric(med_single_ade)
         
  ) %>% 
  group_by(treatment, gene_set_name, mediator) %>%
  # add_tally()
  mutate(across(.cols = "gene1",
                .fns = ~if_else(. <= 0,"Low","High"),
                .names = "{.col}_"))

temp3 <-temp2 %>% 
  group_by(treatment, gene_set_name, mediator  ) %>%
  nest()

exx<-temp3 %>% 
  ungroup %>% 
  mutate(ttest=map(.x=temp3$data, .f=~ wilcox.test(.x$med_single_prop)))

exx <- exx %>%
  mutate(tidy = map(ttest, broom::tidy))

exx <-exx %>%
  unnest(tidy, .drop = T) %>% 
  dplyr::select(treatment, gene_set_name, mediator, data, p.value) %>% 
  mutate(padj.tval=p.adjust(p.value, method="BH"))

#Filtering only those with significant p-values in the t.test
exx %>% mutate(padj.tval=p.adjust(p.value, method="BH"))  %>% dplyr::filter(padj.tval<0.05) %>% print(n=Inf)

ex2hist <-exx %>% unnest(data) #%>% dplyr:: filter(gene_set_name== c("Aging All Genes"))

ex2hist <-ex2hist %>% mutate(pval=case_when(padj.tval<0.05 ~ 0.05,
                                            padj.tval>0.05 ~ 100),
                             pval2=case_when(padj.tval<0.05 ~ "p<0.05",
                                             padj.tval>0.05 ~ "n.s,."))

#'## Mediator smoking

ex2hist2<-ex2hist %>% dplyr::filter(mediator=="currentsmoke") %>% 
  gather(key=effect, value=values, 5:6)
# source the file where you save the geom_split_violin
source("./user_wx/geom_split_violin.R")
ex2hist2 %>% 
  ggplot(aes(x = gene_set_name, y = values, fill = effect)) +
  geom_split_violin()+
  facet_wrap( ~ treatment, scales = "free_x")  +
  # geom_hline(yintercept = 0, color = "red", linetype = "dashed")+
  scale_fill_manual(values=c("goldenrod3", "lightblue"))+
  coord_flip()
