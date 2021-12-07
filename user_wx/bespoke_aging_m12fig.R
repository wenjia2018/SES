bespoke_m12fun = function(example0_with1k) {
  example = example0_with1k %>% 
    hoist(out, out1 = list("result", "example0")) %>%
    select(out1) %>% 
    unnest(out1) %>%
    hoist(out, m = list("result", "m12_fdr", 1, "other", "m")) %>% 
    filter(!map_lgl(m, ~is.null(.x))) %>% 
    mutate(p = m %>% map(~ dplyr::slice(.x, which.min(adj.p.withinunion))) %>% map_dbl(~ .x %>% pluck("adj.p.withinunion"))) %>% 
    arrange(p) %>% 
    dplyr::select(treatment, gene_set_name, controls, p) %>% 
    filter(p<0.05)
}

# all sample

example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_aging_sc5levels_allsample_bespoke.rds")
example0_with1k <- readRDS("~/ses-1/user_wx/m12_with1k_aging_sccont_allsample_bespoke.rds")

res = bespoke_m12fun(example0_with1k) %>% select(-controls) %>% 
  mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
         gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
         # gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI"),
         # gene_set_name = gene_set_name %>% replace(gene_set_name == "Ctra", "CTRA"),
         # gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflame", "Inflamation in CTRA"),
         gene_set_name= case_when(gene_set_name =="Aging Up Cl1" ~ "Innate/Adaptive Immunity",
                                  gene_set_name =="Aging Up Cl2" ~ "Actin Cytoskeleton",#/Focal Adhesion/Tight Junctions",
                                  gene_set_name =="Aging Up Cl3" ~ "Fatty Acid Metabolism",#/Peroxisome Activity",
                                  gene_set_name =="Aging Up Cl4" ~ "Lysosome Metabolism",#/Glycosaminoglycan Degradation",
                                  # gene_set_name =="Aging Down Cl1" ~ "Aging Down Cl1",
                                  gene_set_name =="Aging Down Cl1a" ~ "RNA Metabolism",#/Ribosome Biogenesis/Purine Metabolism",
                                  gene_set_name =="Aging Down Cl1b" ~ "Mitochondrial",
                                  gene_set_name =="Aging Down Cl1c" ~ "DNA Replication",#/Elongation/Mismatch Repair",
                                  gene_set_name =="Aging Down Cl2" ~ "Ribosome",
                                  gene_set_name =="Aging Down Cl3" ~ "Immune Genes",
                                  gene_set_name =="Aging" ~ "Aging Composite",
                                  # gene_set_name == "Aging Cluster Complement" ~ "Aging Cluster Complement",
                                  TRUE ~ gene_set_name),
         gene_set_name = factor(gene_set_name, levels = c(
           
           "Aging Down Cl1",
           "RNA Metabolism",#/Ribosome Biogenesis/Purine Metabolism",
           "Mitochondrial",
           "DNA Replication",#/Elongation/Mismatch Repair",
           "Ribosome",
           "Immune Genes",
           "Innate/Adaptive Immunity", 
           "Actin Cytoskeleton",#/Focal Adhesion/Tight Junctions",
           "Fatty Acid Metabolism",#/Peroxisome Activity", 
           "Lysosome Metabolism",#/Glycosaminoglycan Degradation"
           "Aging Up",
           "Aging Down",
           "Aging Cluster Complement",
           "Aging Composite"
         )) %>% fct_rev
         # gene_set_name = str_c(gene_set_name, "_", p_id)
  )


res  %>%  kableExtra::kable() %>% kableExtra::kable_styling()

res$pval2 = -log10(res$p)
res$pval2[which(res$pval2>-log10(0.0001))] = -log10(0.0001)
res$treatment = as.character(res$treatment)
res$gene_set_name = as.character(res$gene_set_name)




ggplot(res, aes(x = treatment,y  = gene_set_name,  size = pval2)) +
  geom_point(alpha  = 0.5,shape = 21, stroke = 1) +
  theme_bw(base_size = 12) +
  scale_size_continuous(range = c(1,12),
                        name = "Adjusted\n p-value", 
                        limits = c(-log10(0.05), -log10(0.0001)),
                        breaks = c(-log10(0.05),-log10(0.01),-log10(0.001), -log10(0.0001)), 
                        labels = c("p<0.05", "p<0.01","p<0.001","p<0.0001")) + 
  xlab("SES Indicators") + ylab("mRNA Signatures") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title = element_text(size = 11, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=10, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=10,face = "bold" ,family = "Calibri", angle = 30,hjust = 1)) +
  theme(legend.text=element_text(size=10, family = "Calibri")) +
  theme(legend.title = element_text(size = 10, face = "bold", family = "Calibri")) +
  theme(legend.position="right") +
  scale_y_discrete(limits = rev)  +
  guides(color = guide_legend(override.aes = list(size = 7))) +
  guides(fill = guide_legend(override.aes = list(size = 7))) +
  guides(shape = guide_legend(override.aes = list(size = 7))) 
