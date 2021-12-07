library(tidyverse)
p_eqtl = 0.05
bespoke_m7nnfun = function(example) {
  example = example%>% 
    hoist(out, out1 = list("result", "example0")) %>%
    select(out1) %>% 
    unnest(out1) %>%
    hoist(out, p = list("result", "m7_nn", 1, "p")) %>% 
    unnest_longer(p) %>% 
    group_by(treatment) %>% 
    arrange(p) %>% 
    mutate(p.adj = p.adjust(p, method = "fdr")) %>% 
    select(-out) %>% 
    filter(p.adj <0.05)
  
}
example <- readRDS("~/ses-1/user_wx/m7nn_with1k_aging_sc5levels_allsample_bespoke.rds")
with = bespoke_m7nnfun(example)
ex0_m7pca =
  with %>% 
group_by(treatment, gene_set_name) %>% 
  mutate(noPC = n()) %>% 
  mutate(noPC = ifelse(noPC>1, noPC, NA),
         pval2 = case_when(
           # p.adj<0.0001 ~ 100000,
           #                 p.adj<0.001 ~ 25000,
                           p.adj<0.01 ~ 10000,
                           p.adj<0.05 ~ 5000,
                           p.adj>0.05 ~ 0.0000001)) %>% 
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


ggplot(ex0_m7pca, aes(treatment, gene_set_name, size = pval2
                      # fill = `1KI Genes`,
                      # colour = `1KI Genes`
)) +
  geom_point(stroke = 1.5, shape = 21, alpha = 0.4, colour = "darkblue", fill = "navy") +
  geom_text(aes(label=noPC), show.legend=FALSE, size=3, color ="white") +
  # scale_fill_manual(values = c("darkblue", "goldenrod3")) +
  # scale_color_manual(values = c("darkblue", "goldenrod3")) +
  # geom_jitter(height = 0.00000025) +
  # gghighlight(class == "inflam") +
  theme_bw() +
  labs(
    caption = paste("p_eqtl = ", p_eqtl),
    title = "Figure . PCA regression",
    y = "mRNA Signatures PC",
    x = "Race Skincolor") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        # plot.margin=unit(c(1, 1, 0.1, 1), "cm"),
        # axis.text.y = element_text(color=axiscolor),
        text = element_text(size=10, face = "bold")) +
  scale_size_continuous(name = "Adjusted P-value", range = c(0, 21), 
                        limits = c(0.0000001, 10000), breaks = c(0.0000001, 5000, 10000),
                        labels = c("n.s.", "p<0.05", "p<0.01")) +
  scale_alpha(guide = 'none') +
  guides(shape = guide_legend(override.aes = list(size = 10)),
         fill = guide_legend(override.aes = list(size = 8)))   
