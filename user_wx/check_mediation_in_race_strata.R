black <- readRDS("/home/xu/ses-1/user_wx/color3_bespoke_NonHblack_strata_16.03.2021.rds")
white <- readRDS("/home/xu/ses-1/user_wx/color3_bespoke_NonHwhite_strata_15.03.2021.rds")
hispanic <- readRDS("/home/xu/ses-1/user_wx/color3_bespoke_Hispanic_strata_15.03.2021.rds")


data_black <- outm7med(p = p_eqtl, mediators, control = "ancestryPC_ses",  black) %>% mutate(group = "Non Hispanic Black")
data_white <- outm7med(p = p_eqtl, mediators, control = "ancestryPC_ses", white) %>% mutate(group = "Non Hispanic White")
# data_hispanic <- outm7med(p = p_eqtl, mediators, control = "ancestryPC_ses", hispanic) %>% mutate(group = "Hispanic")
ex0_m7pca <- data_black %>%
  rbind(data_white) %>%
  # rbind(data_hispanic) %>%
  filter(gene_set_name %>% str_detect("aging")) %>%
  filter(gene_set_name != "aging_down_cl1_mRNA") %>%
  dplyr::select(treatment, gene_set_name, p_med, med_id, mediator) %>% 
  mutate_at(.vars = vars(c("p_med")), .funs = funs(. %>%  str_remove("<") %>% as.numeric)) %>% 
  group_by(treatment, gene_set_name) %>%
  mutate(n=n(), p = min(p_med*n, 1), p_id = med_id) %>% 
  ungroup %>% 
  filter(p <0.05) %>% 
  mutate(pval=case_when(p<0.0001 ~ 0.0001,
                        p<0.001 ~ 0.001,
                        p<0.01 ~ 0.01,
                        p<0.05 ~ 0.05,
                        p>0.05 ~ 100),
         pval2=case_when(p<0.0001 ~ 100000,
                         p<0.001 ~ 25000,
                         p<0.01 ~ 15000,
                         p<0.05 ~ 10000,
                         p>0.05 ~ 0.0000001),
         treatment= case_when(treatment == "raceethnicity_Hispanic" ~ "Hispanic",
                              treatment =="raceethnicity_NonHblack" ~  "Non Hispanic Black",
                              treatment =="color_byinterviewer3_LightMed"  ~ "Light Medium",  
                              treatment =="color_byinterviewer3_DarkBlack"  ~ "Dark Black")) %>%
  mutate(treatment = factor(treatment, levels = c("Hispanic", "Non Hispanic Black","Light Medium","Dark Black"))) %>%
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
                                  gene_set_name == "Inflam1k" ~ "1KI",
                                  gene_set_name == "Ctra" ~ "CTRA",
                                  gene_set_name == "Inflame" ~ "Inflamation in CTRA",
                                  TRUE ~ gene_set_name),
         gene_set_name = factor(gene_set_name, levels = c("1KI", "CTRA", "Antbintf",  "Interferon", "Inflamation in CTRA",
                                                          
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
                                                          "Aging Cluster Complement",
                                                          "Aging Up",
                                                          "Aging Down",
                                                          "Aging Composite"
         )) %>% fct_rev
         # gene_set_name = str_c(gene_set_name, "_", p_id)
  )
