
omnibus_plot = function(data, p_eqtl) {

  skincolor = m8_present(p = p_eqtl, control = "ancestryPC_ses", data)
  ex0_m8fdr = skincolor %>% 
    filter(gene_set_name %>% str_detect("aging")) %>% 
    filter(gene_set_name!="aging_down_cl1_mRNA") %>% 
    dplyr::select(treatment, gene_set_name, p_omnibus) %>% 
    filter(p_omnibus!="NA") %>% 
    rename(p = p_omnibus) %>% 
    mutate(
      pval2=case_when(
        # p<0.0001 ~ 100000,
        p<0.001 ~ 55000,
        p<0.01 ~ 25000,
        p<0.05 ~ 10000,
        p>0.05 ~ 0.0000001),
      treatment= case_when(treatment == "raceethnicity_Hispanic" ~ "Hispanic",
                           treatment =="raceethnicity_NonHblack" ~  "Non Hispanic Black",
                           treatment =="color_byinterviewer3_LightMed"  ~ "Light Medium",  
                           treatment =="color_byinterviewer3_DarkBlack"  ~ "Dark Black")
    )  %>%
    mutate(treatment = factor(treatment, levels = c("Hispanic", "Non Hispanic Black","Light Medium","Dark Black"
    ))) %>%
    mutate(gene_set_name = gene_set_name %>% str_replace_all("_mRNA",""),
           gene_set_name = gene_set_name %>% str_replace_all("_"," ") %>% str_to_title(),
           # gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflam1k", "1KI"),
           # gene_set_name = gene_set_name %>% replace(gene_set_name == "Ctra", "CTRA"),
           # gene_set_name = gene_set_name %>% replace(gene_set_name == "Inflame", "Inflamation in CTRA"),
           gene_set_name= case_when(gene_set_name =="Aging Up Cl1" ~ "Innate/Adaptive Immunity",
                                    gene_set_name =="Aging Up Cl2" ~ "Actin Regulation",#/Focal Adhesion/Tight Junctions",
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
                                    TRUE ~ gene_set_name
           )) %>% 
    mutate(
           pval2 = ifelse(p<0.05,pval2,NA))
  # need to specify shape, because some shape does not have borders. 
  # reorder y axis
  ex0_m8fdr$gene_set_name <- factor(ex0_m8fdr$gene_set_name, levels = c("1KI", "CTRA", "Antbintf",  "Interferon", "Inflamation in CTRA",
                                                                        "Aging Down Cl1",
                                                                        "Lysosome Metabolism",#/Glycosaminoglycan Degradation"
                                                                        "Fatty Acid Metabolism",#/Peroxisome Activity", 
                                                                        "Actin Regulation",#/Focal Adhesion/Tight Junctions",
                                                                        "Ribosome",
                                                                        "DNA Replication",#/Elongation/Mismatch Repair",
                                                                        "Mitochondrial",
                                                                        "RNA Metabolism",#/Ribosome Biogenesis/Purine Metabolism",
                                                                        "Immune Genes",
                                                                        "Innate/Adaptive Immunity",
                                                                        "Aging Up",
                                                                        "Aging Down",
                                                                        "Aging Cluster Complement",
                                                                        "Aging Composite"
                                                                        
  )) %>% fct_rev
  ex0_m8fdr =
    ex0_m8fdr %>% 
    mutate(gene_set_name = fct_recode(gene_set_name, "Lysosome Metabolism (n=5)"="Lysosome Metabolism",
                                "Innate/Adaptive Immunity (n=68)" = "Innate/Adaptive Immunity",
                                "Fatty Acid Metabolism (n=7)"="Fatty Acid Metabolism",
                                "Actin Regulation (n=9)"="Actin Regulation",
                                "Immune Genes (n=44)"="Immune Genes",
                                "Ribosome (n=12)"="Ribosome",
                                "Mitochondrial (n=29)"="Mitochondrial",
                                "RNA Metabolism (n=40)"="RNA Metabolism",
                                "DNA Replication (n=16)"="DNA Replication",
                                "Aging Up (n=438)"="Aging Up",
                                "Aging Down (n=610)"="Aging Down",
                                "Aging Cluster Complement (n=818)"="Aging Cluster Complement",
                                "Aging Composite (n=1048)"="Aging Composite"))
  # 
  # ex0_m8fdr$gene_set_name <- factor(ex0_m8fdr$gene_set_name,
  #                                   levels = c("1KI", "CTRA", "Antbintf",  "Interferon", "Inflamation in CTRA",                                                                    "Aging", "Aging Down", "Aging Up",
  #                                              "Aging Down Cl1",  "Aging Down Cl1a", "Aging Down Cl1b", "Aging Down Cl1c",
  #                                              "Aging Down Cl2",  "Aging Down Cl3",
  #                                              "Aging Up Cl1",    "Aging Up Cl2",    "Aging Up Cl3",    "Aging Up Cl4"  )) %>% fct_rev
  axiscolor = c("darkblue","darkblue","darkblue","darkblue","grey30", "grey30", "grey30", "grey30", "grey30", "grey30","grey30", "grey30", "grey30")
  panelA =  ex0_m8fdr %>% 
    # mutate(group = ifelse(treatment %>% str_detect("Hispanic"), "Race", "Skin Color")) %>% 
    ggplot(aes(treatment, gene_set_name, size = pval2
               # fill = `1KI Genes`,
               # colour = `1KI Genes`
    )) +
    geom_point( shape = 21,
                # stroke = 1.5,
               # alpha = 0.4,
               # colour = "darkblue", fill = "navy",
               colour = "black",
               fill = "goldenrod3"
               ) +
    # scale_fill_manual(values = c( "goldenrod3", "lightblue","cornflowerblue")) +
    # facet_wrap(~group, scales = "free_x") +
    # scale_fill_manual(values = c("darkblue", "goldenrod3")) +
    # scale_color_manual(values = c("darkblue", "goldenrod3")) +
    # geom_jitter(height = 0.00000025) +
    # gghighlight(class == "inflam") +
    theme_bw() +
    labs(
      # caption = paste("p_eqtl = ", p_eqtl),
      # title = "Figure 1. Associations between Race Ethnicity/Skin Color
      #         and mRNA-Based Disease Signatures, Add Health
      #         (p-values reported, FDR-corrected for whole genome)",
      y = "mRNA Signatures",
      x = "Skin Color") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          # plot.margin=unit(c(1, 1, 0.1, 1), "cm"),
          axis.text.y = element_text(color=axiscolor),
          text = element_text(size=10, face = "bold")) +
    scale_size_continuous(name = "Adjusted P-value", range = c(0, 21), 
                          limits = c(0.0000001, 55000), breaks = c(
                            # 0.0000001,
                            10000, 25000, 55000),
                          labels = c(
                            # "n.s.",
                            "p<0.05", "p<0.01", "p<0.001")) +
    scale_alpha(guide = 'none') +
    guides(shape = guide_legend(override.aes = list(size = 15)),
           fill = guide_legend(override.aes = list(size = 15)),
           size = FALSE
           )  
}

pca_plotting = function(data, p_eqtl, adjmethod) {
  skincolor = outm7pca(p = p_eqtl, control = "ancestryPC_ses", data)
  ex0_m7pca = skincolor %>% 
    filter(gene_set_name %>% str_detect("aging")) %>% 
    filter(gene_set_name!="aging_down_cl1_mRNA") %>% 
    dplyr::select(treatment, gene_set_name, p_pca, p_id) %>% 
    filter(p_pca!="NA")
  if(adjmethod=="fdr") {
    ex0_m7pca = ex0_m7pca %>% 
      group_by(treatment) %>% 
      # using fdr for correction
      mutate(p= p_pca %>% p.adjust("fdr")) %>% 
      ungroup
  } else if(adjmethod =="bonferroni") {
    ex0_m7pca = ex0_m7pca %>% 
      mutate(p = case_when(gene_set_name == "aging_up_cl2_mRNA" ~ p_pca*9*13,
                           gene_set_name == "aging_up_cl3_mRNA" ~ p_pca*7*13,
                           gene_set_name == "aging_up_cl4_mRNA" ~ p_pca*5*13,
                           TRUE ~ p_pca*10*13))
  }
  ex0_m7pca = ex0_m7pca %>%
   # filter(p<0.05) %>% 
    mutate(
      pval2=case_when(
        # p<0.0001 ~ 100000,
        p<0.001 ~ 55000,
        p<0.01 ~ 25000,
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
                                    gene_set_name =="Aging Up Cl2" ~ "Actin Regulation",#/Focal Adhesion/Tight Junctions",
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
                                                            "Lysosome Metabolism",#/Glycosaminoglycan Degradation"
                                                            "Fatty Acid Metabolism",#/Peroxisome Activity", 
                                                            "Actin Regulation",#/Focal Adhesion/Tight Junctions",
                                                            "Ribosome",
                                                            "DNA Replication",#/Elongation/Mismatch Repair",
                                                            "Mitochondrial",
                                                            "RNA Metabolism",#/Ribosome Biogenesis/Purine Metabolism",
                                                            "Immune Genes",
                                                            "Innate/Adaptive Immunity",
                                                            "Aging Up",
                                                            "Aging Down",
                                                            "Aging Cluster Complement",
                                                            "Aging Composite"
           )) %>% fct_rev
           # gene_set_name = str_c(gene_set_name, "_", p_id)
    )
  
  ex0_m7pca =
    ex0_m7pca %>%
    mutate(gene_set_name = fct_recode(gene_set_name, "Lysosome Metabolism (n=5)"="Lysosome Metabolism",
                                      "Innate/Adaptive Immunity (n=68)" = "Innate/Adaptive Immunity",
                                      "Fatty Acid Metabolism (n=7)"="Fatty Acid Metabolism",
                                      "Actin Regulation (n=9)"="Actin Regulation",
                                      "Immune Genes (n=44)"="Immune Genes",
                                      "Ribosome (n=12)"="Ribosome",
                                      "Mitochondrial (n=29)"="Mitochondrial",
                                      "RNA Metabolism (n=40)"="RNA Metabolism",
                                      "DNA Replication (n=16)"="DNA Replication",
                                      "Aging Up (n=438)"="Aging Up",
                                      "Aging Down (n=610)"="Aging Down",
                                      "Aging Cluster Complement (n=818)"="Aging Cluster Complement",
                                      "Aging Composite (n=1048)"="Aging Composite")) %>% 
    group_by(treatment, gene_set_name) %>% 
    mutate(noPC = n()) %>% 
    mutate(noPC = ifelse(noPC>1, noPC, NA)) %>% 
    mutate(group = ifelse(treatment %>% str_detect("Hispanic"), "Race", "Skin Color")) %>% 
    mutate(sample = "Whole Sample") %>% 
    mutate(p_id = ifelse(p<0.05,p_id,NA),
           pval2 = ifelse(p<0.05,pval2,NA))
  # ex0_m8fdr$gene_set_name <- factor(ex0_m8fdr$gene_set_name,
  #                                   levels = c("1KI", "CTRA", "Antbintf",  "Interferon", "Inflamation in CTRA",                                                                    "Aging", "Aging Down", "Aging Up",
  #                                              "Aging Down Cl1",  "Aging Down Cl1a", "Aging Down Cl1b", "Aging Down Cl1c",
  #                                              "Aging Down Cl2",  "Aging Down Cl3",
  #                                              "Aging Up Cl1",    "Aging Up Cl2",    "Aging Up Cl3",    "Aging Up Cl4"  )) %>% fct_rev
  axiscolor = c("darkblue","darkblue","darkblue", "darkblue", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30","grey30", "grey30", "grey30")
  panelB =  ggplot(ex0_m7pca, aes(treatment, gene_set_name, size = pval2,
                                  fill = treatment,
                                  shape = sample
  )) + 
    # geom_point(
    #   # stroke = 1.5, 
    #            shape = 21,
    #            alpha = 0.4, 
    #            # colour = "darkblue",
    #            fill = "navy"
    #            # position = position_jitterdodge()
    #            ) +
    geom_jitter(shape =21, colour = "black", position = position_jitter(width = 0.5, height = 0.15, seed = 1)) +
    # geom_text(aes(label = p_id), 
    #           show.legend=FALSE, 
    #           size=2,
    #           colour = "black",
    #           position = position_jitter(width = 0.5, height = 0.15, seed = 1)) +
    # facet_wrap(~group, scales = "free_x") +
    # geom_text(aes(label=noPC), show.legend=FALSE, size=3, color ="white")+
    # scale_fill_manual(values = c("darkblue", "goldenrod3")) +
    # scale_color_manual(values = c("darkblue", "goldenrod3")) +
    # geom_jitter(height = 0.00000025) +
    # gghighlight(class == "inflam") +
    theme_bw() +
    labs(
      # caption = paste("p_eqtl = ", p_eqtl),
      # title = "Figure . PCA regression",
      fill = "Skin Color",
      y = "mRNA Signatures PC",
      x = "Skin Color") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          # plot.margin=unit(c(1, 1, 0.1, 1), "cm"),
          # set y axis text color
          # axis.text.y = element_text(color=axiscolor),
          # remove y axis title, text and ticks
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          text = element_text(size=10, face = "bold")) +
    scale_fill_manual(values = c("goldenrod1", "goldenrod3", "lightblue","cornflowerblue")) +
    # scale_color_brewer(palette = "Dark2")
    # scale_color_viridis(option = "D", discrete = TRUE)+
    scale_size_continuous(name = "Adjusted P-value", range = c(0, 21), 
                          limits = c(0.0000001, 55000), breaks = c(
                            # 0.0000001,
                            10000, 25000, 55000),
                          labels = c(
                            # "n.s.",
                            "p<0.05", "p<0.01", "p<0.001")) +
    scale_alpha(guide = 'none') +
    guides(shape = guide_legend(override.aes = list(colour = "darkblue", size = 15)),
           fill = guide_legend(override.aes = list(size = 15)),
           colour = guide_legend(override.aes = list(size = 15)),
           size = guide_legend(override.aes = list(colour = "lightgrey")))    
  
  
}


logFC_ploting = function(example, p_eqtl){
  isnt_out_z <- function(x, thres = 3, na.rm = TRUE) {
    abs(x - mean(x, na.rm = na.rm)) <= thres * sd(x, na.rm = na.rm)
  }
  
  
  outm8 = p_eqtl %>%
    f0(data = example) %>% 
    filter(control_set==control) %>% 
    hoist(out, m = list("result", "m8_fdr", 1, "other", "m"))
  
  a = outm8 %>%
    select(p_eqtl, treatment, gene_set_name, m) %>%
    filter(m!="NULL") %>% 
    filter(gene_set_name %>% str_detect("aging")) %>% 
    filter(gene_set_name!="aging_down_cl1_mRNA") %>% 
    mutate_at(.vars = vars("m"),
              .funs = list(~ map(., ~ mutate(., logFC_outlier = isnt_out_z(logFC)) %>% filter(logFC_outlier==TRUE)))) %>% 
    mutate(test = m %>% map(~ wilcox.test(.$logFC, mu = 0, alternative = "two.sided")),
           logFC.t.test.p = test %>% map_dbl(~ .$p.value) %>% p.adjust(method = "fdr") 
           # %>% format(digits = 3, scientific =T)
    ) %>% 
    unnest(m) %>% 
    
    dplyr::select(treatment, gene_set_name, logFC, logFC.t.test.p) %>% 
    mutate(
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
                                    gene_set_name =="Aging Up Cl2" ~ "Actin Regulation",#/Focal Adhesion/Tight Junctions",
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
                                                            "Lysosome Metabolism",#/Glycosaminoglycan Degradation"
                                                            "Fatty Acid Metabolism",#/Peroxisome Activity", 
                                                            "Actin Regulation",#/Focal Adhesion/Tight Junctions",
                                                            "Ribosome",
                                                            "DNA Replication",#/Elongation/Mismatch Repair",
                                                            "Mitochondrial",
                                                            "RNA Metabolism",#/Ribosome Biogenesis/Purine Metabolism",
                                                            "Immune Genes",
                                                            "Innate/Adaptive Immunity",
                                                            "Aging Up",
                                                            "Aging Down",
                                                            "Aging Cluster Complement",
                                                            "Aging Composite"
           )) %>% fct_rev
           # gene_set_name = str_c(gene_set_name, "_", p_id)
    )
  
  axiscolor = c("darkblue","darkblue","darkblue", "darkblue", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30","grey30", "grey30", "grey30")
  a %>% 
    mutate(gene_set_name = fct_recode(gene_set_name, "Lysosome Metabolism (n=5)"="Lysosome Metabolism",
                                      "Innate/Adaptive Immunity (n=68)" = "Innate/Adaptive Immunity",
                                      "Fatty Acid Metabolism (n=7)"="Fatty Acid Metabolism",
                                      "Actin Regulation (n=9)"="Actin Regulation",
                                      "Immune Genes (n=44)"="Immune Genes",
                                      "Ribosome (n=12)"="Ribosome",
                                      "Mitochondrial (n=29)"="Mitochondrial",
                                      "RNA Metabolism (n=40)"="RNA Metabolism",
                                      "DNA Replication (n=16)"="DNA Replication",
                                      "Aging Up (n=438)"="Aging Up",
                                      "Aging Down (n=610)"="Aging Down",
                                      "Aging Cluster Complement (n=818)"="Aging Cluster Complement",
                                      "Aging Composite (n=1048)"="Aging Composite")) %>% 
    mutate(p_sig = ifelse(logFC.t.test.p>0.05, "Non Sig", "Sig")) %>% 
    ggplot(aes( x = logFC, y = gene_set_name,
                color = p_sig,
                fill = p_sig)) +
    geom_violin()+
    facet_wrap( ~  treatment , scales = "free_x", strip.position = "bottom")  +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
    scale_color_manual(values=c("goldenrod3", "lightblue"))+
    scale_fill_manual(values=c("goldenrod3", "lightblue"))+
    theme(
      # legend.position = "none",
      panel.margin.x = unit(1, "lines"),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(color=axiscolor),
      axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
      axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
      plot.title = element_text(size = 20, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
      plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
      strip.text = element_text(size = 7),
      text = element_text(family = "Georgia"))
  
}


source("/home/xu/ses-1/user_wx/geom_split_violin.R")
out_lier_NA <- function(x, thres = 3, na.rm = TRUE) {
  x = ifelse(abs(x - mean(x, na.rm = na.rm)) <= thres * sd(x, na.rm = na.rm), x, NA)
  
}

med_ploting = function(example, p_eqtl){
  isnt_out_z <- function(x, thres = 3, na.rm = TRUE) {
    abs(x - mean(x, na.rm = na.rm)) <= thres * sd(x, na.rm = na.rm)
  }
  
  p_eqtl %>% map_df(outm7med, mediators, control,example)
  
  outm7med = p_eqtl %>%
    f0(data = example) %>% 
    filter(control_set==control) %>% 
    hoist(out, m = list("result", "m7_ob", 1, "mediation")) %>% 
    unnest_longer(m) %>% 
    unnest_longer(m)
  
  a = outm8 %>%
    select(p_eqtl, treatment, gene_set_name, m) %>%
    filter(m!="NULL") %>% 
    filter(gene_set_name %>% str_detect("aging")) %>% 
    filter(gene_set_name!="aging_down_cl1_mRNA") %>% 
    mutate_at(.vars = vars("m"),
              .funs = list(~ map(., ~ mutate(., logFC_outlier = isnt_out_z(logFC)) %>% filter(logFC_outlier==TRUE)))) %>% 
    mutate(test = m %>% map(~ wilcox.test(.$logFC, mu = 0, alternative = "two.sided")),
           logFC.t.test.p = test %>% map_dbl(~ .$p.value) %>% p.adjust(method = "fdr") 
           # %>% format(digits = 3, scientific =T)
    ) %>% 
    unnest(m) %>% 
    
    dplyr::select(treatment, gene_set_name, logFC, logFC.t.test.p) %>% 
    mutate(
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
                                    gene_set_name =="Aging Up Cl2" ~ "Actin Regulation",#/Focal Adhesion/Tight Junctions",
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
                                                            "Lysosome Metabolism",#/Glycosaminoglycan Degradation"
                                                            "Fatty Acid Metabolism",#/Peroxisome Activity", 
                                                            "Actin Regulation",#/Focal Adhesion/Tight Junctions",
                                                            "Ribosome",
                                                            "DNA Replication",#/Elongation/Mismatch Repair",
                                                            "Mitochondrial",
                                                            "RNA Metabolism",#/Ribosome Biogenesis/Purine Metabolism",
                                                            "Immune Genes",
                                                            "Innate/Adaptive Immunity",
                                                            "Aging Up",
                                                            "Aging Down",
                                                            "Aging Cluster Complement",
                                                            "Aging Composite"
           )) %>% fct_rev
           # gene_set_name = str_c(gene_set_name, "_", p_id)
    )
  
  axiscolor = c("darkblue","darkblue","darkblue", "darkblue", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30","grey30", "grey30", "grey30")
  a %>% mutate(gene_set_name = fct_recode(gene_set_name, "Lysosome Metabolism (n=5)"="Lysosome Metabolism",
                                          "Innate/Adaptive Immunity (n=68)" = "Innate/Adaptive Immunity",
                                          "Fatty Acid Metabolism (n=7)"="Fatty Acid Metabolism",
                                          "Actin Regulation (n=9)"="Actin Regulation",
                                          "Immune Genes (n=44)"="Immune Genes",
                                          "Ribosome (n=12)"="Ribosome",
                                          "Mitochondrial (n=29)"="Mitochondrial",
                                          "RNA Metabolism (n=40)"="RNA Metabolism",
                                          "DNA Replication (n=16)"="DNA Replication",
                                          "Aging Up (n=438)"="Aging Up",
                                          "Aging Down (n=610)"="Aging Down",
                                          "Aging Cluster Complement (n=818)"="Aging Cluster Complement",
                                          "Aging Composite (n=1048)"="Aging Composite")) %>% 
    mutate(p_sig = ifelse(logFC.t.test.p>0.05, "Non Sig", "Sig")) %>% 
    ggplot(aes( x = logFC, y = gene_set_name,
                color = p_sig,
                fill = p_sig)) +
    geom_violin()+
    facet_wrap( ~  treatment , scales = "free_x", strip.position = "bottom")  +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
    scale_color_manual(values=c("goldenrod3", "lightblue"))+
    scale_fill_manual(values=c("goldenrod3", "lightblue"))+
    theme(
      # legend.position = "none",
      panel.margin.x = unit(1, "lines"),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(color=axiscolor),
      axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
      axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
      plot.title = element_text(size = 20, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
      plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
      strip.text = element_text(size = 7),
      text = element_text(family = "Georgia"))
  
}

