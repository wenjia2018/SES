
#' ---
#' title: skin color 5 levels with basic control
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE


#' ## Controls, treatment and threshold
#' 
#' ### controls:
#'   basic = 
#' c(
#'  "sex_interv", "famid_fullsib",
#'   "Plate", "AvgCorrelogram100" ,"age_w5",
#'   # over-representation 
#'  "NK.cells.activated",
#'  "T.cells.CD8",
#'  # under-representation
#'  "Macrophages.M0", 
#'  "Macrophages.M2",
#'  "B.cells.naive",
#'  "T.cells.CD4.memory.resting"
#' )
#' 

#' ### treatments:
#' 
#' * skin color 5 levels: white as reference, Light, Med, Dark, Black as treatments

#' 
#' ## omnibus regression p values are corrected genowide

#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(here)
library(tidyverse)
example_skincolor <- readRDS("/home/xu/ses-1/user_wx/example_skincolor5_aging_FE_23.06.2021.rds")
source("/home/xu/ses-1/user_wx/extract_v2.R")
source("/home/xu/ses-1/user_wx/plotting_utils.R")
example_skincolor %>% 
  hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
  hoist(out, ttT = list("result", "m8_fdr", 1, "other", "m")) %>% 
  filter(ttT != "NULL") %>% 
  mutate(gene_sig = map(ttT, ~ dplyr::filter(., adj.P.Val<0.05) %>% pull(gene))) %>% 
  dplyr::select(treatment, gene_set_name, controls, p, gene_sig) %>% 
  dplyr::filter(p<0.05) %>%
  # left_join(race_ge_tfbm, by = c("treatment", "gene_set_name", "controls")) %>% 
  rename(p_omnibus = p) %>% 
  filter(gene_set_name %>%str_detect("aging")) %>% 
  kableExtra::kable() %>%
  kableExtra::kable_styling()



#' ## omnibus regression logFC plotting
#+ echo=F, eval=T, warning=FALSE, message=FALSE



isnt_out_z <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - mean(x, na.rm = na.rm)) <= thres * sd(x, na.rm = na.rm)
}


outm8 = example_skincolor %>%
  filter(controls =="basic") %>% 
  hoist(out, m = list("result", "m8_fdr", 1, "other", "m"))

a = outm8 %>%
  select(treatment, gene_set_name, m) %>%
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
                         treatment =="color_byinterviewer5_Light"  ~ "Light Brown", 
                         treatment =="color_byinterviewer5_Medium"  ~ "Medium Brown",
                         treatment =="color_byinterviewer5_Dark"  ~ "Dark Brown",
                         treatment =="color_byinterviewer5_Black"  ~ "Black")) %>%
  mutate(treatment = factor(treatment, levels = c("Hispanic", "Non Hispanic Black","Light Brown","Medium Brown","Dark Brown","Black"))) %>%
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
                                                          "Lysosome Metabolism",#/Glycosaminoglycan Degradation"
                                                          "Fatty Acid Metabolism",#/Peroxisome Activity", 
                                                          "Actin Cytoskeleton",#/Focal Adhesion/Tight Junctions",
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
  mutate(p_sig = ifelse(logFC.t.test.p>0.05, "P > .05", "P < .05")) %>% 
  ggplot(aes( x = logFC, y = gene_set_name,
              color = p_sig,
              fill = p_sig)) +
  geom_violin()+
  facet_wrap( ~  treatment , scales = "free_x", strip.position = "bottom")  +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
  scale_color_manual(values=c("goldenrod3", "lightblue"))+
  scale_fill_manual(values=c("goldenrod3", "lightblue"))+
  labs(
    y = "mRNA Signatures",
    x = "logFC") +
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
    strip.text = element_text(size = 7)
    # remove text font because when saving to pdf, cannot recognize this font
    # text = element_text(family = "Georgia")
  ) +
  guides(fill = guide_legend(title = "Wilcoxon Test", override.aes = list(size = 12)),
         color = FALSE)




#' ## PCA regression p values uncorrected, but only the cases where the p value is less than 0.05/10 are presented
#' (no p <0.005)
#+ echo=F, eval=T, warning=FALSE, message=FALSE
threshold = 0.05/10
threshold_med = 0.05


example_skincolor %>%
  hoist(out, p = list("result", "m7_ob", 1, "p")) %>% 
  unnest_longer(p) %>%
  dplyr::select(treatment, gene_set_name, controls, p, p_id) %>% 
  dplyr::filter(p < threshold) %>% 
  kableExtra::kable() %>%
  # kableExtra::column_spec(column = 7, width = "16in", width_min="8in") %>% 
  kableExtra::kable_styling() 



