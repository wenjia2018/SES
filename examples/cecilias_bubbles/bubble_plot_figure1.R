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

setwd("/home/share/projects/Mike/ses2020.7.17/ses/user_ms")
example0 = readRDS("/home/share/scratch/example0.rds")

ex0<-example0 %>%
  hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
  filter(controls=="all") %>% 
  dplyr::select(treatment, gene_set_name, p)

ex_0 <- ex0 %>%
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
         treatment= case_when(treatment == "edu_p" ~ "Parental Education",
                   treatment =="income_pp1_log" ~  "Parental Income" ,
                   treatment =="SEI_max_p_w12" ~ "Parental SEI",
                   treatment =="ses_composite_pp1" ~ "Parental SES Composite",
                   treatment =="work_collar_rm_f12" ~ "Mother's Occupation",
                   treatment =="work_collar_rf_f12" ~ "Father's Occupation" ,
                   treatment =="work_collar_ff5" ~ "Occupation",
                   treatment =="edu_max" ~ "Education" ,
                   treatment =="income_hh_ff5" ~ "Income"     ,
                   treatment =="SEI_ff5" ~ "SEI"      ,
                   treatment =="ses_sss_composite" ~"SES Composite 4"  ,
                   treatment =="sss_5" ~ "SSS",
                   treatment =="ses_composite_ff5"  ~"SES Composite 3"))  %>%
  mutate(treatment = factor(treatment, levels = c("Parental SES Composite", "Parental Education","Parental Income",
                                                  "Parental SEI", "Mother's Occupation", "Father's Occupation",
                                                  "SES Composite 3", "SES Composite 4", "Education", "Income",
                                                  "SEI", "SSS","Occupation"
  )))


ggplot(ex_0, aes(treatment, gene_set_name, size = pval2, alpha = 0.4)) +
  geom_point(fill = "red", color="navy", stroke = 1.5, alpha=0.4) +
  theme_bw() +
  labs(title = "Figure 1. Associations between Indicators of Socioeconomic Status 
          and mRNA-Based Disease Signatures, Add Health 
          (p-values reported, FDR-corrected for whole genome)",
       y = "mRNA Signatures",
       x = "SES Indicators") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_size_continuous(name = "P-value", range = c(0, 21), 
                       limits = c(0.0000001, 100000), breaks = c(0.0000001, 10000, 15000, 25000, 100000),
                       labels = c("n.s.", "p<0.05", "p<0.01", "p<0.001", "p<0.0001")) 
                      

