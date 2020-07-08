############################################################
# DEFINE "TREATMENTS", CONTROLS, OUTCOMES (INCLUDING TFBMs) 
############################################################

# TABLE 1
table1 =
  c(
    "CVD_mRNA",
    "diabetes_mRNA",
    "inflam1k_mRNA", "breast_cancer_mRNA",
    "Lupus_mRNA", "Colorectal_mRNA",
    "Rheumatoid_Arthritis_mRNA", "Alzheimers_mRNA",
    "Aortic_Aneurysm_mRNA", "COPD_mRNA",
    "Asthma_mRNA","Hypertension_mRNA",
    "kidney_transplant_tolerance_mRNA"
  )

treatment = c(
  "ses_sss_composite",
  "ses_composite_ff5",
  "sss_5",
  "SEI_ff5",
  "edu_max",
  "income_hh_ff5",
  "work_collar_ff5",
  "work_collar_rf_f12",
  "work_collar_rm_f12",
  "ses_composite_pp1",
  "edu_p",
  "SEI_max_p_w12",
  "income_pp1_log"
)
 
controls = 
  list(
    basic = 
      c(
        "sex_interv", "re", "Plate", "AvgCorrelogram100" ,"age_w5",
        "BirthY", "W5REGION", "pregnant_biow5", 
        "kit_biow5", "tube_biow5",  "FastHrs",
        "travel_biow5",  "months_biow5", "time_biow5"
      ),
    biological = 
      c( 
        "B.cells.naive", "B.cells.memory", "Plasma.cells",
        "T.cells.CD8", "T.cells.CD4.naive", "T.cells.CD4.memory.resting",
        "T.cells.CD4.memory.activated",
        # "T.cells.follicular.helper",
        "T.cells.regulatory..Tregs.", "T.cells.gamma.delta",
        "NK.cells.resting", "NK.cells.activated", "Monocytes", "Macrophages.M0", 
        # "Macrophages.M1",
        "Macrophages.M2", "Dendritic.cells.resting",
        "Dendritic.cells.activated", "Mast.cells.resting",
        # "Mast.cells.activated", # not estimable in limma
        "Eosinophils", "Neutrophils" 
      ) 
  ) %>% 
  c(all = list(unique(unlist(.)))) 

gene_set_name = signature_names %>% append("whole_genome_and_tfbm") 

args = crossing(treatment, gene_set_name, controls)

mediators = 
  c(
  "phys_activ_ff5",
  "bingedrink",
  "currentsmoke",
  "w5bmi"
)

immune_tfbms = 
  c(
    "CEBPG_CREB3L1", "CREB3", "CREB3L1", "IRF2", "IRF3",
    "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "JUN", "NFKB1", 
    "NFKB2", "NR3C1"
  )