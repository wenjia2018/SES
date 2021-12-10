define_treatments_and_controls = function(){ 
  
  ############################################################
  # DEFINE "TREATMENTS", CONTROLS, OUTCOMES (INCLUDING TFBMs) 
  ############################################################
  
  # TABLE 1
  table1 =
    c(
      "CVD_mRNA",
      "diabetes_mRNA",
      "inflam1k_mRNA",
      # "breast_cancer_mRNA",
      # "Lupus_mRNA", "Colorectal_mRNA",
      "Rheumatoid_Arthritis_mRNA", "Alzheimers_mRNA",
      "Aortic_Aneurysm_mRNA", "COPD_mRNA",
      "Asthma_mRNA","Hypertension_mRNA",
      "Depression_mRNA",
      "CKD_mRNA"
      # "kidney_transplant_tolerance_mRNA"
    )
  
  treatment = c(
    # "ses_sss_composite_3_1",
    # "ses_sss_composite_3_3"
    "ses_sss_composite",
    "sss_5",
    "SEI_ff5",
    "edu_max",
    "income_hh_ff5"
    # "ses_composite_ff5",
    # "ses_composite_pp1",
    # "edu_p",
    # "SEI_max_p_w12",
    # "income_pp1_log"
    # "work_collar_ff5",
    # "work_collar_rf_f12",
    # "work_collar_rm_f12",
  )
  
  controls = 
    list(
      basic = 
        c(
          "sex_interv", "re", "Plate", "age_w5",
          "BirthY", "W5REGION", "pregnant_biow5", 
          "kit_biow5", "tube_biow5",  "FastHrs",
          "travel_biow5",  "months_biow5", "time_biow5"
        ),
      basic_less = 
        c( 
          "sex_interv"
          , "re" 
          ,"age_w5"
          #,"BirthY"
          #, "W5REGION"
          , "pregnant_biow5"
          #, "kit_biow5", "tube_biow5"
          ,  "FastHrs"
          #,"travel_biow5",  "months_biow5", "time_biow5"
          , "Plate",
          "H5INFECT", "H5SUBCLN", "H5CRP8"
        ),
      basic_less_ancestry = 
        c( 
          "sex_interv"
          , "re" 
          ,"age_w5"
          #,"BirthY"
          #, "W5REGION"
          , "pregnant_biow5"
          #, "kit_biow5", "tube_biow5"
          ,  "FastHrs"
          #,"travel_biow5",  "months_biow5", "time_biow5"
          , "Plate"
          , "AncestryPC1", "AncestryPC2", "AncestryPC3", "AncestryPC4",
           "H5INFECT", "H5SUBCLN", "H5CRP8"
        )
    ) %>% 
    c(all = list(unique(unlist(.)))) 
  
  gene_set_name = signature_names %>% append("whole_genome_and_tfbm") 
  
  args = crossing(treatment, gene_set_name, controls)
  
  mediators = 
    c(
      # "drink_category"
      "stress_perceived_lm",
      "w5bmi_lm",
      "bills_binary",
      "currentsmoke_binary",
      "insurance_lack_binary"
    )
  
  immune_tfbms = 
    c(
      "CEBPG_CREB3L1", "CREB3", "CREB3L1", "IRF2", "IRF3",
      "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "JUN", "NFKB1", 
      "NFKB2", "NR3C1"
    )
  
  list2env(mget(ls()), .GlobalEnv)
  
}