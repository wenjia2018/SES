define_treatments_and_controls_skincolor_fullsibFE = function(){ 
  
  ############################################################
  # DEFINE "TREATMENTS", CONTROLS, OUTCOMES (INCLUDING TFBMs) 
  ############################################################
  
  # TABLE 1
  # table1 =
  #   c(
  #     "CVD_mRNA",
  #     "inflam1k_mRNA",
  #     "diabetes_mRNA",
  #     "Rheumatoid_Arthritis_mRNA", 
  #     "Alzheimers_mRNA",
  #     "Aortic_Aneurysm_mRNA",
  #     "COPD_mRNA",
  #     "Asthma_mRNA",
  #     "Hypertension_mRNA",
  #     "Depression_mRNA",
  #     "CKD_mRNA"
  #   )
  table1 =
    c(
      "ctra_mRNA",
      "inflame_mRNA",
      "interferon_mRNA",
      "AntBIntF_mRNA",
      # "antibody_mRNA", #only 1 gene
      "inflam1k_mRNA",
      "aging_mRNA",
      "aging_up_mRNA",
      "aging_down_mRNA",
      "aging_cluster_complement_mRNA",
      "aging_down_cl1_mRNA",
      "aging_down_cl1a_mRNA",
      "aging_down_cl1b_mRNA",
      "aging_down_cl1c_mRNA",
      "aging_down_cl2_mRNA",
      "aging_down_cl3_mRNA",
      "aging_up_cl1_mRNA",
      "aging_up_cl2_mRNA",
      "aging_up_cl3_mRNA",
      "aging_up_cl4_mRNA",
      "whole_genome"
    )
  # treatment = c("color_byinterviewer_continuous",
  #               "color_byinterviewer_binary")

  # treatment = c("color_byinterviewer3_DarkBlack",
  #               "color_byinterviewer3_LightMed")
  
    treatment = c(
    # "color_byinterviewer5_White",
    # "color_byinterviewer5_Black",
    "color_byinterviewer5_Dark")
    # "color_byinterviewer5_Medium",
    # "color_byinterviewer5_Light")
  
  # check controls for each analysis!!!!!!!!!
  basic = 
    c(
      "sex_interv", "famid_fullsib",
      "Plate", "AvgCorrelogram100" ,"age_w5",
      # over-representation 
      "NK.cells.activated",
      "T.cells.CD8",
      # under-representation
      "Macrophages.M0", 
      "Macrophages.M2",
      "B.cells.naive",
      "T.cells.CD4.memory.resting",
      "color_byinterviewer5_Black",
      "color_byinterviewer5_Dark",
      "color_byinterviewer5_Medium",
      "color_byinterviewer5_Light"
      
      
       # "color_byinterviewer3_DarkBlack",
       # "color_byinterviewer3_LightMed"
    )
  race = c("raceethnicity_NonHblack",
           "raceethnicity_Hispanic", basic)
  
  ses = c("ses_sss_composite", basic)
  # 
  # ses_race = c("raceethnicity_NonHblack",
  # "raceethnicity_Hispanic", "ses_sss_composite", basic)
  controls = list(basic = basic,
                  # race = race
                  ses = ses
                  # ses_race = ses_race
  )
  
  gene_set_name = signature_names %>% append("whole_genome")

  args = crossing(treatment, gene_set_name, controls)
  
  
  for(i in 1:dim(args)[1]){
    args$controls[i] = args$controls[i] %>% map(setdiff, args$treatment[i])
  }
  
  
  mediators = 
    c(
      "stress_perceived_lm",
      "bills_binary",
      "currentsmoke_binary",
      "w5bmi_lm",
      "insurance_lack_binary",
      "lowbirthweight_binary",
      "high_lowbirth_binary",
      
      "totdiscrim2_binary",  # binary
      "discrim2_binary",   # binary
      "totdiscrim1_category",   # categorical of 4
      # special treatment in mediation
      "totdiscrim1_gamma", #exponential fit glm for mediation
      "totdiscrim2_gamma", #exponential fit glm for mediation
      "countdiscrimwhy_gamma"#exponential fit glm for mediation
      
      # "totdiscrim1_pois", #poisson fit glm for mediation
      # "totdiscrim2_pois", #poisson fit glm for mediation
      # "countdiscrimwhy_pois"#poisson fit glm for mediation
      
      # # "totdiscrim1_category", #ordered logistic fit glm for mediation
      # "totdiscrim2_category", #ordered logistic fit glm for mediation
      # "countdiscrimwhy_category"#ordered logistic fit glm for mediation
    )
  
  immune_tfbms = 
    c(
      "CEBPG_CREB3L1", "CREB3", "CREB3L1", "IRF2", "IRF3",
      "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "JUN", "NFKB1", 
      "NFKB2", "NR3C1"
    )
  
  list2env(mget(ls()), .GlobalEnv)
  
}