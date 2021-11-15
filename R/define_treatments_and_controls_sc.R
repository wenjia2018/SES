define_treatments_and_controls_sc = function(){ 
  
  ############################################################
  # DEFINE "TREATMENTS", CONTROLS, OUTCOMES (INCLUDING TFBMs) 
  ############################################################
  
  # TABLE 1
  table1 =
    c(
      # "CVD_mRNA",
      # "diabetes_mRNA",
      # "inflam1k_mRNA",
      # "Rheumatoid_Arthritis_mRNA", "Alzheimers_mRNA",
      # "Aortic_Aneurysm_mRNA", "COPD_mRNA",
      # "Asthma_mRNA","Hypertension_mRNA",
      # "Depression_mRNA",
      # "CKD_mRNA"
      # "ctra_mRNA",
      # "inflame_mRNA",
      # "interferon_mRNA",
      # "AntBIntF_mRNA",
      # "antibody_mRNA", #only 1 gene
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
      "aging_up_cl4_mRNA"
    )
  
  treatment = c(
    "color_byinterviewer_continuous"
    # "color_byinterviewer5_Black",
    # "color_byinterviewer5_Dark",
    # "color_byinterviewer5_Medium"
    # "color_byinterviewer5_Light"
    # "color_byinterviewer3_DarkBlack",
    # "color_byinterviewer3_LightMed"
  )
  
  controls = 
    list(
      basic = 
        c(
          "sex_interv", "Plate", "age_w5",
          "BirthY", "W5REGION", "pregnant_biow5", 
          "kit_biow5", "tube_biow5",  "FastHrs",
          "travel_biow5",  "months_biow5", "time_biow5"
          # ,
          # "color_byinterviewer5_Black",
          # "color_byinterviewer5_Dark",
          # "color_byinterviewer5_Medium"
          # "color_byinterviewer5_Light"
          # "color_byinterviewer3_DarkBlack",
          # "color_byinterviewer3_LightMed"
        ),
      ancestryPC =
        c(
          "AncestryPC1", "AncestryPC2", "AncestryPC3", "AncestryPC4", "AncestryPC5", "AncestryPC6", "AncestryPC7",
          "AncestryPC8", "AncestryPC9", "AncestryPC10", "AncestryPC11", "AncestryPC12", "AncestryPC13", "AncestryPC14",
          "AncestryPC15", "AncestryPC16", "AncestryPC17", "AncestryPC18", "AncestryPC19", "AncestryPC20"
        )
    ) %>% 
    c(all = list(unique(unlist(.)))) 
  
  gene_set_name = signature_names %>% append("whole_genome_and_tfbm") 
  
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
      # "lowbirthweight_binary",
      # "high_lowbirth_binary",
      
      "totdiscrim2_binary",  # binary
      "discrim2_binary",   # binary
      "totdiscrim1_category" 
    )
  
  immune_tfbms = 
    c(
      "CEBPG_CREB3L1", "CREB3", "CREB3L1", "IRF2", "IRF3",
      "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "JUN", "NFKB1", 
      "NFKB2", "NR3C1"
    )
  
  list2env(mget(ls()), .GlobalEnv)
  
}