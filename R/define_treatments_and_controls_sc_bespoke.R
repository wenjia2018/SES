define_treatments_and_controls_sc_bespoke = function(gene_set_name, ancestryPC_keep){ 
  
  ############################################################
  # DEFINE "TREATMENTS", CONTROLS, OUTCOMES (INCLUDING TFBMs) 
  ############################################################
  
 
  treatment = c(
    # "color_byinterviewer_continuous"
    "color_byinterviewer5_Black",
    "color_byinterviewer5_Dark",
    "color_byinterviewer5_Medium",
    "color_byinterviewer5_Light"
    # "color_byinterviewer3_DarkBlack",
    # "color_byinterviewer3_LightMed"
  )
  
      basic = 
        c(
          "sex_interv", "Plate", "age_w5",
          "BirthY", "W5REGION", "pregnant_biow5", 
          "kit_biow5", "tube_biow5",  "FastHrs",
          "travel_biow5",  "months_biow5", "time_biow5"
          ,
          "color_byinterviewer5_Black",
          "color_byinterviewer5_Dark",
          "color_byinterviewer5_Medium",
          "color_byinterviewer5_Light"
          # "color_byinterviewer3_DarkBlack",
          # "color_byinterviewer3_LightMed"
        )
      ancestryPC = c(ancestryPC_keep, basic)
      controls = list(
        basic = basic,
        ancestryPC = ancestryPC)
  
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