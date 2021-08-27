define_treatments_and_controls_bespoke = function(gene_set_name, ancestryPC_keep){ 
  
  ############################################################
  # DEFINE "TREATMENTS", CONTROLS, OUTCOMES (INCLUDING TFBMs) 
  ############################################################
  IV = c(
    "rs3822214",  "rs12203592", "rs2153271",  "rs11198112", "rs4930263",  "rs2376558",
    "rs10896418",
    "rs10765819", "rs9971729",  "rs642742",   "rs12821256", "rs1407995",  "rs12896399", "rs1805005",
    "rs2228479", "rs1110400",  "rs1805008",  "rs885479"
  )
  # treatment = c(
  #   "raceethnicity_NonHblack",
  #   "raceethnicity_Hispanic"
  # )
  treatment = c(
    "color_byinterviewer3_DarkBlack",
    "color_byinterviewer3_LightMed"
    )

  # treatment = c("color_byinterviewer3_DarkBlack", "color_byinterviewer3_LightMed",
  #               "raceethnicity_NonHblack", "raceethnicity_Hispanic")

  # treatment = c("color_byinterviewer_continuous",
  #               "color_byinterviewer_binary")
  
  # 
  # in hispanic and nonhispanic black race group, use hispanic and (lightmed and white) as reference, construct the following dummies
  # treatment = c("color_byinterviewer3_DarkBlack",
  #               # "color_byinterviewer3_LightMed",
  #               "raceethnicity_NonHblack",
  #               # "raceethnicity_Hispanic",
  #               "raceethnicity__color_byinterviewer3_NonHblack|DarkBlack"
  #               # "raceethnicity__color_byinterviewer3_NonHblack|LightMed"
  #               # "raceethnicity__color_byinterviewer3_Hispanic|DarkBlack",
  #               # "raceethnicity__color_byinterviewer3_Hispanic|LightMed"
  #               )

  # in hispanic and nonhispanic white race group, use hispanic and white as reference, construct the following dummies
  # treatment = c("color_byinterviewer3_DarkBlack",
  #               "color_byinterviewer3_LightMed",
  #               "raceethnicity_NonHwhite",
  #               # "raceethnicity_Hispanic",
  #               "raceethnicity__color_byinterviewer3_NonHwhite|DarkBlack",
  #               "raceethnicity__color_byinterviewer3_NonHwhite|LightMed"
  #               # "raceethnicity__color_byinterviewer3_Hispanic|DarkBlack",
  #               # "raceethnicity__color_byinterviewer3_Hispanic|LightMed"
  # )
  # treatment = c(
  #     # "raceethnicity",
  #     "color_byinterviewer3"
  #     # "raceethnicity__color_byinterviewer3"
  #   )
  
  # treatment = c(
  # "color_byinterviewer3_DarkBlack",
  # "color_byinterviewer3_LightMed"
  # # "color_byinterviewer3_White"
  # )
  # treatment = c("color_byinterviewer3_White", "color_byinterviewer3_LightMed")
  # treatment = c(
  #   # "color_byinterviewer5_White",
  #   "color_byinterviewer5_Black",
  #   "color_byinterviewer5_Dark",
  #   "color_byinterviewer5_Medium",
  #   "color_byinterviewer5_Light")
  
  # treatment = c("color_byinterviewer3",
  #               "immigrat",
  #               "color_byinterviewer3__DarkBlack|Immigrate",
  #               "color_byinterviewer3__DarkBlack|Nonimmigrate",
  #               "color_byinterviewer3__LightMed|Immigrate")
  
  basic = 
    c(
      "sex_interv",
      "Plate", "AvgCorrelogram100" ,"age_w5",
      "BirthY", "W5REGION", 
      "pregnant_biow5",
      "kit_biow5",
      "tube_biow5",
      "FastHrs",
      "travel_biow5",  "months_biow5", "time_biow5",
      "B.cells.naive", "B.cells.memory", "Plasma.cells",
      "T.cells.CD8", "T.cells.CD4.naive", "T.cells.CD4.memory.resting",
      "T.cells.CD4.memory.activated",
      # "T.cells.follicular.helper",
      "T.cells.regulatory..Tregs.",
      ############################################## 
      # singleton
      # https://www.stata.com/statalist/archive/2005-10/msg00594.html
      # https://www.statalist.org/forums/forum/general-stata-discussion/general/1478632-why-f-test-is-missing-could-you-please-help
      ############################################## 
      ######################ATTENTION!############################
      ###### in strata regression remove gamma.delta#####
      # "T.cells.gamma.delta",
      "NK.cells.resting", "NK.cells.activated", "Monocytes", "Macrophages.M0", 
      # "Macrophages.M1",
      "Macrophages.M2", "Dendritic.cells.resting",
      "Dendritic.cells.activated", "Mast.cells.resting",
      # "Mast.cells.activated", # not estimable in limma
      "Eosinophils", "Neutrophils",
      
      # nonhispanic black and hispanic group:color race and interaction in dummy variables
      # "color_byinterviewer3_DarkBlack",
      # "raceethnicity_NonHblack",
      # "raceethnicity__color_byinterviewer3_NonHblack|DarkBlack"

      
      # nonhispanic white and hispanic group:color race and interaction in dummy variables 
      # "color_byinterviewer3_DarkBlack",
      # "color_byinterviewer3_LightMed",
      # "raceethnicity_NonHwhite",
      # # "raceethnicity_Hispanic",
      # "raceethnicity__color_byinterviewer3_NonHwhite|DarkBlack",
      # "raceethnicity__color_byinterviewer3_NonHwhite|LightMed"
      # 
      # 
      
      # color race and interaction in single categorical variables
      # "raceethnicity",
      # "color_byinterviewer3"
      # "raceethnicity__color_byinterviewer3"
      
      # only color and race no interaction
      "color_byinterviewer3_DarkBlack",
      "color_byinterviewer3_LightMed"
      # "raceethnicity_NonHblack",
      # "raceethnicity_Hispanic"

      # strata raceethnicity
      # "color_byinterviewer3_DarkBlack",
      # # "color_byinterviewer3_LightMed",
      # "color_byinterviewer3_White"
      # strata raceethnicity color 5 levels
      # "color_byinterviewer5_White",
      # "color_byinterviewer5_Black",
      # "color_byinterviewer5_Dark",
      # "color_byinterviewer5_Medium",
      # "color_byinterviewer5_Light"
      
      # "rs3822214",  "rs12203592", "rs2153271",  "rs11198112", "rs4930263",  "rs2376558",  "rs10896418",
      # "rs10765819", "rs9971729",  "rs642742",   "rs12821256", "rs1407995",  "rs12896399", "rs1805005",
      # "rs2228479", "rs1110400",  "rs1805008",  "rs885479"
      # immigrate skincolor in race strat interaction
      # "color_byinterviewer3_DarkBlack",
      # "immigrat_Nonimmigrate",
      # "color_byinterviewer3__immigrat_DarkBlack|Immigrate",
      # "color_byinterviewer3__immigrat_DarkBlack|Nonimmigrate",
      # "color_byinterviewer3__immigrat_LightMed|Immigrate"
    )
  ses = c("ses_sss_composite", basic)
  ancestryPC = c(ancestryPC_keep, basic)
  ancestryPC_ses = c(ancestryPC_keep, ses)
  
  controls = list(
    basic = basic,
    ses = ses,
    ancestryPC = ancestryPC,
    ancestryPC_ses = ancestryPC_ses)
  
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
      "totdiscrim1_category"   # categorical of 4
      # special treatment in mediation
      # "totdiscrim1_gamma", #exponential fit glm for mediation
      # "totdiscrim2_gamma", #exponential fit glm for mediation
      # "countdiscrimwhy_gamma",#exponential fit glm for mediation
      # # 
      # "totdiscrim1_pois", #poisson fit glm for mediation
      # "totdiscrim2_pois", #poisson fit glm for mediation
      # "countdiscrimwhy_pois"#poisson fit glm for mediation
      # 

    )
  
  immune_tfbms = 
    c(
      "CEBPG_CREB3L1", "CREB3", "CREB3L1", "IRF2", "IRF3",
      "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "JUN", "NFKB1", 
      "NFKB2", "NR3C1"
    )
  
  list2env(mget(ls()), .GlobalEnv)
  
}