define_treatments_and_controls_snps = function(gene_set_name, ancestryPC_keep){ 
  
  ############################################################
  # DEFINE "TREATMENTS", CONTROLS, OUTCOMES (INCLUDING TFBMs) 
  ############################################################
  
  treatment = c(
    "rs3822214",  "rs12203592", "rs2153271",  "rs11198112", "rs4930263",  "rs2376558",  "rs10896418",
    "rs10765819", "rs9971729",  "rs642742",   "rs12821256", "rs1407995",  "rs12896399", "rs1805005",
    "rs2228479", "rs1110400",  "rs1805008",  "rs885479" 
  )
  
  
  basic = 
    c(
      "sex_interv",
      "Plate", "AvgCorrelogram100" ,"age_w5",
      "BirthY", "W5REGION", "pregnant_biow5", 
      "kit_biow5", "tube_biow5",  "FastHrs",
      "travel_biow5",  "months_biow5", "time_biow5",
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
  ses = c("ses_sss_composite", basic)
  ancestryPC = c(ancestryPC_keep, ses)
  
  controls = list(basic = basic, ses = ses, ancestryPC = ancestryPC)
  
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
      "countdiscrimwhy_gamma",#exponential fit glm for mediation
      
      "totdiscrim1_pois", #poisson fit glm for mediation
      "totdiscrim2_pois", #poisson fit glm for mediation
      "countdiscrimwhy_pois",#poisson fit glm for mediation
      
      "totdiscrim1_category", #ordered logistic fit glm for mediation
      "totdiscrim2_category", #ordered logistic fit glm for mediation
      "countdiscrimwhy_category"#ordered logistic fit glm for mediation
    )
  
  immune_tfbms = 
    c(
      "CEBPG_CREB3L1", "CREB3", "CREB3L1", "IRF2", "IRF3",
      "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "JUN", "NFKB1", 
      "NFKB2", "NR3C1"
    )
  
  list2env(mget(ls()), .GlobalEnv)
  
}