
remove_diseased_subjects_from_datt <- function(datt, gene_set_name, controls){ 
  # IMPURE FUNCTION: NO RETURN
  # ONLY CALLED FOR ITS SIDE EFFECTS ON datt
  # only keep subjects without corresponding diseases
  keep = 
    switch(
      gene_set_name,
      diabetes_mRNA = dat$diabetes != 1,
      CVD_mRNA = dat$heartatk != 1,
      Asthma_mRNA = dat$H5ID6F != 1 & dat$H5ID6FM != 1,
      Hypertension_mRNA = dat$H5ID6C != 1 & dat$H5ID6CM != 1,
      Aortic_Aneurysm_mRNA = dat$H5ID6Q != 1 & dat$H5ID6QM != 1,
      Melanoma_mRNA = dat$H5ID6A != 1 & dat$H5ID6AM != 1,
      Prostate_mRNA = dat$sex_interv != "f",
      breast_cancer_mRNA = dat$sex_interv == "f",
      breast_cancer_up_mRNA = dat$sex_interv == "f",
      breast_cancer_down_mRNA = dat$sex_interv == "f",
      rep(TRUE, dim(datt)[1])
    )
  datt <- datt[keep, ]
  
  # remove covariates
  if(str_detect(gene_set_name, "breast|Prostate") && is.element("sex_interv", controls)) datt <- select(datt, -sex_interv)
  if(str_detect(gene_set_name, "Prostate") && is.element("pregnant_biow5", controls)) datt <- select(datt, -pregnant_biow5)
  
  return(datt = datt)
}