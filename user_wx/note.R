

# Among non hispanic black and hispanic race and Skin Color interaction model, use hispanic and lightmed as reference, construct the following dummies
# code skincolor white also as 0 as there is ony one 1 white in nonhispanic black group
treatment = c("color_byinterviewer3_DarkBlack", 
              "raceethnicity_NonHblack",
              "raceethnicity__color_byinterviewer3_NonHblack|DarkBlack")

# color3_and_race3_NonHblackhispanic_bespoke_31.03.2021.rds



# Among non hispanic white and hispanic race and Skin Color interaction model, use hispanic and white as reference

treatment = c("color_byinterviewer3_DarkBlack",
              "color_byinterviewer3_LightMed",
              "raceethnicity_NonHwhite",
              # "raceethnicity_Hispanic",
              "raceethnicity__color_byinterviewer3_NonHwhite|DarkBlack",
              "raceethnicity__color_byinterviewer3_NonHwhite|LightMed"
              # "raceethnicity__color_byinterviewer3_Hispanic|DarkBlack",
              # "raceethnicity__color_byinterviewer3_Hispanic|LightMed"
)

# color3_bespoke_29.03.2021 with loose method for tfbm m11 
# color3_bespoke_28.03.2021 <- readRDS("~/ses-1/user_wx/color3_bespoke_28.03.2021.rds")
# color3_bespoke_29.03.2021 <- readRDS("~/ses-1/user_wx/color3_bespoke_29.03.2021.rds")