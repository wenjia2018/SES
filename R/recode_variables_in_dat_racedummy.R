recode_variables_in_dat_racedummy <- function() {
  # IMPURE FUNCTION: NO RETURN
  # only recode variables by special requirements
  # mediator binary(0,1 coded) or continuous
  # treatment binary(do not need to be coded as 0 and 1 but will get a warning message)
  # or continuous
  pData(dat) <- pData(dat) %>%
    # mutate_at(.vars = vars(matches("income|composite|SEI")),
    #           .funs = list(~ .x %>% ntile(4) %>% as.factor())
    # ) %>%
    mutate(raceethnicity = re %>%
             fct_recode(
               NonHwhite = "1",
               # white nonhispanic
               NonHblack = "2",
               # black nonhispanic
               NULL = "3",
               # asian nonhispanic
               NULL = "4",
               # other nonhispanic
               Hispanic = "5"
             ),
           raceethnicity2 = re %>%
             fct_recode(
               NonHwhite = "1",
               # white nonhispanic
               NonHblack = "2",
               # black nonhispanic
               NULL = "3",
               # asian nonhispanic
               NULL = "4",
               # other nonhispanic
               NULL = "5"
             ),
           raceethnicity3 = re %>%
             fct_recode(
               NonHwhite = "1",
               # white nonhispanic
               NULL = "2",
               # black nonhispanic
               NULL = "3",
               # asian nonhispanic
               NULL = "4",
               # other nonhispanic
               Hispanic = "5"
             )) %>%
    mutate_at(
      .vars = vars(matches("^edu_p$|^edu_max$")),
      .funs = list(~ .x %>%
                     factor() %>%
                     fct_collapse(
                       "high or less" = "high",
                       "more than high" = c("votec", "college", "post")
                     ))
    ) %>%
    mutate_at(
      .vars = vars("work_collar_rf_f12", "work_collar_rm_f12", "work_collar_ff5"),
      .funs = list(~ .x %>%
                     factor() %>%
                     fct_collapse(
                       "blue" = c("blue_collar", "farming", "service"),
                       "white" = "white_collar"
                     ))
    ) %>%
    mutate_at(
      .vars = vars("phys_activ_ff5"),
      .funs = list(~ .x %>%
                     factor() %>%
                     fct_collapse(
                       "1" = c("intens", "moderate"),
                       "0" = "none"
                     ))
    ) %>% 
    fastDummies::dummy_cols(select_columns = c("raceethnicity")) %>% 
    mutate(
      totdiscrim1_category = case_when(
        totdiscrim1 %in% c(0, 1, 2, 3) ~ totdiscrim1,
        totdiscrim1 %in% c(4, 5, 6) ~ 4),
      countdiscrimwhy = ifelse(countdiscrimwhy=="no_discrim", NA, countdiscrimwhy) %>% as.numeric(),
      totdiscrim2_binary = ifelse(totdiscrim2 > 0, 1, 0),
      discrim2_binary = ifelse(discrim2 > 0, 1, 0),
      # add 1 for all in oder to fit a gamma distribution in glm (exponential has to be fitted by gamma first then
      # set the dispersion=1)
      totdiscrim1_gamma = totdiscrim1 + 0.001, #exponential fit glm for mediation, add 0.001 to avoid 0
      totdiscrim2_gamma = totdiscrim2 + 0.001, #exponential fit glm for mediation, add 0.001 to avoid 0
      countdiscrimwhy_gamma = countdiscrimwhy + 0.001#exponential fit glm for mediation, add 0.001 to avoid 0
      
    ) %>% 
    mutate_at(vars(c("stress_perceived",
                "w5bmi")), .funs = list(lm = ~ .x)) %>% 
    mutate_at(vars(c("bills",
                "currentsmoke",
                "insurance_lack",
                "lowbirthweight",
                "high_lowbirth")), .funs = list(binary = ~ .x)) %>% 
    mutate_at(vars(c("totdiscrim1",
                     "totdiscrim2",
                     "countdiscrimwhy")), .funs = list(pois = ~ .x)) %>% 
    mutate(H3IR17 = ifelse(re %in% c(1,2,5), H3IR17, NA), #only use white, black, hispanic
      color_byinterviewer_continuous =  case_when(
      H3IR17==5 ~ 0,
      H3IR17==4 ~ 1,
      H3IR17==3 ~ 2,
      H3IR17==2 ~ 3,
      H3IR17==1 ~ 4
    ),
    color_byinterviewer_binary =  case_when(
      H3IR17==5 ~ 0,
      H3IR17==4 ~ 0,
      H3IR17==3 ~ 0,
      H3IR17==2 ~ 1,
      H3IR17==1 ~ 1
    ),
    keep = .$famid_fullsib %in% .$famid_fullsib[duplicated(.$famid_fullsib)],
    # https://stackoverflow.com/questions/16905425/find-duplicate-values-in-r
    famid_fullsib = ifelse(keep == TRUE, famid_fullsib, NA) %>% as.factor,
      color_byinterviewer3 = H3IR17 %>%
             as.character() %>% 
             as.factor %>% 
             fct_collapse(
               DarkBlack = c("1", "2"),
               LightMed = c("3", "4"),
               White = "5") %>%
             relevel(ref = "White"),
    color_byinterviewer5 = H3IR17 %>%
      as.character() %>% 
      as.factor %>% 
      fct_recode(
        White = "5",
        Light = "4",
        Medium = "3",
        Dark = "2",
        Black = "1"
      ) %>% 
      relevel(ref = "White")
             ) %>% fastDummies::dummy_cols(select_columns = c("color_byinterviewer3","color_byinterviewer5")) 
  # %>% 
    # dplyr::select(-starts_with("AncestryPC")) 
  # %>% 
  #   dplyr::left_join(custom_PCA %>% select(-fid))
  # mutate_at(.vars = vars("sss_5"),
  #           .funs = list(~ .x %>% factor))

  dat <<- dat
}