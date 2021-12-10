recode_variables_in_dat <- function() {
  # IMPURE FUNCTION: NO RETURN
  # only recode variables by special requirements
  # mediator binary(0,1 coded) or continuous
  # treatment binary(do not need to be coded as 0 and 1 but will get a warning message)
  # or continuous
  pData(dat) <- pData(dat) %>%
    # mutate_at(.vars = vars(matches("income|composite|SEI")),
    #           .funs = list(~ .x %>% ntile(4) %>% as.factor())
    # ) %>%
    dplyr::mutate(
      H5TO12 = ifelse(H5TO11==0, 0, H5TO12),
      drink_category = case_when(H5TO12 == 0 ~ "never",
                                 H5TO12 %in% c(1,2) ~ "occasion",
                                 H5TO12 %in% c(3,4,5,6) ~ "severe"),
      ses_sss_composite_3 = ses_sss_composite %>% ntile(3),
      raceethnicity = re %>%
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
    mutate_at(vars(c("stress_perceived",
                       "w5bmi")), .funs = list(lm = ~ .x)) %>% 
    # mutate_at(vars(c("")), .funs = list(category = ~ .x)) %>% 
    mutate_at(vars(c("bills",
                     "bingedrink",
                     "currentsmoke",
                     "insurance_lack",
                     "lowbirthweight",
                     "high_lowbirth")), .funs = list(binary = ~ .x)) %>% 
    
    mutate_at(vars(c("H5SUBCLN",
                     "H5CRP8",
                     "H5INFECT")), .funs = list(~ ifelse(.x %in% c(0,1,2,3), .x, NA))) %>% 
  fastDummies::dummy_cols(select_columns = c("ses_sss_composite_3")) %>% 
  mutate_at(.vars = vars("sss_5"),
            .funs = list(~ .x %>% as.numeric()))
  dat <<- dat
}
