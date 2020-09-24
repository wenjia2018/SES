recode_variables_in_dat <- function(){ 
  # IMPURE FUNCTION: NO RETURN
  # only recode variables by special requirements
  # mediator binary(0,1 coded) or continuous
  # treatment binary(do not need to be coded as 0 and 1 but will get a warning message)
  # or continuous
  pData(dat) = pData(dat) %>%
    # mutate_at(.vars = vars(matches("income|composite|SEI")),
    #           .funs = list(~ .x %>% ntile(4) %>% as.factor())
    # ) %>% 
    mutate_at(.vars = vars(matches("^edu_p$|^edu_max$")),
              .funs = list(~ .x %>%
                             factor %>%
                             fct_collapse('high or less' = "high",
                                          'more than high' = c("votec", "college","post")))
    ) %>%
    mutate_at(.vars = vars("work_collar_rf_f12", "work_collar_rm_f12", "work_collar_ff5"),
              .funs = list(~ .x %>% factor %>%
                             fct_collapse('blue' = c("blue_collar","farming","service"),
                                          'white' = "white_collar"))) %>%
    mutate_at(.vars = vars("phys_activ_ff5"),
              .funs = list(~ .x %>% factor %>%
                             fct_collapse("1" = c("intens", "moderate"),
                                          "0" = "none")))
    # mutate_at(.vars = vars("sss_5"),
    #           .funs = list(~ .x %>% factor))
  dat <<- dat
  
}