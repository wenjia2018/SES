# for categorical treatment with more than 2 levels
if(0){
  ttT = NULL
  
  ttT <- names %>% set_names() %>% map(~
                                         lmFit(y, X) %>%
                                         eBayes %>%
                                         tidy_topTable(of_in = .x))
  
  ttT$omni<-
    lmFit(y, X) %>%
    eBayes %>%
    tidy_topTable(of_in = treatment)
}