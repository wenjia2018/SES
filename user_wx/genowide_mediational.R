fit_m97 = function(mediator){
  y = dat %>% Biobase::exprs() %>% t
  X = dat %>% Biobase::pData() %>%
    dplyr::select(all_of(controls), all_of(treatment), mediator) %>% 
    rename(treatment = treatment,
           mediator = mediator) %>% 
    as_tibble()
  keep = X %>% complete.cases()
  X = X[keep, ]
  y = y[keep, ]
  
  Y = y %>% array_branch(2)
  # gene outcome is dependent measure
  out.y = Y[1:10] %>% map(~ lm(.x ~ ., data = X %>% cbind(.x)))

  # mediate is dependent measure      
  if (is.numeric(mediator)){
    
    out.m = lm(mediator ~., data = X )
    
  }else if(mediator %>% table() %>% length()==2){
    X$mediator =X$mediator %>% as.character() %>%  as.numeric()
    
    out.m = glm(mediator ~., family = "binomial", data = X)
  }else if(datt$mediator %>% table() %>% length()>2){
    X$mediator = X$mediator %>% as.factor()
    # has met some numerical problem with the starting value for optimization
    # therefore not stable, might have again for other mediators outcome combination.
    # https://stackoverflow.com/questions/28916377/r-error-with-polr-initial-value-in-vmmin-is-not-finite
    
    out.m = MASS::polr(mediator ~., data = X, Hess=TRUE)
}

  out<- out.y %>% map(~ mediation::mediate(out.m, .x, treat = "treatment", mediator = "mediator"))
  
  