# simulation
library(ggformula)
library(tidyverse)

sim = function(){
  
  age = 
    function(t){
      z  = rnorm(n = 1, mean = t, sd = 1) # E(z|t) latent biological age evolves over chronological time t
      L  = c(-2, -.2, .2, 2) 
      mu = z * L # conditional mean E(x|z,t) = E(x|z), note conditional independence from t
      
      x = 
        rnorm(n = length(mu), mean = mu, sd = 1) %>% 
        set_names(str_c("x", 1:length(mu)))
      
      return(out = c(t = t, z = z, x))
    }
  
  t    = rgamma(100, 20, 1/2) 
  IN = map_df(t, age)  
  
  ############################################################
  # IN-SAMPLE
  ############################################################
  
  m = lm(t ~ x1 + x2 + x3 + x4, data = IN)
  
  IN = mutate(IN, 
              pred = predict(m), 
              aa   = pred - t,
              r    = resid(lm(t ~ pred))) 
  
  plot(IN$r, IN$aa) # accelerated aging == residuals
  
  IN %>% GGally::ggpairs()
  
  ############################################################
  # OUT OF SAMPLE
  ############################################################
  
  OUT = map_df(t, age) 
  OUT = mutate(OUT, 
               pred = predict(m, new = OUT),
               aa   = pred - t,              # peter piper
               r    = resid(lm(t ~ pred)))   # simons
  
  OUT %>% GGally::ggpairs()  
  plot(OUT$r, OUT$aa) # accelerated aging ?== residuals
  
  
  return(out =
           list(  
             a_in  = lm(aa ~ z, data = IN),
             b_in  = lm(r  ~ z, data = IN),
             c_in  = lm(aa ~ t, data = IN),
             d_in  = lm(r  ~ t, data = IN),
             a_out = lm(aa ~ z, data = OUT),
             b_out = lm(r  ~ z, data = OUT),
             c_out = lm(aa ~ t, data = OUT),
             d_out = lm(r  ~ t, data = OUT),
             undoubted = lm(r ~ pred, IN) # must be zero
           )
  )
  
}

q = purrr::rerun(1000, sim())

# FACT: correlation between residuals and explanatory variables is zero as a result of —at least in the case of linear regression— how the algorithm works
tmp = map(q, "undoubted") %>% map_df(broom::tidy) %>% filter(term == "pred") 
tmp %>% gf_histogram(~p.value)
tmp %>% gf_histogram(~estimate)

tmp = map(q, "d_in") %>% map_df(broom::tidy) %>% filter(term == "t") 
