library(tidyverse)
k = 10
n = 32
fit1 = lm(disp ~. , mtcars)
fit1  = fit1 %>% summary
rs1 = fit1$r.squared 
fit2 = lm(disp ~. , mtcars %>% select(-drat))
fit2 = fit2 %>% summary
rs2 = fit2$r.squared
f = (rs1 - rs2)/((1-rs1) / (n - k))
qf(0.95, 1, n - k)
pf(f, 1, n - k, lower.tail = FALSE)
