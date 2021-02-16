
# m7
fit <- lm(mpg ~ . ,data=mtcars, x = TRUE)
partialFtest(fit[["x"]], mtcars$mpg, c(1:2))
# m8
partialFtest(X %>% as.matrix(), y[1,], c(1:4))