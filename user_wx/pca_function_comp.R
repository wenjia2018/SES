x = iris[,1:4]

p1= princomp(x)
p1 %>% summary
p1

p2 = psych::principal(x, nfactors = 4, rotate = "none")
summary(p2)
p2

p3= princomp(scale(x))
p3 %>% summary
p3
