data("CigarettesSW", package = "AER")
CigarettesSW$rprice <- with(CigarettesSW, price/cpi)
CigarettesSW$rincome <- with(CigarettesSW, income/population/cpi)
CigarettesSW$tdiff <- with(CigarettesSW, (taxs - tax)/cpi)

fm <- ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff,
            data = CigarettesSW, subset = year == "1995")
summary(fm)
summary(fm, vcov = sandwich, df = Inf, diagnostics = TRUE)

# https://cran.r-project.org/web/packages/ivreg/vignettes/Diagnostics-for-2SLS-Regression.html
# https://www.r-bloggers.com/2013/09/detecting-weak-instruments-in-r/
library(AER)

var = data %>% colnames
var = var[!var=="ses_sss_composite"]
ex = var[1:33]
en = var[34]
ins = var[35:52]
p1 = c(ex,en) %>% str_c(collapse=" + ")
p2 = c(ex, ins) %>% str_c(collapse=" + ")
form = str_c("ses_sss_composite"," ~ ", p1 ,"|", p2)
fm_skincolor = ivreg(form , data = data)



# MR
extract_t = function(x, var) broom::tidy(x) %>% filter(str_detect(term, var))
lm1 = lm(log(rprice) ~ log(rincome) + tdiff,
         data = CigarettesSW, subset = year == "1995") %>% extract_t("tdiff")
lm2 = lm(log(packs) ~  log(rincome)  + tdiff,
         data = CigarettesSW, subset = year == "1995")%>% extract_t("tdiff")
MRInputObject <- mr_input(bx = lm1$estimate,
                          bxse = lm1$std.error,
                          by = lm2$estimate,
                          byse = lm2$std.error)
IVWObject <- mr_ivw(MRInputObject,
                    model = "default",
                    robust = FALSE,
                    penalized = FALSE,
                    correl = FALSE,
                    weights = "simple",
                    psi = 0,
                    distribution = "normal",
                    alpha = 0.05)
IVWObject %>% summary
