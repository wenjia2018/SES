# https://rpubs.com/malshe/214303
# https://stats.idre.ucla.edu/r/dae/multinomial-logistic-regression/
#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(here)
library(tidyverse)
library(rlang)
library(skimr)
library(gtools)
library(limma)
library(recipes)
library(parsnip)
library(workflows)
library(Biobase)
library(nnet)
library(broom)
walk(dir(path = here("R"),full.names = TRUE), source)
# choose normalization methods for downstream analysis
tmm = TRUE
rle = FALSE
log2cpm = FALSE
# which PCA to perform
oblimin = FALSE
nn = TRUE
# explicitly assign ncomp as the smallest number of table signatures gene numbers
ncomp = 10
# for doing genowide DE analysis only
normalization_bydesign = TRUE
# specify if subjects with disease shall be removed
remove_diseased_subjects = TRUE
load_data(reconciled = FALSE, remove_inflam = FALSE)
dt_color_snp <- readRDS("/home/share/dna_ancestry/dna/dt_color_snp.rds")
snp <- pData(dat) %>%
  dplyr::select(AID, 
                # re,
                H3IR17) %>% 
  left_join(dt_color_snp)
#' ## snp from table2 and table 3
#' * in table2 and table3 "rs10896418", "rs10765819","rs9971729","rs1407995","rs12896399" 
#' are left blank for column skin color. blank means unknown.
#' 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
unknown = c("rs10896418", "rs10765819","rs9971729","rs1407995","rs12896399") 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
snp %>% 
  select(-AID, -H3IR17) %>%
  mutate_at(.vars = vars(starts_with("rs")),
            .funs = list(~ ifelse(.=="00", NA, .))) %>% 
  pivot_longer(., cols = colnames(.), names_to = "Var", values_to = "Val") %>% 
  ggplot(aes(x = Val)) +
  geom_bar(stat = "count") +
  facet_wrap(~Var, nrow = 6, ncol =3)

dt <- snp %>%  
  mutate(
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
    # raceethnicity = re %>%
    #               fct_recode(
    #                 white = "1",
    #                 # white nonhispanic
    #                 black = "2",
    #                 # black nonhispanic
    #                 NULL = "3",
    #                 # asian nonhispanic
    #                 NULL = "4",
    #                 # other nonhispanic
    #                 hispanic = "5"
    #               ) %>% relevel(ref = "white")
  ) %>% 
  mutate_at(.vars = vars(starts_with("rs")),
            .funs = list(~ ifelse(.=="00", NA, str_c("_", .)) %>% as.factor))

dt = dt %>% dplyr::select(-AID, -H3IR17)
#' ## multinomial regression:  color_byinterviewer ~ snps, white as reference
#+ echo=F, eval=T, warning=FALSE, message=FALSE, results = FALSE
# Run the Multinomial logistic regression
model_full3 <- nnet::multinom(color_byinterviewer3 ~ rs3822214 + rs12203592 + rs2153271 + rs11198112 + 
                                 rs4930263 + rs2376558 + rs10896418 + rs10765819 + rs9971729 + 
                                 rs642742 + rs12821256 + rs1407995 + rs12896399 + rs1805005 + 
                                 rs2228479 + rs1110400 + rs1805008 + rs885479, data = dt)

# model_full %>% summary()
# coefs <- coef(model_full)
# exp_coefs <- exp(coefs)


# tidy(model_full) %>% 
#   kableExtra::kable() %>%
#   kableExtra::kable_styling()
#+ echo=F, eval=T, warning=FALSE, message=FALSE
#' ### color_byinterviewer as a 3 category factor: White, LightMed, DarkBlack
#+ echo=F, eval=T, warning=FALSE, message=FALSE
temp3 = car::Anova(model_full3)
temp3

#' ### likelihood ratio test to see the overall effect of all the color related snps
#+ echo=F, eval=T, warning=FALSE, message=FALSE
library(lmtest)
keep = dt %>% complete.cases()
dt_null = dt[keep,]
model_null3 = nnet::multinom(color_byinterviewer3 ~ 1, data = dt_null)
lrtest(model_full3,model_null3)


#+ echo=F, eval=T, warning=FALSE, message=FALSE
temp3 = temp3 %>% tidy %>% mutate(sig = stars.pval(p.value))
colnames(temp3) = str_c(colnames(temp3), "_3level")
sig3=temp3 %>% filter(p.value_3level<0.05) %>% pull("term_3level")
# sig3

#' ## skin color 3 level regression coefficient(only for significant snps, and taking the exponential of the original coefficient)
#' for example : $\displaystyle ln\left(\frac{P(skin color = DarkBlack)}{P(skin color = White)}\right) = b_{0} + b_{1}(snp=02) + b_{2}(snp=03)$

#' taking exponential of both sides
#' $\displaystyle \frac{P(skin color = DarkBlack)}{P(skin color = White)} = exp^\left(b_{0} + b_{1}(snp=02) + b_{2}(snp=03)\right)$

#' switching from snprs12203592=01 to snprs12203592=02, the odds for being darkblack vs. white is 2.148585e+09 
#' 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m3 = exp(coef(model_full3)) %>% as.data.frame() 
sig3 %>% map( ~ dplyr::select(m3, contains(.))) %>% bind_cols() %>% t %>% 
  kableExtra::kable() %>% kableExtra::kable_styling()

#' ### color_byinterviewer as a 5 category factor: white, light brown, medium brown, dark brown, black
#+ echo=F, eval=T, warning=FALSE, message=FALSE, results = FALSE


model_full5 <- nnet::multinom(color_byinterviewer5 ~ rs3822214 + rs12203592 + rs2153271 + rs11198112 + 
                                rs4930263 + rs2376558 + rs10896418 + rs10765819 + rs9971729 + 
                                rs642742 + rs12821256 + rs1407995 + rs12896399 + rs1805005 + 
                                rs2228479 + rs1110400 + rs1805008 + rs885479, data = dt)

#+ echo=F, eval=T, warning=FALSE, message=FALSE

temp5 = car::Anova(model_full5)
temp5

temp5 = temp5 %>% tidy %>% mutate(sig = stars.pval(p.value))
colnames(temp5) = str_c(colnames(temp5), "_5level")

sig5 = temp5 %>% filter(p.value_5level<0.05) %>% pull("term_5level")
# sig5

#' ## skin color 5 level regression coefficient(only for significant snps, and taking the exponential of the original coefficient)
#' for example : $\displaystyle ln\left(\frac{P(skin color = Black)}{P(skin color = White)}\right) = b_{0} + b_{1}(snp=02) + b_{2}(snp=03)$

#' taking exponential of both sides
#' $\displaystyle \frac{P(skin color = Black)}{P(skin color = White)} = exp^\left(b_{0} + b_{1}(snp=02) + b_{2}(snp=03)\right)$
#+ echo=F, eval=T, warning=FALSE, message=FALSE
m5 = exp(coef(model_full5)) %>% as.data.frame() 
sig5 %>% map( ~ dplyr::select(m5, contains(.))) %>% bind_cols() %>% t %>% 
  kableExtra::kable() %>% kableExtra::kable_styling()


# ## comparison
# + echo=F, eval=T, warning=FALSE, message=FALS
# temp3 = temp3 %>% tidy %>% mutate(sig = stars.pval(p.value))
# colnames(temp3) = str_c(colnames(temp3), "_3level")
# temp5 = temp5 %>% tidy %>% mutate(sig = stars.pval(p.value))
# colnames(temp5) = str_c(colnames(temp5), "_5level")
# temp = temp3 %>%
#   bind_cols(temp5) 
# temp %>%
#   kableExtra::kable() %>%
#   kableExtra::kable_styling() 




# Unknown snps
#+ echo=F, eval=T, warning=FALSE, message=FALSE
# temp %>% tidy %>% filter(term %in% unknown)


## model comparision

#+ echo=F, eval=T, warning=FALSE, message=FALSE
# keep = dt %>% complete.cases()
# dt_complete = dt[keep, ]
# model_reduced <- nnet::multinom(color_byinterviewer3 ~ 1, data = dt_complete)
# 
# 
# anova(model_reduced, model_full)

