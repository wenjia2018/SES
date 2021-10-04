library(tidyverse)
library(data.table)
library(mma)
datt <- readRDS("/home/share/scratch/mma_example.rds")
gene_set = "d1"

controls = c(
        "sex_interv", "re", "Plate", "age_w5",
        "BirthY", "W5REGION", "pregnant_biow5", 
        "kit_biow5", "tube_biow5",  "FastHrs",
        "travel_biow5",  "months_biow5", "time_biow5"
      )
mediators = 
  c(
    "stress_perceived_lm",
    "w5bmi_lm",
    "bills_binary",
    "currentsmoke_binary",
    "insurance_lack_binary"
  )


x = datt %>% dplyr::select(all_of(mediators), all_of(controls))
contmed =  which(colnames(x)%like%"_lm")
binmed = which(colnames(x)%like%"_binary")
med_index = 1:length(mediators)
y = datt %>% dplyr::select(!!gene_set)
pred = datt %>%dplyr::select(treatment)
# m1 <- mma::mma(x = x, y = y, pred = pred, contmed = contmed, binmed = binmed, jointm = list(n = 1, j1 = med_index),
               # n2 = 10, alpha = 0.05 )




xx = datt %>% dplyr::select(treatment, all_of(mediators), all_of(controls))
vars = colnames(xx)



library(medflex)
imp_formula = as.formula(str_c(gene_set, "~", vars %>% str_c(collapse = " + ")))
impdata = neImpute(imp_formula,
                   family = gaussian,
                   nMed = 5,
                   data = datt %>% na.omit()
)

model_formula=as.formula(str_c(gene_set, "~", "treatment0", "+", "treatment1", "+",
                               controls %>% str_c(collapse = " + ")))
expData = impdata
neMod=neModel(model_formula,
              family = gaussian,
              expData = impdata,
              nBoot = 10
)
datt = datt %>% na.omit()
library(multimediate)
M1reg = lm(as.formula(str_c("w5bmi_lm ~ treatment +", controls %>% str_c(collapse = " + "))), datt)
M2reg = lm(as.formula(str_c("stress_perceived_lm ~ treatment +", controls %>% str_c(collapse = " + "))), datt)
M3reg = glm(as.formula(str_c("bills_binary ~ treatment +", controls %>% str_c(collapse = " + "))), datt, family = binomial("logit"))
M4reg = glm(as.formula(str_c("currentsmoke_binary ~ treatment +", controls %>% str_c(collapse = " + "))), datt, family = binomial("logit"))
M5reg = glm(as.formula(str_c("insurance_lack_binary ~ treatment +", controls %>% str_c(collapse = " + "))), datt, family = binomial("logit"))

Yreg = lm(as.formula(str_c(gene_set, "~", vars %>% str_c(collapse = " + "))), datt)
med.analysis = multimediate::multimediate(lmodel.m = list(M1reg, M2reg, M3reg, M4reg, M5reg),
                          correlated = FALSE,
                          model.y = Yreg,
                          treat = "treatment",
                          treat.value = 1,
                          control.value = 0,
                          J=10,
                          conf.level=0.95)
