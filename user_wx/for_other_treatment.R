library(here)
# dev.off()
source(here("cluster/startup.R"))
set.seed(1234)
# rravi_path = "/home/rravi/Projects/PLOS/"
counts <- readRDS(str_c(here(),"/Res/WGCNA.rds"))
signatures <- readRDS("/home/share/preprocessing/preprocessed_two_batches/all.batches.expression.set.tmm_waves_01.09.2021_signature.rds")

signatures = signatures[[1]]
# sel_signatures = c(6,7,8,21,32,39,46,35,38,64,65)
sel_signatures = c(47)
sel_signatures = c(6,7,8,21,32,39,46,35,38,64,65)
## Get pheno data for all subjects
# df = readRDS("/home/share/preprocessing/preprocessed_two_batches/all.batches.expression.set.rawcount_waves_05.08.2021.rds")
df = readRDS("/home/share/preprocessing/preprocessed_two_batches/all.batches.expression.set.tmm_waves_01.09.2021.rds")
pheno = pData(df)
pheno = pheno[pheno$AID %in% colnames(counts),]
all.equal(colnames(counts), pheno$AID)
pheno = 
  pheno %>% mutate(
    raceethnicity = re %>%
      fct_recode(
        NonHwhite = "1",
        # white nonhispanic
        NonHblack = "2",
        # black nonhispanic
        NULL = "3",
        # asian nonhispanic
        NULL = "4",
        # other nonhispanic
        Hispanic = "5"
      ),
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
      ) )%>% fastDummies::dummy_cols(select_columns = c("color_byinterviewer3","color_byinterviewer5")) 
## DE Analysis
# Method #3: LIMMA VOOM
treatment = c(          
  # "color_byinterviewer3_DarkBlack",
  # "color_byinterviewer3_LightMed"
  "color_byinterviewer5_Black",
  "color_byinterviewer5_Dark",
  "color_byinterviewer5_Medium"
)

controls = 
  list(
    basic = 
      c(
        # "color_byinterviewer3_DarkBlack",
        # "color_byinterviewer3_LightMed",
        "color_byinterviewer5_Black",
        "color_byinterviewer5_Dark",
        "color_byinterviewer5_Medium",
        "sex_interv","age_w5",
        "BirthY", "W5REGION", "pregnant_biow5", 
        "kit_biow5", "tube_biow5",  "FastHrs",
        "travel_biow5",  "months_biow5", "time_biow5", "Plate"
      ),
    ancestryPC =
      c(
        "AncestryPC1", "AncestryPC2", "AncestryPC3", "AncestryPC4"
        # ,
        # "AncestryPC5", "AncestryPC6", "AncestryPC7", 
        # "AncestryPC8", "AncestryPC9", "AncestryPC10", "AncestryPC11", "AncestryPC12", "AncestryPC13", "AncestryPC14", 
        # "AncestryPC15", "AncestryPC16", "AncestryPC17", "AncestryPC18", "AncestryPC19", "AncestryPC20"
      )
  ) %>% 
  c(all = list(unique(unlist(.)))) 

controls = controls$basic
subData = pheno %>%
  filter(raceethnicity =="NonHblack" &
           color_byinterviewer5 !="White" ) %>% 
  dplyr::select(AID, all_of(treatment), all_of(controls))
subData = droplevels(subData)

full_clus_res <- readRDS("~/ses-1/Res/Clustering_aging/1k/full_clus_res.rds")
average_expr_all <- readRDS("~/ses-1/Res/Clustering_aging/1k/average_expr.rds")

full_clus_res <- readRDS("~/ses-1/Res/Clustering_aging/without_1k/full_clus_res.rds")
average_expr_all <- readRDS("~/ses-1/Res/Clustering_aging/without_1k/average_expr.rds")

full_clus_res <- readRDS("~/ses-1/Res/Clustering/1k/full_clus_res.rds")
average_expr_all <- readRDS("~/ses-1/Res/Clustering/1k/average_expr.rds")



res_sig = c()
res_treatment = c()
res_clus = list()
res_tab =list()
names = c()


for (i in 1:length(sel_signatures)) {

    
     hclust_analysis = full_clus_res[[i]]
     avg_matrix = average_expr_all[[i]] 
     avg_matrix = avg_matrix[, subData$AID]    
    for (j in 1:length(treatment)) {
      rhs = str_c(c(treatment[[j]],controls %>% str_subset(treatment[[j]], negate = TRUE)), collapse = " + ")
      model_formula = str_c(" ~ ",rhs) %>% as.formula()
      design  = model.matrix(model_formula, data = subData)
      
      fit = lmFit(avg_matrix, design)
      fit = eBayes(fit, trend = T)
      res = topTable(fit, coef = treatment[[j]], n= Inf)
      res = res[order(rownames(res)),]
      res_sig = c(res_sig,names(signatures)[[sel_signatures[[i]]]])
      res_treatment = c(res_treatment, treatment[[j]])
      res_clus[[length(res_clus)+1]] = hclust_analysis
      res_tab[[length(res_tab)+1]] =  res
      
    }
    gc()
    rm("ncounts")
} 

full_clus_res = setNames(full_clus_res, names)
average_expr_all = setNames(average_expr_all, names)


# save new dataset for treatment and subset specific case which is different from the case used to generate the clusters

saveRDS(res_clus, str_c(here("Res/Clustering_aging/1k/"), "res_clus.rds"))
saveRDS(res_sig, str_c(here("Res/Clustering_aging/1k/"), "res_sig.rds"))
saveRDS(res_treatment, str_c(here("Res/Clustering_aging/1k/"), "res_treatment.rds"))
saveRDS(res_tab, str_c(here("Res/Clustering_aging/1k/"), "res_tab.rds"))
saveRDS(full_clus_res, str_c(here("Res/Clustering_aging/1k/"), "full_clus_res.rds"))



# use for mediation.R to get significant 

saveRDS(average_expr_all, str_c(here("Res/Clustering_aging/1k/"), "average_expr.rds"))