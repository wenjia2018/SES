library(here)

source(here("R/startup.R"))

## Read raw counts for all batches
raw = readRDS("/home/share/preprocessing/from_Brandt/all.batches.expression.set.070121.Rds")
pheno = pData(raw)
raw = exprs(raw)
genes = rownames(raw)

## Get pheno data for all subjects
df = readRDS("/home/share/preprocessing/preprocessed_two_batches/all.batches.expression.set.rawcount_waves_05.08.2021.rds")
pheno = pData(df)
all.equal(colnames(raw), pheno$AID)

## Convert gene names to HGNC
genedat = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters ="ensembl_gene_id", values = genes, mart = sapiens_ensembl)
rawgenedata = data.frame(Ensembl = genes, stringsAsFactors = F)
rawgenedata$HGNC = NA
for (i in 1:length(rawgenedata$Ensembl)) {
  match = which(genedat$ensembl_gene_id==rawgenedata$Ensembl[i])
  if  (!is_empty(match)) {
    rawgenedata$HGNC[i] = genedat$hgnc_symbol[match[1]]
  }
}
keepgenes = which(rawgenedata$HGNC!="" | rawgenedata$HGNC!=NA)
raw = raw[keepgenes,]
genes = genes[keepgenes]
rawgenedata = rawgenedata[keepgenes,]

## Count exploration and data filtering
keepgenes = filterByExpr(raw, min.count = 10, min.pro = 0.10)
raw = raw[keepgenes,]
genes = genes[keepgenes]
rawgenedata = rawgenedata[keepgenes,]
rownames(raw) = rawgenedata$HGNC

mean_counts <- apply(raw, 1, mean)        
variance_counts <- apply(raw, 1, var)
df <- data.frame(mean_counts, variance_counts)
df_add = df[c(1,2),]
df_add$mean_counts = c(1000000000,1)
df_add$variance_counts = c(1,1000000000)

tiff(str_c(here("temp/"), "Mean_Variance_All_batches_filt.tiff"), units="px", width=(3*650), height=(3*650), res=300)

t1 = ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  geom_blank(data = df_add,aes(x=mean_counts, y=variance_counts)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_abline(intercept = 0, slope = 1, color="red") +
  theme_bw(base_size = 12) + 
  xlab("Mean") + ylab("Variance") + ggtitle("") + 
  theme(plot.title = element_text(size = 14, face = "bold", family = "Calibri", hjust = 0.5)) +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text = element_text(size=12, family = "Calibri")) 
t1
dev.off()

## DE Analysis
# Method #3: LIMMA VOOM
treatment = c("ses_sss_composite","sss_5","SEI_ff5",
              "edu_max","income_hh_ff5")

controls = 
  list(
    basic = 
      c(
        "sex_interv", "re", "Plate", "AvgCorrelogram100" ,"age_w5",
        "BirthY", "W5REGION", "pregnant_biow5", 
        "kit_biow5", "tube_biow5",  "FastHrs",
        "travel_biow5",  "months_biow5", "time_biow5"
      ),
    biological = 
      c( 
        "B.cells.naive", "B.cells.memory", "Plasma.cells",
        "T.cells.CD8", "T.cells.CD4.naive", "T.cells.CD4.memory.resting",
        "T.cells.CD4.memory.activated",
        # "T.cells.follicular.helper",
        "T.cells.regulatory..Tregs.", "T.cells.gamma.delta",
        "NK.cells.resting", "NK.cells.activated", "Monocytes", "Macrophages.M0", 
        # "Macrophages.M1",
        "Macrophages.M2", "Dendritic.cells.resting",
        "Dendritic.cells.activated", "Mast.cells.resting",
        # "Mast.cells.activated", # not estimable in limma
        "Eosinophils", "Neutrophils" 
      ) 
  ) %>% 
  c(all = list(unique(unlist(.)))) 

controls = controls$basic
subData = pheno %>% dplyr::select(ses_sss_composite,  all_of(controls))
non_missing = complete.cases(subData) & subData$sex_interv=="f"
rhs = str_c(c("ses_sss_composite",controls[-c(1)]), collapse = " + ")
model_formula = str_c(" ~ ",rhs) %>% as.formula()
design  = model.matrix(model_formula, data = subData[non_missing, ])

## Filter for the new subject class(only females)
keepgenes = filterByExpr(raw[,non_missing], min.count = 10, min.pro = 0.10)
dge = DGEList(counts = raw[keepgenes,non_missing])
ori_lib_size = dge$samples$lib.size
annoori = var(apply(dge$counts, 2, median))

dge_rle = calcNormFactors(dge, method = "RLE")

#Library sizes
eff_lib_size_rle = dge_rle$samples$norm.factors * colSums(raw[keepgenes,non_missing])

# Voom transformation

v_rle = voom(dge_rle, design, plot = T)

annorle = var(apply(v_rle$E, 2, median))

libdata = data.frame(ID  = colnames(v_rle$E), Value = colSums(v_rle$E)/1e6, stringsAsFactors = F)
libdata = libdata[order(libdata$Value),]
libdata$ID = factor(libdata$ID, levels =  libdata$ID)

tiff(str_c(here("temp/"), "RLE_batches3_LibSize_voom.tiff"), units="px", width=(3*850), height=(3*700), res=300)

t1 = ggplot(libdata, aes(x=ID, y=Value)) +
  geom_line(group = 1, color = "#01665e", size = 1.5) +
  theme_bw(base_size = 12) +
  xlab("Samples") + ylab("Library Size (in millions)") + ggtitle("Sequencing Depth - RLE") + 
  ylim(c(0,0.08)) +
  theme(plot.title = element_text(size = 14, face = "bold", family = "Calibri", hjust = 0.5)) +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=12, family = "Calibri")) + 
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(panel.grid = element_blank())
t1

dev.off()

## DE

fit_rle = lmFit(v_rle, design) 
fit_rle = eBayes(fit_rle)

res_rle = topTable(fit_rle, coef = "ses_sss_composite", n= Inf) %>%
  rownames_to_column(var = "gene")
res_rle = res_rle[res_rle$adj.P.Val<0.05,]

