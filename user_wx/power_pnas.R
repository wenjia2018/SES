
library(tidyverse)
library(ssizeRNA)
library(edgeR)
library(Biobase)
data(hammer.eset)
## load hammer dataset (Hammer, P. et al., 2010)

counts <- exprs(hammer.eset)[, phenoData(hammer.eset)$Time == "2 weeks"]
counts <- counts[rowSums(counts) > 0,]
trt <- hammer.eset$protocol[which(hammer.eset$Time == "2 weeks")] 

mu <- apply(counts[, trt == "control"], 1, mean)  
## average read count in control group for each gene

d <- DGEList(counts)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
disp <- d$tagwise.dispersion      
## dispersion for each gene


## fixed fold change
fc <- 1.1
size <- ssizeRNA_vary(mu = mu, disp = disp, fc = fc, pi0 = 0.8, m = 30,
                      maxN = 2500, replace = FALSE)

size$ssize

# prepare our sample size to others
ns = c(ourstudy = 4000, secondLargest = 1033, thirdLargest = 268)


# from other paper example

tpm <- read_tsv('https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2019_gene_TPM_03.tsv.gz')
metadata <- read_tsv('https://hpc.nih.gov/~mcgaugheyd/eyeIntegration/2019_metadata_03.tsv.gz')

group <- metadata %>% 
  filter(sample_accession %in% names(tpm)) %>% 
  select(sample_accession, Tissue) %>% 
  unique() %>% 
  pull(Tissue)
d <- DGEList(tpm[,2:ncol(tpm)], group = group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)

power <- c()
for (i in seq(5,100,5)){
  power <- c(power, check.power(nGenes = 37000, pi0 = 0.6, i, d$counts[,1], d$tagwise.dispersion, 2, up = 0.5,  replace = TRUE, fdr = 0.05, sims = 5)$pow_qvalue_ave)}

power_data <- as_tibble(cbind(seq(5,100,5), power))

