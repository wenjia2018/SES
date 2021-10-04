library(sparsepca)
signatures <- readRDS("/home/share/preprocessing/preprocessed_two_batches/all.batches.expression.set.tmm_waves_01.09.2021_signature.rds")
dat <- readRDS("/home/share/preprocessing/preprocessed_two_batches/all.batches.expression.set.tmm_waves_01.09.2021.rds")
genes <- t( Biobase::exprs(dat) )[ , signatures$outcome_set$ctra_mRNA ]

sparse.pca = sparsepca::spca(genes, k = 5, scale = TRUE)
sparse.pca2 = sparsepca::spca(scale(genes), k = 5)
summary(sparse.pca)
summary(sparse.pca2)
pca <- prcomp( scale(genes) )

summary(pca)


pca2 = psych::principal(genes, nfactors = 43, rotate = "none", scores = TRUE)
summary(pca2)
centergenes = scale(genes, scale =F)
scores = centergenes %*% sparse.pca$loadings 