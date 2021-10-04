library("BiocManager")
# BiocManager::install("mdqc")
library("mdqc")
library("sparsepca")
select.lambda <- function( lambda, weight = 0.99 ){
  sapply( lambda, function( lam ){
    spa <- spca( rescale2( genes ), alpha = lam )
    weight * summary(spa)[4,] + (1-weight) * colMeans( spa$loadings == 0 )
  })
}
rescale2 <- function(x) ( x - mean(x) ) / sd(x)
signatures <- readRDS("/home/share/preprocessing/preprocessed_two_batches/all.batches.expression.set.tmm_waves_01.09.2021_signature.rds")
dat <- readRDS("/home/share/preprocessing/preprocessed_two_batches/all.batches.expression.set.tmm_waves_01.09.2021.rds")
genes <- t( Biobase::exprs(dat) )[ , signatures$outcome_set$ctra_mRNA ]
n <- 200
p <- 20
genes <- rescale2( genes[1:n, 1:p] )
#============================
# Exploratory data analysis:
#============================
pairs(genes[, 1:5])
plot( genes[, 1] ~ genes[, 2] )
plot( genes[, 3] ~ genes[, 4] )
hist(genes[, ceiling(runif(1) * p )])
#============================
# Dimension reduction:
#============================
# PCA:
pca <- prcomp( genes )
# Sparse PCA:
lam <- seq( 0.1, 20, length.out = 30 ) / 10000
crit <- select.lambda( lam, weight = 0.99 )
plot( crit[1, ] ~ lam, type = "l" )
lam.opt <- lam[ which.max(crit[1, ] ) ]
spa <- spca( rescale2( genes ), alpha = lam.opt )
colMeans( spa$loadings != 0 )
spa$loadings
# comparisons:
lapply( p:1, function(j){
  plot( pca$x[, j], main = j )
  lines( spa$scores[, j], col = "red" )
})
lapply( p:1, function(j){
  plot( spa$scores[, j] ~ pca$x[, j] , main = j )
})
# Explained variance:
plot(summary(spa)[4, ], type = "l")
plot(summary(pca)$importance[3, ], type = "l")

