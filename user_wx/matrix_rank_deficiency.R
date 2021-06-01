# design matrix rank deficiency

limma::is.fullrank(X)
limma::nonEstimable(X)
ordinal::drop.coef(X)


rankifremoved <- sapply(1:ncol(X), function (x) qr(X[,-x])$rank)
which(rankifremoved == max(rankifremoved))