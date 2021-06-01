library(tidyverse)
library(Biobase)
load_data(reconciled = FALSE, remove_inflam = FALSE)
gene_set_name ="aging_cluster_complement_mRNA"
gene_set = pluck(signatures, "outcome_set", gene_set_name)
genes = t(exprs(dat[gene_set, ]))
data = genes %>% as.tibble()
# out = sparsepca::spca(data, k=10, alpha=1e-3, beta=1e-3, center = TRUE, scale = FALSE, verbose=0)
# doing svd for a data matrix
out = svd(data)
data_re = out$u %*% diag(out$d) %*% t(out$v)
