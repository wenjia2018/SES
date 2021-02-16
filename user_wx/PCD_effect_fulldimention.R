
fit_old = function(gene_set_name, p_eqtl){
    eigenvalue <- data.table::fread(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl,".omni.pca.eigenval"))
    eigenvec <- data.table::fread(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl,".omni.pca.eigenvec"))
    # header =FALSE other wise the first snp will be col names
    kept <- data.table::fread(str_c("/home/share/dna_ancestry/dna/keep", gene_set_name, "_", p_eqtl, ".txt"), header = FALSE)
    
    ag = PCDimension::AuerGervini(eigenvalue$V1, dd=c(dim(eigenvec)[1], dim(kept)[1]))
    d = PCDimension::agDimension(ag)
   return(d)

}

fit_new = function(gene_set_name, p_eqtl){
  eigenvalue <- data.table::fread(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl,".omni.pca.eigenval"))
  eigenvec <- data.table::fread(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl,".omni.pca.eigenvec"))
  # header =FALSE other wise the first snp will be col names
  kept <- data.table::fread(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl,".omni.bim"))
  ag = PCDimension::AuerGervini(eigenvalue$V1, dd=c(dim(eigenvec)[1], dim(kept)[1]))
  d = PCDimension::agDimension(ag)
  return(d)
  
}

table1 =
  c(
    "ctra_mRNA",
    "inflame_mRNA",
    "interferon_mRNA",
    "AntBIntF_mRNA", 
    # "antibody_mRNA", #only 1 gene
    "inflam1k_mRNA",
    "aging_mRNA",
    "aging_up_mRNA",
    "aging_down_mRNA"
  )

p_eqtl = c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)

args_eqtl = crossing(table1, p_eqtl)
debugonce(fit_old)
old = args_eqtl %>%  mutate(out = pmap(list(gene_set_name = table1, p_eqtl = p_eqtl), safely(fit_old)))
debugonce(fit_new)

new = args_eqtl %>%  mutate(out = pmap(list(gene_set_name = table1, p_eqtl = p_eqtl), safely(fit_new)))
