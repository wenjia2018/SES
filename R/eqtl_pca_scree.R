  eqtl_pca_scree = function(gene_set_name, p_eqtl){
    signatures = readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_17.11.2020_signature.rds")
    gene_sets  = signatures$outcome_set %>% pluck(gene_set_name)
    G_list     = readRDS("/home/share/data_input/genename_15062020.rds")
    # specific gens of our interest
    gene       = G_list %>% filter(hgnc_symbol %in% gene_sets) %>% pull(ensembl_gene_id)
    # eqtl catelog
    dat = readRDS("/home/share/dna_ancestry/eqtl/whole.blood.rds")
    # filter Snps according to some p value threshold p_eqtl
    dat = dat[dat$pvalue < p_eqtl, ]
    dat = dat %>% filter(gene_id %in% gene) 
    
    # save Snps which are assoicated with specific genes passing p_eqtl
    data.table::fwrite(unique(dat$rsid) %>% as.data.frame(),
                       file = str_c("/home/share/dna_ancestry/dna/keep", gene_set_name, "_", p_eqtl, ".txt"),
                       col.names = F)
    
    
    # keep the snps from data set
    # system("plink --bfile ./dna/omni_joined.freeze3.sharedMarkers --extract ./dna/keep.txt --make-bed --out ./dna/kept.omni")
    system(str_c("plink --bfile /home/share/dna_ancestry/dna/omni_joined.freeze3.sharedMarkers --extract /home/share/dna_ancestry/dna/keep", gene_set_name, "_", p_eqtl,".txt --make-bed --out /home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl,".omni"))
    # run PCAs from plink file
    # system("plink --bfile ./dna/kept.omni --pca --out ./dna/kept.omni.pca")
    system(str_c("plink --bfile /home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl, ".omni --pca --out /home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl,".omni.pca"))
    # output will contain a kept......pca.eigenvec file with family ID, AID, then 20 PCs
    eigenvalue <- data.table::fread(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl,".omni.pca.eigenval"))
    eigenvec <- data.table::fread(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl,".omni.pca.eigenvec"))
    kept <- data.table::fread(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl,".omni.bim"))
    
    ag = PCDimension::AuerGervini(eigenvalue$V1, dd=c(dim(eigenvec)[1], dim(kept)[1]))
    d = PCDimension::agDimension(ag)
    ancestryPC_full = dput(str_c("AncestryPC", c(1:20)))
    custom_PCA = read.table(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl,".omni.pca.eigenvec"),
                            col.names = c("fid", "AID", ancestryPC_full))
    
    custom_PCA %>% saveRDS(str_c("/home/share/dna_ancestry/dna/custom_PCA_", gene_set_name, "_", p_eqtl,".rds"))
    
    ancestryPC = ancestryPC_full[1:d]
    
  } 
  
  
  
  
  # ancestryPC_keep = eqtl_pca_scree(gene_set_name, p_eqtl)