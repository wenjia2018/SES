# utiliy functions in race snp PCA process, extraction etc.


get_PC_dim <- function(gene_set_name, p_eqtl) {
  eigenvalue <- data.table::fread(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl, ".omni.pca.eigenval"))
  eigenvec <- data.table::fread(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl, ".omni.pca.eigenvec"))
  kept <- data.table::fread(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl, ".omni.bim"))
  ag <- PCDimension::AuerGervini(eigenvalue$V1, dd = c(dim(eigenvec)[1], dim(kept)[1]))
  d <- PCDimension::agDimension(ag)
  ancestryPC <- str_c("AncestryPC", c(1:d))
  return(ancestryPC)
}


get_PC_dim_general <- function(gene_set_name) {
  eigenvalue <- data.table::fread(str_c("/home/share/dna_ancestry/dna/gene2snp/", gene_set_name,  ".omni.pca.eigenval"))
  eigenvec <- data.table::fread(str_c("/home/share/dna_ancestry/dna/gene2snp/", gene_set_name,  ".omni.pca.eigenvec"))
  kept <- data.table::fread(str_c("/home/share/dna_ancestry/dna/gene2snp/", gene_set_name,  ".omni.bim"))
  ag <- PCDimension::AuerGervini(eigenvalue$V1, dd = c(dim(eigenvec)[1], dim(kept)[1]))
  d <- PCDimension::agDimension(ag)
  ancestryPC <- str_c("AncestryPC", c(1:d))
  return(ancestryPC)
}




# eqtl_pca is used to pre calculate all the pca as it takes at least 20 minitues for calculating one case and
# it does not work within the pipeline, possibly because of the long waiting time
eqtl_pca <- function(gene_set_name, p_eqtl) {
  signatures <- readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_17.11.2020_signature.rds")
  gene_sets <- signatures$outcome_set %>% pluck(gene_set_name)
  G_list <- readRDS("/home/share/data_input/genename_15062020.rds")
  # specific gens of our interest
  gene <- G_list %>%
    filter(hgnc_symbol %in% gene_sets) %>%
    pull(ensembl_gene_id)
  # eqtl catelog
  dat <- readRDS("/home/share/dna_ancestry/eqtl/whole.blood.rds")
  # filter Snps according to some p value threshold p_eqtl
  dat <- dat[dat$pvalue < p_eqtl, ]
  dat <- dat %>% filter(gene_id %in% gene)

  # save Snps which are assoicated with specific genes passing p_eqtl
  data.table::fwrite(unique(dat$rsid) %>% as.data.frame(),
    file = str_c("/home/share/dna_ancestry/dna/keep", gene_set_name, "_", p_eqtl, ".txt"),
    col.names = F
  )

  # keep the snps from data set
  # system("plink --bfile ./dna/omni_joined.freeze3.sharedMarkers --extract ./dna/keep.txt --make-bed --out ./dna/kept.omni")
  system(str_c("plink --bfile /home/share/dna_ancestry/dna/omni_joined.freeze3.sharedMarkers --extract /home/share/dna_ancestry/dna/keep", gene_set_name, "_", p_eqtl, ".txt --make-bed --out /home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl, ".omni"))
  # run PCAs from plink file
  # system("plink --bfile ./dna/kept.omni --pca --out ./dna/kept.omni.pca")
  system(str_c("plink --bfile /home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl, ".omni --pca --out /home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl, ".omni.pca"))
  # output will contain a kept......pca.eigenvec file with family ID, AID, then 20 PCs

  custom_PCA = read.table(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl,".omni.pca.eigenvec"),
                          col.names = c("fid", "AID", "AncestryPC1", "AncestryPC2", "AncestryPC3", "AncestryPC4", "AncestryPC5",
                                        "AncestryPC6", "AncestryPC7", "AncestryPC8", "AncestryPC9", "AncestryPC10", "AncestryPC11",
                                        "AncestryPC12", "AncestryPC13", "AncestryPC14", "AncestryPC15", "AncestryPC16", "AncestryPC17",
                                        "AncestryPC18", "AncestryPC19", "AncestryPC20"))
  custom_PCA %>% saveRDS(str_c("/home/share/dna_ancestry/dna/custom_PCA_", gene_set_name, "_", p_eqtl,".rds"))
  }
