library(tidyverse)
library(snpGeneSets)
signatures <- readRDS("/home/share/preprocessed_two_batches/dt_batches1_2_steve_waves_17.11.2020_signature.rds")
snp_inhouse <- data.table::fread("/home/share/dna_ancestry/dna/omni_joined.freeze3.sharedMarkers.bim") %>% pull("V2")
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

gene_sets=signatures$outcome_set[table1]


snps = gene_sets %>%
  map(gene2snp, isGeneID = FALSE)

out = snps %>% map(~ .$snp[.$snp %in% snp_inhouse])
out %>% saveRDS("/home/share/dna_ancestry/dna/gene2snp/table1snp.rds")

# check how many snps are in house
inhouse = out %>% map(length)
fromPackge = snps %>% map(~ .$snp %>% length)
inhouse %>% cbind(fromPackge)

out = readRDS("/home/share/dna_ancestry/dna/gene2snp/table1snp.rds")
gene2snpPCA = function(gene_set_name){
  # keep the snps for each signature in .txt
  data.table::fwrite(out %>% pluck(gene_set_name) %>% as.data.frame(),
                     file = str_c("/home/share/dna_ancestry/dna/gene2snp/", gene_set_name, ".txt"),
                     col.names = F)
  
  # make plink file containing only the snps obtained through gene2snp
  system(str_c("plink --bfile /home/share/dna_ancestry/dna/omni_joined.freeze3.sharedMarkers --extract /home/share/dna_ancestry/dna/gene2snp/", gene_set_name, ".txt --make-bed --out /home/share/dna_ancestry/dna/gene2snp/", gene_set_name, ".omni"))
  # run PCAs from plink file
  system(str_c("plink --bfile /home/share/dna_ancestry/dna/gene2snp/", gene_set_name, ".omni --pca --out /home/share/dna_ancestry/dna/gene2snp/", gene_set_name, ".omni.pca"))
  # output will contain a kept......pca.eigenvec file with family ID, AID, then 20 PCs
  
  
  custom_PCA = read.table(str_c("/home/share/dna_ancestry/dna/gene2snp/", gene_set_name,".omni.pca.eigenvec"),
                          col.names = c("fid", "AID", "AncestryPC1", "AncestryPC2", "AncestryPC3", "AncestryPC4", "AncestryPC5",
                                        "AncestryPC6", "AncestryPC7", "AncestryPC8", "AncestryPC9", "AncestryPC10", "AncestryPC11",
                                        "AncestryPC12", "AncestryPC13", "AncestryPC14", "AncestryPC15", "AncestryPC16", "AncestryPC17",
                                        "AncestryPC18", "AncestryPC19", "AncestryPC20"))
  
  custom_PCA %>% saveRDS(str_c("/home/share/dna_ancestry/dna/gene2snp/", gene_set_name,".rds"))
}
table1 %>% map(gene2snpPCA)



