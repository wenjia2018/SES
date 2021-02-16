#' check if there are any overlap or intersection of the 18 color realted snps and our disease gene set realted snps obtained from eqtl by serie of eqtl_p values
#+ echo=F, eval=T, warning=FALSE, message=FALSE

library(tidyverse)

check_in = function(gene_set_name, p_eqtl){
  kept_eqtl <- data.table::fread(str_c("/home/share/dna_ancestry/dna/kept", gene_set_name, "_", p_eqtl,".omni.bim")) %>% pull("V2")
  kept_spatial <- data.table::fread(str_c("/home/share/dna_ancestry/dna/gene2snp/", gene_set_name, ".omni.bim")) %>% pull("V2")
  index = kept_spatial %in% kept_eqtl
  venn.plot = VennDiagram::venn.diagram(list(eqtl = kept_eqtl, spatial = kept_spatial), filename = NULL, main = str_c(gene_set_name, p_eqtl))
  grid.draw(venn.plot)
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

# color_snp_in_eqtl = 
  args_eqtl %>% 
  mutate(out = pmap(list(gene_set_name = table1, p_eqtl = p_eqtl), safely(check_in)))

# out = color_snp_in_eqtl %>% 
#   hoist(out, snps = "result") %>% 
#   mutate(Num = snps %>% map_int(length))
# 
# dt = out %>% 
#   mutate(p_eqtl = p_eqtl %>% format(scientific = TRUE)) %>% 
#   select(1, 2, 3,5) %>% 
#   filter(Num>0)
# 
# dt %>% 
#   kableExtra::kable() %>%
#   kableExtra::kable_styling()
