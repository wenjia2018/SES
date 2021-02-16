
# this is the first attemp to organize eqtl_pca_scree.R in the pipeline with 
# define_treatments_and_controls_pergene but it does not work as the PCA calculating time is too long
# then spliting this function to two parts: precalculating PCA and get the dimenstion of PC
fit_bespoke = function(gene_set_name, p_eqtl){
define_treatments_and_controls_pergene(gene_set_name, p_eqtl)
custom_PCA <- readRDS(str_c("/home/share/dna_ancestry/dna/custom_PCA_", gene_set_name, "_", p_eqtl,".rds"))
custom_PCA = custom_PCA %>% select(-fid) %>% mutate(AID = AID %>% as.character())
recode_variables_in_dat_racedummy()
print(abbreviations)
funcs = str_subset(abbreviations$shorthand, "^m")
funcs = funcs %>% str_subset("m[7-8]")
comp = 10
fit_pca_util = partial(fit_pca_util, ncomp = ncomp) # specify n_perm
# debugonce(model_fit)
example0 =
args %>%
mutate(out = pmap(., safely(model_fit), funcs),
control_set = names(controls))
example1 =
args %>%
mutate(gene_set_name = "whole_genome_and_tfbm") %>%
mutate(out = pmap(., safely(model_fit), funcs),
control_set = names(controls))
return(list(example0 = example0, example1 = example1))
}



test = function(a,b,c){
  temp1= a+b
  temp1 %>% saveRDS("./user_wx/sum.rds")
  sum = readRDS("./user_wx/sum.rds")
  environment()
  return(c + d +sum)
  
}

test1 = function(a,b){
  d = a*b+c
  list2env(list(d=d),environment(test))
}
# debugonce(test)
test(1,23,3)


