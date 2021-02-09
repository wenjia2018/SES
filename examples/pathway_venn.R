############################################################
# after obtained analysis results containing mm8_fdr
# pathway presenting
# please specify the example object and the signature you
# want to look at
############################################################
example = readRDS("/home/share/projects/aging/example0_aging_clusters_pc5.rds")
gene_signature = "aging_mRNA"
# create a folder in your current directory to save all the intermediate files and report generated at last from
# webgestalt, you could name your new folder differently
tempfolder = "temp_webgestalt"
dir.create(tempfolder)

# create input for the gsea_webgestalt function
# ttT is the results extracted from m8_fdr
args = 
  example %>% 
  hoist(out, ttT = list("result", "m8_fdr", 1, "other", "m")) %>% 
  filter(gene_set_name == gene_signature) %>% # for a specific signature 
  select(treatment, ttT) %>% 
  mutate(file_output = str_c(getwd(), "/", tempfolder))

# function gsea_webgestalt calculate functional pathway using preferred database:
# pathway_KEGG and pathway_Reactome for example
# you could change the argument in the function 
gsea_webgestalt = function(treatment, ttT, file_output, enrichMethod = "GSEA", enrichDatabase="pathway_Reactome"){
  rankFile = str_c(file_output,"/", treatment,".rnk")
  
  ttT %>%
    dplyr::select(gene, logFC) %>% 
    write.table(file = rankFile, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  
  enrichResult = WebGestaltR::WebGestaltR(interestGeneFile = rankFile,
                                          interestGeneType = "genesymbol",
                                          enrichMethod = enrichMethod,
                                          organism = "hsapiens",
                                          enrichDatabase = enrichDatabase,
                                          fdrThr = 1,
                                          minNum=5,
                                          perNum = 1000,
                                          outputDirectory = file_output)
}

# get pathway for FRD < 0.05
# be noted that the results are slightly different every time you run it, as it is based on permutations
# you could change perNum to enlarge the number of permutation
complete_tables = args %>%
  pmap(gsea_webgestalt) %>% 
  map(~ .x %>% filter(FDR < 0.05)) %>%
  set_names(args$treatment)

# extract sets of pathway for each treatment
gsea_genesetnames = 
  complete_tables %>%
  map(~ .$geneSet)

# union of all the pathway
complete_tables =
  complete_tables %>%
  reduce(rbind) %>%
  select(geneSet, description) %>% 
  unique()

# calculate the complement of the pathway for each treatment
gsea_genesetnames_complement = map(gsea_genesetnames, ~setdiff(complete_tables$geneSet, .x)) # the complement of these reactome terms

# the universe of reactome terms considered here is "all reactome sets deemed significantly related to at least one ses predictor"
# there are then 2^5 intersections in the venn diagram (i.e. intersections over the 5 ses predictors, with each being the complement or not)

universe= 
  list(gsea_genesetnames_complement,
       gsea_genesetnames) %>% 
  transpose() 

get_venn_cell = 
  function(P){ 
    
    # Get the intersection for each collection of "literals" of ABCDE (a "literal"
    # of A is either A or the complement of A, i.e. A') for each row in matrix
    # of indexes P. These fill the 2^D cells of the Ven diagram.
    map(1:dim(P)[1], 
        # pluck the corresponding subset 
        ~ map2(names(universe),
               P[.x, ], 
               ~ pluck(universe, .x, .y + 1)
        ) 
    ) %>% 
      # interestion over all
      map(reduce, intersect)
  }

crossing(!!!set_names(rerun(length(universe), 0:1), 
                      names(universe))) %>% 
  mutate(geneSet = get_venn_cell(.)) %>% 
  unnest(geneSet) %>% 
  left_join(complete_tables, by = "geneSet")  %>% 
  knitr::kable()



venn::venn(gsea_genesetnames, ilabels = TRUE, 
           zcolor = "style", ellipse = FALSE,
           opacity = 0.15)