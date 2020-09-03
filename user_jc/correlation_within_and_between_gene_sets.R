############################################################
# CALL AN EXAMPLE FILE FIRST TO ENSURE dat ETC ARE IN SCOPE
############################################################

############################################################
# CORRELATION WITHIN AND BETWEEN  GENE SETS
############################################################
 
# get the normalized mRNA for each gene set
vals = map(sigs, ~t(exprs(dat[.x])))

# calculate the average correlation within gene-set
vals %>% map_df(compose(mean, cor)) 

# calculate average cross correlation across gene sets (poor man's cannonical correlation or similar)
crossing(x = vals, y = vals) %>% 
  mutate(x = names(x), y = names(y), b = pmap_dbl(., compose(mean, cor))) %>%
  pivot_wider(values_from = b, names_from = x)