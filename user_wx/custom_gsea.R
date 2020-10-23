# create own genesets for GSEA

sets = signatures$outcome_set[table1[1:11]]
sets$CVD_mRNA %>% as.data.frame() %>% t
a = plyr::ldply(sets, rbind) %>% add_column(link = "//", .after = ".id")
write.table(a, file ="./user_wx/disease_sets.gmt", quote = FALSE,
            col.names = FALSE, row.names = FALSE, na = "",
            sep = "\t")


# http://www.webgestalt.org/results/1603207340/#