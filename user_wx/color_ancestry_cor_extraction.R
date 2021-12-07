#+ echo=F, eval=T, warning=FALSE, message=FALSE,comment=NA
# extraction of the results

color_ancestry_cor <- readRDS("~/ses-1/user_wx/color_ancestry_cor.rds")

temp = 
  color_ancestry_cor %>% 
  hoist(out, anova_ordinal = list("result", "anova_ordinal")) %>% 
  hoist(out, anova_unordered = list("result", "anova_unordered")) %>% 
  hoist(out, parallel_test = list("result", "parallel_test")) %>% 
  hoist(out, cor = list("result", "cor_results")) 

for (i in 1 : dim(temp)[1]) {
  
  cat(" ############################################################","\n",
      "gene_set_name is", temp[i, ]$table1,"\n",
      "############################################################","\n")
  cat("correlation","\n","\n")
  
  print(temp[i,]$cor)
  cat("\n","\n")
  
}

for (i in 1 : dim(temp)[1]) {
  
  cat(" ############################################################","\n",
      "gene_set_name is", temp[i, ]$table1,"\n",
      "############################################################","\n")
  cat("correlation","\n","\n")
  
  print(temp[i,]$cor)
  cat("\n","\n")
  
  cat("ordered categorical regression Parallel assumption test ","\n","\n")
  print(temp[i,]$parallel_test)
  cat("ordered categorical regression anova","\n","\n")
  print(temp[i,]$anova_ordinal)
  cat("unordered categorical regression anova","\n","\n")
  print(temp[i,]$anova_unordered)
  cat("\n","\n")
  
}


