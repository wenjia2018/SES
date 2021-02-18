devtools::install_github("ctesta01/QualtricsTools")
devtools::install_github("datastorm-open/DependenciesGraphs")

library(DependenciesGraphs)
library(QualtricsTools)

library(dbr)
"package:dbr" %>%
  envirDependencies() %>%
  plot()
