
# check the missingness in a database dt
b = a %>% map(~ which(is.na(.x), arr.ind=TRUE))
cc = b[lengths(b) != 0]
cc %>% venn::venn(zcolor = "style", ellipse = FALSE, opacity = 0.15, ilcs = 1, box = FALSE, sncs = 1)
