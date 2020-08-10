celltype_cibersort <- function(treatment, controls) {

  cell_types <- c(
    "B.cells.naive", "B.cells.memory", "Plasma.cells",
    "T.cells.CD8", "T.cells.CD4.naive", "T.cells.CD4.memory.resting",
    "T.cells.CD4.memory.activated",
    "T.cells.follicular.helper",
    "T.cells.regulatory..Tregs.", "T.cells.gamma.delta",
    "NK.cells.resting", "NK.cells.activated", "Monocytes", "Macrophages.M0",
    "Macrophages.M1",
    "Macrophages.M2", "Dendritic.cells.resting",
    "Dendritic.cells.activated", "Mast.cells.resting",
    # "Mast.cells.activated", # all zeros
    "Eosinophils", "Neutrophils"
  )

  term <- c(controls, treatment)
  data <- pData(dat) %>% dplyr::select(all_of(term)) %>% rename(treatment = treatment)
  rhs <- str_c(data %>% colnames(), collapse = " + ")
  keep <- data %>%
    complete.cases() # keep only the complete cases
  nsubs <- sum(keep) # number of complete cases

  ########################################################
  # DEFINE DEPENDENT VARIABLE SETS AND LOAD INTO THE GLOBAL ENV
  ########################################################

  outcomes <- pData(dat) %>%
    dplyr::select(all_of(cell_types)) %>%
    compositions::acomp() %>%
    compositions::ilr()
  outcomes = outcomes[keep, ]

  m <- lm(str_c("outcomes ~ ", rhs) %>% as.formula(), data = data[keep, ])
  Manova <- extract_m1(m)
  # partial effect
  include_list <- coef(m) %>% rownames() %>% str_subset("treatment")
  b <- coef(m)[include_list, ] %>% compositions::ilrInv()
  names(b) <- cell_types

  out <- list(b = b, Manova = Manova)
  
  return(out)
}
