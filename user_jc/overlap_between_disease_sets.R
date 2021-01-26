#' ---
#' title: Examples
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#' Set global options 
#+ setup, warning=FALSE, message=FALSE
# knitr::opts_chunk$set(echo = FALSE)

set.seed(123)
library(here)
library(tidyverse) 
library(Biobase) 

walk(dir(path = here("R"), full.names = TRUE), source) 

############################################################
# LOAD DATA, DEFINE VARIABLES, RECODE VARIABLES
############################################################

load_data(reconciled = FALSE, remove_inflam = FALSE)
define_treatments_and_controls()

sigs = signatures$outcome_set[c(table1, "ctra_mRNA")]

get_intersection  = 
  function(intersection_order, op = intersect) {
    # GET ALL INTERSECTIONS OF A GIVEN ORDER
    # y = combinations of sets
    # nm = names of such sets
    # z = intersection of such sets
    
    tibble(y = combn(sigs, intersection_order, simplify = FALSE),
           nm = map(y, names),
           z =  map(y, reduce, op)) %>%
      arrange(-lengths(z)) %>% 
      filter(lengths(z) > 0)
  }

prettify = . %>% 
  select(-y) %>% 
  knitr::kable()

# overlap between the disease sets and the 1k.
a = sigs[-3] # less 1k 
b = sigs[3] # 1k
crossing(a,b) %>% 
  pmap_dbl(function(a,b) length(intersect(a,b))/length(a)) %>% 
  enframe(name = "signature", value = "proportion which overlaps with 1k")

# how many n-way intersections?
map(2:5, compose(prettify, get_intersection))

# Inflammation1k intersects with all other disease sets, by at least one gene (rarely the same gene)
map(2:5, get_intersection) %>% map(pluck("nm")) %>% map(unlist) %>% map(table) %>% map(sort) %>% map(enframe)

# Only 1 gene belongs to as many as 4 signatures (TCF7). 

# non-monogomous genes (appearing in at 2 or more signatures)
map(2:5, get_intersection) %>% map(pluck, "z") %>% map(unlist) %>% map(unique) %>% lengths()

# how many n-way intersections does each gene belong to
map(2:5, get_intersection) %>% 
  map(pluck, "z") %>%
  map(compose(fct_count, factor, unlist)) %>%
  map(arrange, -n) %>% 
  walk(print, n = Inf)

# GOC "genes of co-morbidity". 
sigs %>% lengths()
ns = lengths(sigs)
map(2:5, get_intersection) %>% map(transmute, nm, n = lengths(z)) %>% map(unnest_wider, "nm")

# size of intersection and the size of the disease sets
map(2, get_intersection) %>% 
  map(transmute, nm, n = lengths(z)) %>% 
  map(unnest_wider, "nm") %>% 
  map(mutate, n1 = ns[`...1`], n2 = ns[`...2`]) %>% 
  map(mutate, p1 = n/n1, p2 = n/n2) %>%
  map(print, n = Inf)
