
#' ---
#' title: creating new signatures from DE geno wide regression
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' <!-- rmarkdown::render("supervised_play/nice_code.R") -->
#' <!-- [See here.](http://brooksandrew.github.io/simpleblog/articles/render-reports-directly-from-R-scripts/) -->
#+ warning=FALSE, message=FALSE


#+ echo=F, eval=T, warning=FALSE, message=FALSE
source("/home/xu/ses-1/user_wx/extract_v2.R")
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

control = "ancestryPC_ses"
p_eqtl <- c(0.05, 1e-2)

#' ### no sig genes for race 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
race <- readRDS("~/ses-1/user_wx/race_bespoke_DE_01.03.2021.rds")
a = p_eqtl %>% map(outttT, control, race)
a

#' ### sig genes for skin color of darkblack 
#+ echo=F, eval=T, warning=FALSE, message=FALSE
skincolor <- readRDS("~/ses-1/user_wx/skincolor3_bespoke_DE_01.03.2021.rds")



darkblack005 = outttT(0.05, control, skincolor)$gene_sig[[1]]
darkblack001 = outttT(0.01, control, skincolor)$gene_sig[[1]]
darkblack005intersec001 = intersect(darkblack001, darkblack005)


venn::venn(list(darkblack005 = darkblack005, darkblack001 = darkblack001),
           ilabels = TRUE, 
           zcolor = "style", ellipse = FALSE,
           opacity = 0.15)

signatures$outcome_set$darkblack005_mRNA = darkblack005
signatures$outcome_set$darkblack001_mRNA = darkblack001
signatures$outcome_set$darkblack005intersec001_mRNA = darkblack005intersec001

# gsea method
#' ### gsea pysiological function
#+ echo=F, eval=T, warning=FALSE, message=FALSE
if(0){
  tempfolder = "temp_webgestalt"
 
   example_race = outttT(0.05, control, race) %>% 
    rbind(outttT(0.01, control, race))
 
   example_skincolor = outttT(0.05, control, skincolor) %>% 
    rbind(outttT(0.01, control, skincolor))
 
   example = example_race %>% 
    rbind(example_skincolor) %>% 
    mutate(treatment = str_c(treatment, "_", p_eqtl)) %>% 
    select(treatment, ttT) 
 
   args = 
    example %>% 
    mutate(file_output = str_c(getwd(), "/", tempfolder))
  
  complete_tables = args %>%
    pmap(gsea_webgestalt) %>%
    set_names(args$treatment)
  
  
  complete_tables %>% saveRDS("/home/xu/ses-1/user_wx/race_skincolor_gsea.rds")
  
}

complete_tables = readRDS("/home/xu/ses-1/user_wx/race_skincolor_gsea.rds")

complete_tables = complete_tables %>% 
  map(~ .x %>% filter(FDR < 0.05))


a = map2(complete_tables, names(complete_tables), ~ mutate(.x, treatment = .y)) %>% 
  bind_rows() %>% 
  select(treatment, everything()) %>% 
  select(treatment, description, FDR, userId)
a %>% kableExtra::kable() %>% kableExtra::kable_styling()