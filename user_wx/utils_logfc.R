isnt_out_z <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - mean(x, na.rm = na.rm)) <= thres * sd(x, na.rm = na.rm)
}

enf <-function(trimmedLogFCs){
  r = seq(0.1, 0.9, by = 0.1)
  n = length(trimmedLogFCs)
  mu = mean(trimmedLogFCs)
  stdv = sd(trimmedLogFCs)
  t = map_dbl(r, function(rho) mu/sqrt(stdv^2 * (1/n + rho/(1- rho))))
  thr_left = qt(0.025, n-1)
  thr_right = qt(0.975, n-1)
  
  p = enframe(t, value = "t.statistic") %>% 
    transmute(cor = r,
              # sig = ifelse(between(t.statistic, thr_left, thr_right), "Non Sig", "Sig"),
              p = 2*pt(-abs(t.statistic), df = n-1))
}

DE_logFC_ploting = function(example, caption_text, cor_coef=0.1){
 a = example %>% 
   mutate(dim = map_int(.$m, ~dim(.)[1]),
          gene_set_name = str_c(gene_set_name, "(", dim, ")")) %>% 
    mutate_at(.vars = vars("m"),
              .funs = list(~ map(., ~ mutate(., logFC_outlier = isnt_out_z(logFC)) %>% filter(logFC_outlier==TRUE)))) %>% 
    # mutate(test = m %>% map(~ wilcox.test(.$logFC, mu = 0, alternative = "two.sided")),
    #        logFC.t.test.p = test %>% map_dbl(~ .$p.value) %>% p.adjust(method = "fdr") 
    #        # %>% format(digits = 3, scientific =T)
    # ) %>% 
    mutate(test_statistics = m %>% map(~ enf(.$logFC))) %>% 
    unnest(test_statistics) %>% 
    group_by(cor) %>% 
    mutate(p.adj = p.adjust(p,method= "fdr")) %>% 
    ungroup %>% 
    unnest(m) %>% 
    filter(cor==cor_coef) %>% 
    dplyr::select(treatment, gene_set_name, logFC, p.adj) 
   
  axiscolor = c("grey30", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30","grey30", "grey30", "grey30")
  a %>% 
    mutate(p_sig = ifelse(p.adj>0.05, "Non Sig", "Sig")) %>% 
    ggplot(aes( x = logFC, y = gene_set_name,
                color = p_sig,
                fill = p_sig)) +
    geom_violin()+
    facet_wrap( ~  treatment , scales = "free_x", strip.position = "bottom")  +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
    scale_color_manual(values=c("goldenrod3", "lightblue"))+
    scale_fill_manual(values=c("goldenrod3", "lightblue"))+
    theme(
      # legend.position = "none",
      panel.margin.x = unit(1, "lines"),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(color=axiscolor),
      axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
      axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
      plot.title = element_text(size = 20, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
      plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
      strip.text = element_text(size = 7),
      text = element_text(family = "Georgia"))+
    labs(caption = paste("correlation = ", cor_coef, ",", caption_text))
  
}

DE_logFC_ploting_notest = function(example, caption_text){
  a = example %>% 
    # mutate(dim = map_int(.$m, ~dim(.)[1]),
    #        gene_set_name = str_c(gene_set_name, "(", dim, ")")) %>% 
    ungroup %>% 
    unnest(m) %>% 
    dplyr::select(treatment, gene_set_name, logFC) 
  
  axiscolor = c("grey30", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30","grey30", "grey30", "grey30")
  a %>% 
    ggplot(aes(x = logFC, y = gene_set_name)) +
    geom_violin(fill = "lightblue", color ="lightblue")+
    facet_wrap( ~  treatment , scales = "free_x", strip.position = "bottom") +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    scale_color_manual(values=c("goldenrod3", "lightblue")) +
    scale_fill_manual(values=c("goldenrod3", "lightblue")) +
    theme(
      # legend.position = "none",
      panel.margin.x = unit(1, "lines"),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(color=axiscolor),
      axis.title.y = element_text(margin = margin(r = 20), color = "grey70"),
      axis.title.x = element_text(margin = margin(t = 20), color = "darkslategrey"),
      plot.title = element_text(size = 20, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 7, color = "darkslategrey", margin = margin(b = 25)),
      plot.caption = element_text(size = 7, margin = margin(t = 10), color = "grey70", hjust = 0),
      strip.text = element_text(size = 7),
      text = element_text(family = "Georgia"))+
    labs(y = "mRNA Signatures",
         caption = paste(caption_text))
  
}
