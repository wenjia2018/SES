ta <- example0_with1k %>%
  hoist(out, p = list("result", "m8_fdr", 1, "p")) %>% 
  hoist(out, t = list("result", "m8_fdr", 1, "detail","t")) %>% 
  filter(controls=="all") %>% 
  dplyr::filter(treatment %in% c("ses_sss_composite", "edu_max", "income_hh_ff5", "SEI_ff5",  "sss_5")) %>% 
  dplyr::select(treatment, gene_set_name, p,t)

par(mfrow=c(1,2))

g1= ggplot(mtcars
       ,aes(x=hp, y=disp, size=mpg, fill=factor(gear)))+
  geom_point(shape=21)+
  theme_bw()+
  scale_size_continuous(range=c(2,15))+
  guides(shape = guide_legend(override.aes = list(size = 10)))

g2 =ggplot(mtcars
       ,aes(x=hp, y=disp, size=mpg, fill=factor(gear)))+
  geom_point(shape=21)+
  theme_bw()+
  scale_size_continuous(range=c(2,15))

g3 = ggplot(mtcars
            ,aes(x=hp, y=disp, size=mpg, fill=factor(gear)))+
  geom_point(shape=21)+
  theme_bw()+
  scale_size_continuous(range=c(2,15))+
  guides(shape = guide_legend(override.aes = list(size = 10)),
         fill = guide_legend(override.aes = list(size = 10)))
ggarrange(g1,g3)

xy.df = mtcars
xy.list_row <- as.list(as.data.frame(t(xy.df)))
xy.list_row <- as.list(as.data.frame(t(xy.df)))
