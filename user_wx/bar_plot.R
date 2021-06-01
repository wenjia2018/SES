# https://gist.github.com/Emaasit/bec9825315fae54964c4
# https://sebastiansauer.github.io/percentage_plot_ggplot2_V2/
# https://stackoverflow.com/questions/37008705/ggplot-bar-chart-of-percentages-over-groups
# https://aosmith.rbind.io/2020/07/09/ggplot2-override-aes/
# https://stackoverflow.com/questions/15059093/ggplot2-adjust-the-symbol-size-in-legends
# https://stackoverflow.com/questions/52341385/how-to-automatically-adjust-the-width-of-each-facet-for-facet-wrap/52422707
load_data(reconciled = FALSE, remove_inflam = FALSE)
pData(dat) <- pData(dat) %>%
  # mutate_at(.vars = vars(matches("income|composite|SEI")),
  #           .funs = list(~ .x %>% ntile(4) %>% as.factor())
  # ) %>%
  mutate(raceethnicity = re %>%
           fct_recode(
             NonHwhite = "1",
             # white nonhispanic
             NonHblack = "2",
             # black nonhispanic
             NULL = "3",
             # asian nonhispanic
             NULL = "4",
             # other nonhispanic
             Hispanic = "5"
             # hispanic
           ) %>%
           relevel(ref = "NonHwhite"),
         color_byinterviewer3 = H3IR17 %>%
           as.character() %>% 
           as.factor %>% 
           fct_collapse(
             DarkBlack = c("1", "2"),
             LightMed = c("3", "4"),
             White = "5") %>%
           relevel(ref = "White"))

keep = (!is.na(pData(dat)$color_byinterviewer3) & !is.na(pData(dat)$raceethnicity)) %>% ifelse(is.na(.), FALSE, .)
dat <- dat[, keep]



ggplot(data = data, aes(x = color_byinterviewer3,  fill = raceethnicity)) + 
  geom_bar(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..]), position="dodge") +
  geom_text(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..], label=scales::percent(..count../tapply(..count.., ..x.. ,sum)[..x..]) ),
            stat="count", position=position_dodge(0.9), vjust=1) +
  geom_text(stat='count', aes(y=..count../tapply(..count.., ..x.. ,sum)[..x..], label=..count..), position=position_dodge(0.9), vjust=-0.5)+
  ylab('Percent of Population, %') +
  scale_y_continuous(labels = scales::percent)



# ggplot(data=mtcars, aes(cyl))+
#   geom_bar(aes(fill=as.factor(gear), y = (..count..)/sum(..count..)), 
#            position="dodge")+
#   geom_text(aes( y=..count../(..count..)/sum(..count..), label=scales::percent((..count..)/sum(..count..))),
#             stat="count", position=position_dodge(0.9), vjust=1)
# 
# ggplot(mtcars, aes(x=as.factor(cyl), fill=as.factor(gear)))+
#   geom_bar(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..]), position="dodge" ) +
# geom_text(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..], label=scales::percent(..count../tapply(..count.., ..x.. ,sum)[..x..]) ),
#           stat="count", position=position_dodge(0.9), vjust=1)
