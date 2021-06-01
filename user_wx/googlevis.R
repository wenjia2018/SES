###############
# Bubble Plot #
###############
library(data.table)
library(ggplot2)
library(ggrepel)
bubble_data <-fread("https://raw.githubusercontent.com/zhendata/Medium_Posts/c007346db1575aca391a6623c87bb5a31a60b365/bubble_plot_merged_city_data.csv",sep=",")
bubble_plot <- ggplot(bubble_data, aes(x = Unemployment_Rate, y = Home_Price/1000)) + 
  geom_point(aes(size = Population, fill = Total_Crime),shape=21) +
  # Create 'Bubble' by assigning size a variable #
  scale_fill_continuous(low = "#33FFFF", high ="#FF6699" ) +
  scale_size_area(max_size = 20)+
  # Select bubble color scale and bubble maximum size #
  geom_text_repel(aes(label = City),nudge_x = 0,nudge_y = 0.75,size = 6) +
  # Use geom_text_repel to repel the labels away from each other #
  theme_bw()+
  # Use white background instead of the default grey one #
  ggtitle("Best Cities in US to Live in") +
  labs(x = "Unemployment Rate%", y = "Home Price",
       size = "Population",fill="Crime") +
  theme(plot.title = element_text(size=25, hjust = 0.5),
        axis.title=element_text(size=20, face = "bold"),
        axis.text=element_text(size=15)) +
  # Style title and axis #
  scale_y_continuous(name="Home Price", breaks = seq(0, 1500, by=250), 
                     labels=c("0", "250K", "500K", "750K", "1000k",    "1250k", "1500K"))
# Make y-axis more readable by replacing scientific number by "K" #
print(bubble_plot)
################
# Motion Chart #
################
library(data.table)
library(googleVis)
motion_data <-fread("https://raw.githubusercontent.com/zhendata/Medium_Posts/c007346db1575aca391a6623c87bb5a31a60b365/motion_chart_merged_city_data.csv",sep=",")
motion_chart <- gvisMotionChart(motion_data,idvar = "City", timevar = "Year",    xvar = "Unemployment Rate", sizevar="Population")
plot(motion_chart)
# R automatically opens a tab in the browser for you
# The flash player needs to be enabled in browser