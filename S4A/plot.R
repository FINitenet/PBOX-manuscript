library(openxlsx)
library(tidyverse)
library(ggh4x)
library(ggsci)
library(ggbeeswarm)

df <- read.xlsx("FigS3A.xlsx")

df$group <- as.factor(df$group)
df$X1 <- factor(df$X1, levels = c("inflorescences","seedlings"))
df$X2 <- factor(df$X2, levels = c("Public sRNA-seq data","Conventional sRNA-seq","PBOX-sRNA-seq"))




p <- ggplot(df) +
 aes(x = X1, y = count, colour = group) +
 geom_beeswarm(size = 2, cex = 5, show.legend = FALSE) +
 scale_color_igv() +
 scale_y_continuous(limits = c(0,100),breaks = seq(0,100,20),expand = c(0,0))+
 stat_summary(fun = median, 
               geom = 'crossbar', width = 0.5, linewidth = 0.1, color = 'black') +
 theme_classic()+
 facet_nested(~ X2 + X1,
                 scales = "free",
                 switch = "x")+
  theme(ggh4x.facet.nestline = element_line(colour = "black"),
        plot.subtitle = element_text(size = 9),
    plot.caption = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    plot.title = element_text(size = 9),
    legend.position = "none",
    strip.background = element_rect(
      color = "white", fill = "white"),
    panel.grid = element_blank(),
    strip.text.x = element_text(size = 9)) +
  labs(x = NULL, y = "Relative percentage(%)")
p
ggsave("FigS3A.pdf",p, width = 7,height = 6, units = "cm")


df %>% filter(X2 == "Public sRNA-seq data" & X1 == "inflorescences") %>% 
ggplot() +
  aes(x = X1, y = count, colour = group) +
  geom_beeswarm(size = 1.5, cex = 5, show.legend = FALSE) +
  scale_y_continuous(limits = c(0,100),breaks = seq(0,100,20),expand = c(0,0))+
  stat_summary(fun = median, 
               geom = 'crossbar', width = 0.5, linewidth = 0.1, color = 'black') +
  theme_classic()+
  theme(ggh4x.facet.nestline = element_line(colour = "black"),
        plot.subtitle = element_text(size = 9),
        plot.caption = element_text(size = 9),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        plot.title = element_text(size = 9),
        legend.position = "none",
        strip.background = element_rect(
          color = "white", fill = "white"),
        panel.grid = element_blank(),
        strip.text.x = element_text(size = 9)) +
  labs(x = NULL, y = "Relative percentage(%)")










