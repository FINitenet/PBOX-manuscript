library(tidyverse)
library(ggh4x)
library(lemon)
library(patchwork)


meta <- read.table("/bios-store1/chenyc/Project/TRM-sRNA-seq/QC/2_srna_analysis/02.Mapping_18_28/mapping_results_ShortStack.txt", sep = '\t', header = TRUE)

df1 <- data.frame()
files1 <- list.files(path = "doc/",pattern = "*20200818N_lncRNA_seq_len_5bias.csv")
for (file in files1) {
  sample <- gsub("_lncRNA_seq_len_5bias.csv", "", file)
  tmp <- read.csv(paste0("doc/", file))
  tmp$RPM <- tmp$count/meta[meta$Sample == sample, 3]*1000000
  tmp <- tmp[,-3]
  df1 <- rbind(df1, tmp)
}

df2 <- data.frame()
files2 <- list.files(path = "doc/",pattern = "*20201010N_lncRNA_seq_len_5bias.csv")
for (file in files2) {
  sample <- gsub("_lncRNA_seq_len_5bias.csv", "", file)
  tmp <- read.csv(paste0("doc/", file))
  tmp$RPM <- tmp$count/meta[meta$Sample == sample, 3]*1000000
  tmp <- tmp[,-3]
  df2 <- rbind(df2, tmp)
}

df3 <- rbind(df1,df2)

df3 <- df3 %>% mutate(
  rep = case_when(
    grepl("20200818N", library) ~ "rep1",
    grepl("20201010N", library) ~ "rep2",
    TRUE ~ "others"
  ),
  method = case_when(
    grepl("LCK_|LerCK_", library) ~ "Ler CK",
    grepl("LOX_|LerOX_", library) ~ "Ler OX",
    grepl("LAP_|LerAP_", library) ~ "Ler PBA",
    grepl("LOXAP_|LerOXAP_", library) ~ "Ler PBOX+OX",
    
    grepl("hCK_|henCK_", library) ~ "hen CK",
    grepl("hOX_|henOX_", library) ~ "hen OX",
    grepl("hAP_|henAP_", library) ~ "hen PBA",
    grepl("hOXAP_|henOXAP_", library) ~ "hen PBOX+OX",
    TRUE ~ "others"
  )
)

df4 <- df3 %>% 
  group_by(length, bias, method) %>%
  summarize(avg_count = mean(RPM, na.rm = TRUE), .groups = "drop")


lnc <- ggplot(df4) +
  aes(x = length, y = avg_count, fill = bias) +
  geom_col() +
  scale_fill_brewer(palette = "Accent", direction = 1) +
  scale_x_continuous(expand = c(0,0.5),breaks = seq(18,30,1))+
  theme_bw() +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 9),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        strip.background = element_blank(),
        strip.text = element_text(size = 9))+
  labs(title = "lncRNA")+
  facet_rep_wrap(~ method, scales='free_x',nrow = 2)


ggsave("S1L.pdf",lnc,width = 8, height = 4)


