# Fig1C ------------------------------
df <- read.table("Fig1/summary_18-28nt.txt",header = TRUE)
df_long <- gather(df, sample, value, names(df)[2]:names(df)[17])
df_long$Geneid <- gsub("NoFeatures|otherRNA", "unassigned", df_long$Geneid)
df_long$Geneid <- gsub("miRNA_primary_transcript", "miRNA", df_long$Geneid)
df_long <- aggregate(df_long$value, by=list(sample = df_long$sample, type = df_long$Geneid),sum)
df_long <- df_long %>% mutate(
    rep = case_when(
      grepl("20200818N$", sample) ~ "rep1",
      grepl("20201010N$", sample) ~ "rep2",
      TRUE ~ "others"
    ),
    method = case_when(
      grepl("LCK_|LerCK_", sample) ~ "Ler CK",
      grepl("LOX_|LerOX_", sample) ~ "Ler OX",
      grepl("LAP_|LerAP_", sample) ~ "Ler PBA",
      grepl("LOXAP_|LerOXAP_", sample) ~ "Ler PBOX+OX",
      
      grepl("hCK_|henCK_", sample) ~ "hen CK",
      grepl("hOX_|henOX_", sample) ~ "hen OX",
      grepl("hAP_|henAP_", sample) ~ "hen PBA",
      grepl("hOXAP_|henOXAP_", sample) ~ "hen PBOX+OX",
      
      TRUE ~ "others"
    )
  )

df_long$rep <- factor(df_long$rep, levels = c("rep1","rep2"))
df_long$type <- factor(df_long$type, levels = rnatype_level)
df_long$method <- factor(df_long$method, levels = method_level1)
df_long$sample <- factor(df_long$sample, levels = c(
  "LCK_20200818N","LerCK_20201010N", "LAP_20200818N", "LerAP_20201010N","LOX_20200818N","LerOX_20201010N", "LOXAP_20200818N", "LerOXAP_20201010N",
  "hCK_20200818N","henCK_20201010N", "hAP_20200818N", "henAP_20201010N","hOX_20200818N","henOX_20201010N", "hOXAP_20200818N", "henOXAP_20201010N"))

p <- ggplot(df_long) +
  aes(x = sample, fill = type, weight = x) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = rnatype_colour) +
  theme_classic()+
  labs(x = NULL,y = "Relative percentage(%)",fill = "RNA type")+
  scale_y_continuous(labels = scales::percent_format(suffix = ""),breaks = seq(0,1,0.1),expand = c(0,0),position = "right")+
  scale_x_discrete(labels = rep(c("CK_rep1","CK_rep2", "PBA_rep1","PBA_rep2", "OX_rep1","OX_rep2", "PBA\n+OX\nrep1","PBA\n+OX\nrep2"),2),expand = c(0,0.5))+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        #axis.line.x = element_line(size = 1),#轴线粗细
        axis.text.y = element_text(size = 12),
        #axis.line.y = element_line(size = 1),#轴线粗细
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))
ggsave("Fig1/fig1c_v1.pdf",p, width = 10,height = 6)
