# setting  ------------------------------
setwd("D:/^^Project/TRM-sRNA-seq/publish results/Data Source")
options(scipen = 999)
library(tidyverse)
library(patchwork)

rnatype_colour <- c(unassigned = "#a6cee3",
                    snoRNA = "#33a02c",snRNA = "#fb9a99",rRNA = "#1f78b4",tRNA = "#b2df8a",
                    lncRNA = "#e31a1c",protein_coding = "#fdbf6f",
                    transposable_element = "#ff7f00",`hc-siRNA` = "#cab2d6",`phasi/tasiRNA` = "#6a3d9a",miRNA = "#ffff99")
rnatype_level <- c("unassigned",
                   "snoRNA","snRNA","rRNA","tRNA",
                   "lncRNA","protein_coding",
                   "transposable_element","hc-siRNA","phasi/tasiRNA","miRNA")

method_level <- c("Conventional sRNA library","PBA-sRNA-library","OX-sRNA-library", "PBOX-sRNA-library")

method_level1 <- c("Ler CK","Ler PBA","Ler OX","Ler PBA+OX",
                   "hen CK", "hen PBA","hen OX","hen PBA+OX")


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

p1 <- ggplot(df_long) +
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
ggsave("Fig1/fig1c_v1.pdf",p1, width = 10,height = 6)


p2 <- ggplot(df_long) +
  aes(x = method, fill = type, colour = rep, weight = x) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = rnatype_colour)+
  scale_color_hue(direction = 1) +
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
ggsave("Fig1/fig1c_v2.pdf",p2, width = 8,height = 6)


# Fig2D ------------------------------

S2.list <- list.files(path = "Fig2D/",pattern = "*.summary")

for (i in S2.list) {
  sample <- gsub(".summary","",i)
  tmp <- read.table(paste0("Fig2D/",i))
  tmp$V2 <- gsub("others","unassigned",tmp$V2)
  tmp$V2 <- gsub("miRNA_primary_transcript","miRNA",tmp$V2)
  tmp$id <- paste(tmp$V1, tmp$V2, sep = "-")
  tmp$V2 <- factor(tmp$V2,levels = rnatype_level)
  assign(sample, tmp)
}

result_matrix <- sapply(rnatype_level, function(x) paste(14:35, x, sep = "-"))
result_vector <- as.vector(result_matrix)
meta <- data.frame(id = result_vector)

rep1 <- list.files(path = "Fig2D/", pattern = "*_1.summary")
rep2 <- list.files(path = "Fig2D/", pattern = "*2.summary")

meta1 <- meta
for ( i in rep1) {
  sample <- gsub(".summary","",i)
  tmp <- get(sample)
  tmp <- tmp %>% select(id, V4) 
  names(tmp)[2] <- sample
  meta1 <- full_join(meta1, tmp)
}

meta2 <- meta
for ( i in rep2) {
  sample <- gsub(".summary","",i)
  tmp <- get(sample)
  tmp <- tmp %>% select(id, V4) 
  names(tmp)[2] <- sample
  meta2 <- full_join(meta2, tmp)
}

cor_df <- data.frame()
# 循环比较两个数据框中的每对列，并计算相关性
for (col1 in colnames(meta1)[-1]) {
  for (col2 in colnames(meta2)[-1]) {
    # 计算相关性
    result <- cor.test(meta1[[col1]], meta2[[col2]])
    correlation <- round(result$estimate, 8)
    p_value <- round(result$p.value, 8)
    
    # 将相关性系数和p值添加到数据框中
    cor_df <- rbind(cor_df, c(col1, col2, correlation, p_value))
  }
}

# 重新命名列名
colnames(cor_df) <- c("Column_1", "Column_2", "Correlation", "pvalue")


for (col1 in colnames(meta1)[-1]) {
  col2 <- gsub("_1$","_2",col1)
  col <- gsub("_1$","",col1)
  meta[, col] <- (meta1[[col1]]+meta2[[col2]])/2
}


meta <- separate(meta, id, into = c("length", "type"), sep = "-")
meta$type <- gsub("hc$","hc-siRNA",meta$type)
meta$type <- factor(meta$type,levels = rnatype_level)
meta$length <- as.numeric(meta$length)


C1.P <- meta %>% select(1,2,CK_0) %>% rename("V1" = length, "V2" = type, "V4" = CK_0) %>% 
  ggplot() +
    aes(x = V1, fill = V2, weight = V4) +
    geom_bar() +
    scale_fill_manual(values = rnatype_colour) +
    theme_classic()+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          #axis.line.x = element_line(size = 1),#轴线粗细
          axis.text.y = element_text(size = 12),
          #axis.line.y = element_line(size = 1),#轴线粗细
          plot.title = element_text(size = 12,hjust = 0.5),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))+
    labs(x = "",y = "Relative percentage(%)",fill = "RNA type",title = "Conventional sRNA library")+
    annotate("text",label="R==0.9979~p<0.001", parse=T, x=32, y=55)+
    scale_x_continuous(breaks = seq(10, 35, 1),expand = c(0,0.5))+
    scale_y_continuous(breaks = seq(0,60,10), limits = c(0,60), expand = c(0,0.5))

T1.P <- meta %>% select(1,2,OXAP_0) %>% rename("V1" = length, "V2" = type, "V4" = OXAP_0) %>% 
  ggplot() +
  aes(x = V1, fill = V2, weight = V4) +
  geom_bar() +
  scale_fill_manual(values = rnatype_colour) +
  theme_classic()+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        #axis.line.x = element_line(size = 1),#轴线粗细
        axis.text.y = element_text(size = 12),
        #axis.line.y = element_line(size = 1),#轴线粗细
        plot.title = element_text(size = 12,hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+
  labs(x = "",y = "",fill = "RNA type",title = "TRMRNA-seq sRNA library")+
  annotate("text",label="R==0.9999~p<0.001", parse=T, x=32, y=55)+
  scale_x_continuous(breaks = seq(10, 35, 1),expand = c(0,0.5))+
  scale_y_continuous(breaks = seq(0,60,10), limits = c(0,60), expand = c(0,0.5))

C2.P <- meta %>% select(1,2,CK_10) %>% rename("V1" = length, "V2" = type, "V4" = CK_10) %>% 
  ggplot() +
  aes(x = V1, fill = V2, weight = V4) +
  geom_bar() +
  scale_fill_manual(values = rnatype_colour) +
  theme_classic()+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        #axis.line.x = element_line(size = 1),#轴线粗细
        axis.text.y = element_text(size = 12),
        #axis.line.y = element_line(size = 1),#轴线粗细
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+
  labs(x = "",y = "",fill = "RNA type")+
  annotate("text",label="R==0.9914~p<0.001", parse=T, x=32, y=55)+
  scale_x_continuous(breaks = seq(10, 35, 1),expand = c(0,0.5))+
  scale_y_continuous(breaks = seq(0,60,10), limits = c(0,60), expand = c(0,0.5))

C3.P <- meta %>% select(1,2,CK_15) %>% rename("V1" = length, "V2" = type, "V4" = CK_15) %>% 
  ggplot() +
  aes(x = V1, fill = V2, weight = V4) +
  geom_bar() +
  scale_fill_manual(values = rnatype_colour) +
  theme_classic()+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        #axis.line.x = element_line(size = 1),#轴线粗细
        axis.text.y = element_text(size = 12),
        #axis.line.y = element_line(size = 1),#轴线粗细
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+
  labs(x = "",y = "",fill = "RNA type")+
  annotate("text",label="R==0.9926~p<0.001", parse=T, x=32, y=55)+
  scale_x_continuous(breaks = seq(10, 35, 1),expand = c(0,0.5))+
  scale_y_continuous(breaks = seq(0,60,10), limits = c(0,60), expand = c(0,0.5))

T2.P <- meta %>% select(1,2,OXAP_10) %>% rename("V1" = length, "V2" = type, "V4" = OXAP_10) %>% 
  ggplot() +
  aes(x = V1, fill = V2, weight = V4) +
  geom_bar() +
  scale_fill_manual(values = rnatype_colour) +
  theme_classic()+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        #axis.line.x = element_line(size = 1),#轴线粗细
        axis.text.y = element_text(size = 12),
        #axis.line.y = element_line(size = 1),#轴线粗细
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+
  labs(x = "",y = "",fill = "RNA type")+
  annotate("text",label="R==0.9961~p<0.001", parse=T, x=32, y=55)+
  scale_x_continuous(breaks = seq(10, 35, 1),expand = c(0,0.5))+
  scale_y_continuous(breaks = seq(0,60,10), limits = c(0,60), expand = c(0,0.5))

T3.P <- meta %>% select(1,2,OXAP_15) %>% rename("V1" = length, "V2" = type, "V4" = OXAP_15) %>% 
  ggplot() +
  aes(x = V1, fill = V2, weight = V4) +
  geom_bar() +
  scale_fill_manual(values = rnatype_colour) +
  theme_classic()+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        #axis.line.x = element_line(size = 1),#轴线粗细
        axis.text.y = element_text(size = 12),
        #axis.line.y = element_line(size = 1),#轴线粗细
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+
  labs(x = "",y = "",fill = "RNA type")+
  annotate("text",label="R==0.9998~p<0.001", parse=T, x=32, y=55)+
  scale_x_continuous(breaks = seq(10, 35, 1),expand = c(0,0.5))+
  scale_y_continuous(breaks = seq(0,60,10), limits = c(0,60), expand = c(0,0.5))

S2 <- (C1.P+T1.P +grid::textGrob('T0: RIN=9.47',hjust = 2.5))/(C2.P+T2.P+grid::textGrob('T10: RIN=6.53',hjust = 2.5))/(C3.P+T3.P+grid::textGrob('T15: RIN=4.15',hjust = 2.5))+plot_layout(guides = 'collect')
S2
ggsave("Fig2D/fig2d_merge.pdf", S2, width = 17, height = 8)

# Fig3CD  ------------------------------
S3.list <- list.files(path = "Fig3CD/",pattern = "*.summary")

for (i in S3.list) {
  sample <- gsub(".summary","",i)
  tmp <- read.table(paste0("Fig3CD/",i))
  tmp$V2 <- gsub("otherRNA","unassigned",tmp$V2)
  tmp$V2 <- gsub("miRNA_primary_transcript","miRNA",tmp$V2)
  tmp$id <- paste(tmp$V1, tmp$V2, sep = "-")
  tmp$V2 <- factor(tmp$V2,levels = rnatype_level)
  assign(sample, tmp)
}

result_matrix <- sapply(rnatype_level, function(x) paste(10:35, x, sep = "-"))
result_vector <- as.vector(result_matrix)
meta <- data.frame(id = result_vector)

rep1 <- list.files(path = "Fig3CD/", pattern = "*R1.summary")
rep2 <- list.files(path = "Fig3CD/", pattern = "*R2.summary")

meta1 <- meta
for ( i in rep1) {
  sample <- gsub(".summary","",i)
  tmp <- get(sample)
  tmp <- tmp %>% select(id, V4) 
  names(tmp)[2] <- sample
  meta1 <- full_join(meta1, tmp)
}

meta2 <- meta
for ( i in rep2) {
  sample <- gsub(".summary","",i)
  tmp <- get(sample)
  tmp <- tmp %>% select(id, V4) 
  names(tmp)[2] <- sample
  meta2 <- full_join(meta2, tmp)
}

cor_df <- data.frame()
# 循环比较两个数据框中的每对列，并计算相关性
for (col1 in colnames(meta1)[-1]) {
  for (col2 in colnames(meta2)[-1]) {
    # 计算相关性
    result <- cor.test(meta1[[col1]], meta2[[col2]])
    correlation <- round(result$estimate, 8)
    p_value <- round(result$p.value, 8)
    
    # 将相关性系数和p值添加到数据框中
    cor_df <- rbind(cor_df, c(col1, col2, correlation, p_value))
  }
}

# 重新命名列名
colnames(cor_df) <- c("Column_1", "Column_2", "Correlation", "pvalue")


for (col1 in colnames(meta1)[-1]) {
  col2 <- gsub("R1$","R2",col1)
  col <- gsub("-R1$","",col1)
  meta[, col] <- (meta1[[col1]]+meta2[[col2]])/2
}


meta <- separate(meta, id, into = c("length", "type"), sep = "-")
meta$type <- gsub("hc$","hc-siRNA",meta$type)
meta$type <- factor(meta$type,levels = rnatype_level)
meta$length <- as.numeric(meta$length)



ck <- ggplot(meta) +
  aes(x = length, y = `IP-CK`, fill = type) +
  geom_col() +
  scale_fill_manual(values = rnatype_colour) +
  theme_classic()+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        #axis.line.x = element_line(size = 1),#轴线粗细
        axis.text.y = element_text(size = 12),
        #axis.line.y = element_line(size = 1),#轴线粗细
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.placement="outside",
        strip.background = element_rect(color = "white", fill = "white"),
        strip.text.x = element_text(size = 12))+
  labs(x = "Conventional AGO1-RIP sRNA library",y = "Relative percentage(%)",fill = "RNA type")+
  annotate("text",label="R==0.9969~p<0.001", parse=T, x=32, y=45)+
  scale_x_continuous(breaks = seq(10, 35, 1),expand = c(0,0.5))+
  scale_y_continuous(limits = c(0,50), breaks = seq(0,45,5),expand = c(0,0.5))

oxap <- ggplot(meta) +
  aes(x = length, y = `IP-OXAP`, fill = type) +
  geom_col() +
  scale_fill_manual(values = rnatype_colour) +
  theme_classic()+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        #axis.line.x = element_line(size = 1),#轴线粗细
        axis.text.y = element_text(size = 12),
        #axis.line.y = element_line(size = 1),#轴线粗细
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.placement="outside",
        strip.background = element_rect(color = "white", fill = "white"),
        strip.text.x = element_text(size = 12))+
  labs(x = "AGO1-RIP PBOX sRNA library",y = "Relative percentage(%)",fill = "RNA type")+
  annotate("text",label="R==0.98311~p<0.001", parse=T, x=32, y=45)+
  scale_x_continuous(breaks = seq(10, 35, 1),expand = c(0,0.5))+
  scale_y_continuous(limits = c(0,50), breaks = seq(0,45,5),expand = c(0,0.5))

S3 <- ck+oxap+ plot_annotation(tag_levels =  list(c('C', 'D')))+ plot_layout(guides = 'collect')
S3
ggsave("Fig3CD/fig3cd_merge.pdf",S3, width = 15, height = 4)

# Fig4E ------------------------------
buble_plot1 <- function(matrix.path, title, yanse) {
  df <- read.table(matrix.path, skip = 1)
  #df <- read.table("4_163.results/67_OXAP_R1_11x11/ath-MIR156a-3p_star.txt", skip = 1)
  sum <- sum(df)
  df <- sqrt(df / sum(df))*0.75 
  names(df) <- c(10:0)
  y0 <- c(10:0)
  df <- cbind(y0, df)
  df_long <- gather(df, x0, R, "10":"0")
  df_long$x0 <- as.integer(df_long$x0)
  
  ggplot(df_long) +
    geom_circle(aes(x0 = x0, y0 = y0, r = R), show.legend = FALSE, colour = NA, fill = yanse, na.rm = TRUE) +
    scale_x_reverse(limits = c(11, -1), breaks = seq(0, 10), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1, 11), breaks = seq(0, 10), expand = c(0, 0)) +
    # scale_size_area()+
    coord_fixed() +
    theme_bw() +
    theme(
      panel.grid.minor = element_line(colour = NA),
      panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
      panel.background = element_rect(fill = NA),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 9, face = "bold"),
      axis.text = element_text(size = 9, face = "bold", colour = "black"),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 7, face = "bold", hjust = 0.5),
      plot.background = element_rect(colour = NA),
      legend.position = "none") +
    xlab("Length of trimming") +
    ylab("Length of tailing") +
    ggtitle(title, paste0("total reads: ",sum))
}

for (i in seq_along(files)) { 
  miRNA_title <- gsub("\\.txt", "", files[i])
  miRNA_id <- gsub("_star", "\\*", miRNA_title)
  print(paste0("Processing(", i, "/", length(files), "):", miRNA_id))
  P1 <- buble_plot1(matrix.path = paste0("4_163.results/C1_178_240117N_L001/", files[i]), title = "C_rep1 TRM-sRNA-seq", yanse = "#9b4b8f")
  P2 <- buble_plot1(matrix.path = paste0("4_163.results/C2_178_240117N_L004/", files[i]), title = expression(bold(~bolditalic("C_rep2")~"TRM-sRNA-seq")), yanse = "#9b4b8f")
  P3 <- buble_plot1(matrix.path = paste0("4_163.results/h1_178_240117N_L001/", files[i]), title = expression(bold(~bolditalic("hntp_rep1")~"TRM-sRNA-seq")), yanse = "#9b4b8f")
  P4 <- buble_plot1(matrix.path = paste0("4_163.results/h2_178_240117N_L001/", files[i]), title = expression(bold(~bolditalic("hntp_rep2")~"TRM-sRNA-seq")),yanse = "#9b4b8f")
  merge <- (P1 + P2) / (P3 + P4) + plot_annotation(title = miRNA_id, theme = theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold")))
  # merge <- (P1 + P2 + P3) + plot_annotation(title = miRNA_id, theme = theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold")))
  ggsave(paste0("5_GMC_analysis/plot_bubble/", miRNA_title, ".pdf"), merge, width = 6, height = 6, units = "in", dpi = 300)
}

# Fig4F ------------------------------
every_nth <- function(x, nth, empty = TRUE, inverse = FALSE) 
{
  if (!inverse) {
    if(empty) {
      x[1:nth == 1] <- ""
      x
    } else {
      x[1:nth != 1]
    }
  } else {
    if(empty) {
      x[1:nth != 1] <- ""
      x
    } else {
      x[1:nth == 1]
    }
  }
}

plot_seqlogo <- function(mtfile,title,switch,xmax){
  mt <- read.table(mtfile)
  mx <- consensusMatrix(mt$V1)
  if(ncol(mx)!=30){
    mx.tmp <- matrix(c(0,0,0,0),ncol = 30-ncol(mx),nrow = nrow(mx))
    mx <- cbind(mx,mx.tmp)
  }
  sum <- colSums(mx)[1]
  custom_breaks = seq(0, 30, 1)
  if(switch == "1"){
    ggseqlogo(mx, method = "custom", seq_type = "rna", scales = "none") +
      geom_rect(aes(xmin = 0, xmax = xmax+0.5, ymin = -Inf, ymax = Inf), fill = "#555555", alpha = 0.1) +
      theme_bw()+
      theme(panel.grid.minor = element_line(colour = NA),
            panel.grid.major = element_line(colour = NA),
            panel.background = element_rect(fill = NA),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title = element_text(size = 9, face = "bold"),
            axis.text = element_text(size = 9, face = "bold", colour = "black"),
            plot.title = element_text(size = 9, face = "bold"),
            plot.background = element_rect(colour = NA),
            legend.position = "none") +
      scale_x_continuous(breaks = seq(1, 30, 1), labels = c(1,every_nth(custom_breaks, 5, inverse = TRUE)[-c(1,2)]), expand = c(0, 0))+
      ggtitle("Genome Sequence")
  }else{
    ggseqlogo(mx,method="custom",seq_type="rna",scales = "none")+
      geom_rect(aes(xmin = 0, xmax = xmax+0.5, ymin = -Inf, ymax = Inf), fill = "#555555", alpha = 0.1) +
      theme_bw()+
      theme(panel.grid.minor = element_line(colour = NA),
            panel.grid.major = element_line(colour = NA),
            panel.background = element_rect(fill = NA),
            axis.title = element_text(size = 9, face = "bold"),
            axis.text = element_text(size = 9, face = "bold", colour = "black"),
            plot.title = element_text(size = 9, face = "bold"),
            plot.background = element_rect(colour = NA),
            legend.position = "none") +
      scale_x_continuous(breaks = seq(1, 30, 1), labels = c(1,every_nth(custom_breaks, 5, inverse = TRUE)[-c(1,2)]), expand = c(0, 0))+
      ggtitle(title)
  }
}

for(x in 1:length(miRNA_names)){
   pdf.name <- gsub("txt","pdf",miRNA_names[x])
   miRNA <- gsub(".txt","",miRNA_names[x])
   miRNA <- gsub("_star", "*", miRNA)
   ma <- miRNA_ma[which(miRNA_ma$single == miRNA),2]
   len <- miRNA_length[which(miRNA_length$V1 == miRNA),4]
   print(paste0("Processing(",x,"/",length(miRNA_names),"):",miRNA))
   a <- plot_seqlogo(mtfile = paste0("/bios-store1/chenyc/DataBase_TRM-sRNA-seq/QC/tailing_and_trimming/5_GMC_analysis/doc_seqlogo/LerCK_20201010N_ACAGTGAT_S1_L001_R1_001/", miRNA_names[x]), title = expression(bold("L."~bolditalic("er")~" conventional sRNA-seq")), switch = 2, xmax = len)
   b <- plot_seqlogo(mtfile = paste0("/bios-store1/chenyc/DataBase_TRM-sRNA-seq/QC/tailing_and_trimming/5_GMC_analysis/doc_seqlogo/LerOXAP_20201010N_ACTTGAAT_S4_L002_R1_001/", miRNA_names[x]),title = expression(bold("L."~bolditalic("er")~" TRM-sRNA-seq")),switch = 2, xmax = len)

   c <- plot_seqlogo(mtfile = paste0("/bios-store1/chenyc/DataBase_TRM-sRNA-seq/QC/tailing_and_trimming/5_GMC_analysis/doc_seqlogo/henCK_20201010N_GATCAGAT_S5_L002_R1_001/", miRNA_names[x]), title = expression(bold(~bolditalic("hen1-2")~" conventional sRNA-seq")), switch = 2, xmax = len)
   d <- plot_seqlogo(mtfile = paste0("/bios-store1/chenyc/DataBase_TRM-sRNA-seq/QC/tailing_and_trimming/5_GMC_analysis/doc_seqlogo/henOXAP_20201010N_CTTGTAAT_S8_L002_R1_001/", miRNA_names[x]), title = expression(bold(~bolditalic("hen1-2")~" TRM-sRNA-seq")), switch = 2, xmax = len)
   e <- plot_seqlogo(mtfile = paste0("/bios-store1/chenyc/scripts/Tailing_Trimming/miRNA_template/", miRNA_names[x]), switch = 1 ,title = NULL, xmax = len)
   merge <- a/b/c/d/e + plot_annotation(title = paste0(miRNA, "   ",ma) , theme = theme(plot.title = element_text(size = 9, hjust = 0.5, face = "bold")))
   ggsave(merge,file = paste0("5_GMC_analysis/plot_seqlogo/",pdf.name),width = 5, height = 8, units = "in")
}

# Fig5AB ------------------------------
method_group <- list(
  lck = c("LCK_20200818N","LerCK_20201010N"),
  lpba = c("LAP_20200818N", "LerAP_20201010N"),
  lox = c("LOX_20200818N","LerOX_20201010N"),
  lpbox = c("LOXAP_20200818N", "LerOXAP_20201010N"),
  hck = c("hCK_20200818N","henCK_20201010N"),
  hpba = c("hAP_20200818N", "henAP_20201010N"),
  hox = c("hOX_20200818N","henOX_20201010N"),
  hpbox = c("hOXAP_20200818N", "henOXAP_20201010N"))

for (sample in names(method_group)) {
  
  group <- method_group[[sample]]
  df1 <- read.table(paste0("Fig5AB/",group[1],"_trimmed_priority.summary"), header = TRUE, sep = '\t',strip.white = TRUE)
  df2 <- read.table(paste0("Fig5AB/",group[2],"_trimmed_priority.summary"), header = TRUE, sep = '\t',strip.white = TRUE)
  df <- full_join(df1, df2, by = c("length","feature"))
  df$mean <- rowMeans(df[, c("Percentage.x", "Percentage.y")])
  df$feature <- gsub("otherRNA","unassigned",df$feature)
  df$feature <- gsub("miRNA_primary_transcript","miRNA",df$feature)
  df$feature <- gsub("phasi_tasiRNA","phasi/tasiRNA",df$feature)
  df$feature <- factor(df$feature, levels = rnatype_level)
  
  cor_res <- cor.test(df$Percentage.x,df$Percentage.y)
  
  if (sample == "hck") {
    xlab = expression(italic("hen1-2")~" Conventional sRNA-seq")
  }else if(sample == "hpba"){
    xlab = expression(italic("hen1-2")~" PBA sRNA-seq")
  }else if(sample == "hox"){
    xlab = expression(italic("hen1-2")~" OX sRNA-seq")
  }else if(sample == "hpbox"){
    xlab = expression(italic("hen1-2")~" PBOX sRNA-seq")
  }else if(sample == "lck") {
    xlab = expression("L."~italic("er")~" Conventional sRNA-seq")
  }else if(sample == "lpba"){
    xlab = expression("L."~italic("er")~" PBA sRNA-seq")
  }else if(sample == "lox"){
    xlab = expression("L."~italic("er")~" OX sRNA-seq")
  }else if(sample == "lpbox"){
    xlab = expression("L."~italic("er")~" PBOX sRNA-seq")
  }

  tmp <- ggplot(df) +
    aes(x = length, fill = feature, weight = mean) +
    geom_bar() +
    scale_fill_manual(values = rnatype_colour) +
    theme_classic()+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          #axis.line.x = element_line(size = 1),#轴线粗细
          axis.text.y = element_text(size = 12),
          #axis.line.y = element_line(size = 1),#轴线粗细
          plot.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))+
    labs(x = xlab,y = "Relative percentage(%)",fill = "RNA type")+
    annotate("text",label=paste0("R = ", round(cor_res$estimate,4), " p<0.001"), x=32, y=45)+
    scale_x_continuous(breaks = seq(10, 35, 1),expand = c(0,0.5))+
    scale_y_continuous(breaks = seq(0,50,5), limits = c(0,50),expand = c(0,0.5))
  assign(sample, tmp)

}

ler <- lck+lpba+lox+lpbox+plot_layout(guides = 'collect')+plot_annotation(tag_levels = "A")
ler
ggsave("Fig5AB/ler_merge.pdf",width = 15, height = 8)

hen <- hck+hpba+hox+hpbox+plot_layout(guides = 'collect')+plot_annotation(tag_levels = "A")
hen
ggsave("Fig5AB/hen_merge.pdf",width = 15, height = 8)




