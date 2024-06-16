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
  
          axis.text.y = element_text(size = 12),
  
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

        axis.text.y = element_text(size = 12),

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

        axis.text.y = element_text(size = 12),

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

        axis.text.y = element_text(size = 12),

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

        axis.text.y = element_text(size = 12),

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

        axis.text.y = element_text(size = 12),

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