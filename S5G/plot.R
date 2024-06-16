library(openxlsx)
library(tidyverse)
library(scales)
library(ggbreak)
library(readxl)
library(writexl)

setwd("D:/^^Project/TRM-sRNA-seq/publish results/Data Source/R14")


meta <- read.table("doc/mapping_results_bowtie_cleandata.txt", header = TRUE)

filelist <- list.files(path = "doc/", pattern = "*tail_1_nt_duplicated_for_publish.xls",include.dirs = TRUE,recursive = TRUE)
file <- filelist[1]
for (file in filelist) {
  samplename <- strsplit(file, "/")[[1]][1]
  df <- read.xlsx(paste0("doc/",file),sheet = 2)
  df <- df %>% arrange(desc(SUM))
  
  tmp <- df[1:20, ]
  tmp <- tmp  %>% 
    mutate(tail_1  = (A+C+G+T)/meta[meta$Sample == samplename, 3]*1000000,
           RPM = SUM/meta[meta$Sample == samplename, 3]*1000000,
           ) %>% select(ID, RPM, tail_1, 'miRNA_tail_1/miRNA_total')
  tmp$ID <- reorder(tmp$ID, desc(tmp$RPM))
  
  scale.max <- ceiling(max(tmp$RPM))
  scale.min <- floor(max(tmp$`miRNA_tail_1/miRNA_total`))
  
  #p <- 
  
  
  
  p <- ggplot(tmp, aes(x=ID)) +
    geom_bar( aes(y=RPM), stat="identity", size=.1) + 
    geom_bar( aes(y=tail_1), stat="identity", size=.1,fill = "red") + 
    geom_line( aes(y=rescale(`miRNA_tail_1/miRNA_total`,c(0,scale.max)),group = 1), linewidth=0.7) +
    scale_y_continuous(
      # Features of the first axis
      name = "Reads Count",
      expand = c(0,0),
      # Add a second axis and specify its features
      sec.axis = sec_axis(~rescale(.,c(0,scale.min)), name="miRNA tail 1/miRNA total\nRelative percentage(%)")
    ) +
    scale_y_break(
    c(ceiling(tmp[3,2]),ceiling(tmp[3,2])+100),#截断位置及范围
        space = 0.3,#间距大小
        scales = 1.5,expand = FALSE)+ #上下显示比例，大于1上面比例大，小于1下面比例大
    theme_bw()+
    theme(axis.text = element_text(size = 9,
      hjust = 1, vjust = 0), axis.text.x = element_text(size = 9,
      vjust = 0, angle = 90), axis.text.y = element_text(size = 9))+
    labs(x = NULL)

  
  ggsave(paste0(samplename, "_tail1_percentage.pdf"), p, width = 10, height = 15)
}


for (file in filelist) {
  df <- read.xlsx(paste0("doc/",file),sheet = 2)
  df <- df %>% mutate(RPM = SUM/meta[meta$Sample == samplename, 3]*1000000) %>% 
    select(1,2,3,RPM, everything() )
  write.xlsx(df, paste0("doc/",file), sheet_name = "tail_1_summary" , overwrite = TRUE)
}
