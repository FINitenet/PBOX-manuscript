library(tidyverse)
library(openxlsx)
library(ggpubr)
library(lemon)
library(ggrastr)
library(GGally)
library(scales)
library(patchwork)


meta <- read.table("/bios-store1/chenyc/Project/TRM-sRNA-seq/QC/2_srna_analysis/02.Mapping_10_35/mapping_results_ShortStack.txt", header = TRUE, sep = '\t')

mypath <- "/bios-store1/chenyc/Project/TRM-sRNA-seq/QC/2_srna_analysis/"

multmerge <- function(mypath) {
  filenames <- list.files(path = mypath, pattern = "*counts.summary", full.names = TRUE)
  datalist <- lapply(filenames, function(x) {
    fread(file = x)
  })
  Reduce(function(x, y) {
    merge(x, y, by = c("Geneid","biotype"), all = T)
  }, datalist)
}

count <- multmerge(mypath)

count <- count %>% select(1, 2, starts_with("L")) %>% filter_if(is.numeric,any_vars(.>10))

setDT(count)
for (i in 3:length(count)) {
  col_name <- paste0(names(count)[i], "_rpm")  # 构造新列名
  rpm_values <- count[[names(count)[i]]] / as.numeric(meta[match(names(count)[i], meta$Sample), 3]) * 10^6
  set(count, j = col_name, value = rpm_values)  # 通过引用添加新列
}

groups <- list(
  CK = c("LCK_20200818N_rpm","LerCK_20201010N_rpm"),
  PBA = c("LAP_20200818N_rpm", "LerAP_20201010N_rpm"),
  OX = c("LOX_20200818N_rpm", "LerOX_20201010N_rpm"),
  PBOX = c("LOXAP_20200818N_rpm", "LerOXAP_20201010N_rpm")
)

count <- count %>% dplyr::mutate(CK_avg = rowMeans(select(., groups[["CK"]])),
                                 PBA_avg = rowMeans(select(., groups[["PBA"]])),
                                 OX_avg = rowMeans(select(., groups[["OX"]])),
                                 PBOX_avg = rowMeans(select(., groups[["PBOX"]])))



df_avg <- count %>% 
  select(1,biotype,ends_with("_avg")) %>%
  filter_if(is.numeric, any_vars(. > 1)) %>% 
  mutate(across(all_of(ends_with("_avg")), ~ log10(.+1), .names = "{.col}_log"))

df_avg <- df_avg %>%
  mutate(
    biotype = case_when(
      biotype == "protein_coding" ~ "protein coding",
      biotype == "otherRNA" ~ "unassigned",
      biotype == "miRNA_primary_transcript" ~ "miRNA",
      biotype == "hc-siRNA" ~ "hc-siRNA\ntransposable element",
      biotype == "phasi" ~ "phasi/tasiRNA",
      biotype == "tasi" ~ "phasi/tasiRNA",
      biotype == "transposable_element" ~ "hc-siRNA\ntransposable element",
      TRUE ~ biotype
    )
  )
names(df_avg)[7:10] <- c("Conventional sRNA-seq","PBA-sRNA-seq","OX-sRNA-seq","PBOX-sRNA-seq")

df_avg <- as.data.frame(df_avg)

combined_df <- data.frame(
  Column1 = numeric(),
  Column2 = numeric(),
  Combined1 = character(),
  Combined2 = character(),
  stringsAsFactors = TRUE
)

# 生成两两列的组合并存储在新的数据框中
for (i in 7:(ncol(df_avg))) {
  for (j in 7:(ncol(df_avg))) {
    combined <- paste(colnames(df_avg)[i], "vs", colnames(df_avg)[j], sep = "_")
    tmp <- data.frame(Column1 = df_avg[,i],Column2 = df_avg[,j],Combined1 = colnames(df_avg)[i],Combined2 = colnames(df_avg)[j])
    combined_df <- rbind(combined_df, tmp)
  }
}

rnatype_level1 <- c("unassigned",
                    "snoRNA","snRNA","rRNA","tRNA",
                    "lncRNA","protein coding",
                    "hc-siRNA\ntransposable element","phasi/tasiRNA","miRNA")

rnatype_colour1 <- c(unassigned = "#3694ff",
                     snoRNA = "#33a02c",snRNA = "#fb9a99",rRNA = "#1f78b4",tRNA = "#cab2d6",
                     lncRNA = "#e31a1c",`protein coding` = "#fdbf6f",
                     `hc-siRNA\ntransposable element` = "#e0e0e0",`phasi/tasiRNA` = "#8f7ebb",miRNA = "#ee8f55")

methods <- c("Conventional sRNA-seq","PBA-sRNA-seq","OX-sRNA-seq","PBOX-sRNA-seq")
sirna <- factor(c("hc-siRNA\ntransposable element","phasi/tasiRNA","miRNA"), levels = c("hc-siRNA\ntransposable element","phasi/tasiRNA","miRNA"))
srna <- factor(c("snoRNA","snRNA","rRNA","tRNA"), levels = c("snoRNA","snRNA","rRNA","tRNA"))
rna <- factor(c("unassigned", "lncRNA","protein coding"),levels = c("unassigned", "lncRNA","protein coding"))


sirna_plot <- function(tmp, x1, y1, vs){
  ggplot(tmp, aes(x = x1, y = y1,color = biotype))+
    geom_point_rast(size = 0.2)+
    geom_point_rast(data = subset(tmp, biotype == "phasi/tasiRNA"),size = 0.5)+
    geom_point_rast(data = subset(tmp, biotype == "miRNA"),size = 0.5,)+
    geom_abline(slope = 1,intercept = 0, lty="dashed", color = "grey")+
    scale_color_manual(values = rnatype_colour1)+
    scale_x_continuous(limits = c(0,5), breaks = seq(-1,5,1),expand = c(0,0.5))+
    scale_y_continuous(limits = c(0,5), breaks = seq(-1,5,1),expand = c(0,0.5))+
    stat_cor(color = "black", aes(label = after_stat(r.label)))+
    # labs(y = ylab, x = xlab, color=NULL, title = vs)+
    labs(y = ylab, x = xlab, color=NULL)+
    theme_bw()+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(size = 12,hjust = 0.5),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          panel.grid = element_blank(),
          panel.border = element_rect(fill=NA),
    )+
    guides(color = guide_legend(reverse = T))+
    coord_fixed()
}

rna_plot <- function(tmp, x1, y1, vs){
  ggplot(tmp, aes(x = x1, y = y1,color = biotype))+
    geom_point_rast(size = 0.2)+
    geom_abline(slope = 1,intercept = 0, lty="dashed", color = "grey")+
    scale_color_manual(values = rnatype_colour1)+
    scale_x_continuous(limits = c(0,5), breaks = seq(-1,5,1),expand = c(0,0.5))+
    scale_y_continuous(limits = c(0,5), breaks = seq(-1,5,1),expand = c(0,0.5))+
    stat_cor(color = "black", aes(label = after_stat(r.label)))+
    # labs(y = ylab, x = xlab, color=NULL,title = vs)+
    labs(y = ylab, x = xlab, color=NULL)+
    theme_bw()+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(size = 12,hjust = 0.5),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          panel.grid = element_blank(),
          panel.border = element_rect(fill=NA),
    )+
    coord_fixed()
}


for (type in c("sirna","srna","rna")) {
  tmp_p <- wrap_elements(grid::textGrob(type,hjust = 0.01))
  for (i in 1:(length(methods) - 1)) {
    for (j in (i+1):length(methods)) {
      l1 <- methods[i]
      l2 <- methods[j]
      
      xlab = bquote("log"[10] ~ "RPM("*.(l1)*")")
      ylab = bquote("log"[10] ~ "RPM("*.(l2)*")")
      
      vs <- paste(l1,"vs",l2)
      
      tmp <- df_avg %>% filter(biotype %in% get(type)) %>% select(biotype, x1 = all_of(l1), y1 = all_of(l2))
      
      if (type == "sirna") {
        tmp$biotype <- factor(tmp$biotype, levels = rnatype_level1)
        p1 <- sirna_plot(tmp, x1, y1, vs)
        tmp_p <- tmp_p + p1
      }else{
        p1 <- rna_plot(tmp, x1, y1 , vs)
        tmp_p <- tmp_p + p1
      }
    }
  }
  assign(paste0(type, "_p"), tmp_p + plot_layout(ncol = 7))
  ggsave(paste0("ler_methods_",type,"_0509.pdf"), tmp_p + plot_layout(guides = "collect", ncol = 7), width = 22, height = 6)
}

part <- sirna_p/srna_p/rna_p + plot_layout(guides = "collect")
ggsave("ler_methods_part_0509.pdf", part, width = 25, height = 15)
