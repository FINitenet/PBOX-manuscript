setwd("D:/^^Project/TRM-sRNA-seq/publish results/Data Source/Supplemental table")
library(openxlsx)
library(readxl)
library(tidyverse)

headerstyle <- createStyle(textDecoration = "bold", halign = "center")

# S3 -------
excel_file <- "Table S3 Size distribution and genomic classification of 14-35 nt mapped reads of defrosted sample.xlsx"
xlsx <- readxl::excel_sheets(excel_file)

# 创建一个Excel工作簿
wb <- createWorkbook()

# 循环处理每个工作表
for (sheet_name in xlsx) {
  # 读取当前工作表的数据
  df <- readxl::read_excel(excel_file, sheet = sheet_name)
  addWorksheet(wb, sheet_name)
  # 在这里可以对df进行进一步的处理
  df1 <- pivot_wider(df[-4], names_from = "RNA categories", values_from = Raw_read_count)
  df1[is.na(df1)] <- 0
  
  writeData(wb, sheet_name, "Raw_read_count", startRow = 1, startCol = 2)
  mergeCells(wb, sheet_name, cols = 2:12, rows = 1)
  addStyle(wb, sheet_name,headerstyle, rows = 1, cols = 2)
  writeData(wb, sheet_name, df1, startRow = 2, startCol = 1)
  
  df2 <- pivot_wider(df[-3], names_from = "RNA categories", values_from = Percentage)
  df2[is.na(df2)] <- 0
  writeData(wb, sheet_name, "Percentage", startRow = nrow(df1) + 4, startCol = 2)
  mergeCells(wb, sheet_name, cols = 2:12, rows = nrow(df1) + 4)
  addStyle(wb, sheet_name,headerstyle, rows = nrow(df1) + 4, cols = 2)
  writeData(wb, sheet_name, df2, startRow = nrow(df1) + 5, startCol = 1)
}

# 保存Excel文件
saveWorkbook(wb, paste0("c_",excel_file), overwrite = TRUE)

# S5 -------
excel_file <- "Table S5 Size distribution and genomics classification of 10-35 nt mapped reads from AGO1-RIP-sRNA-seq.xlsx"
xlsx <- readxl::excel_sheets(excel_file)

# 创建一个Excel工作簿
wb <- createWorkbook()

# 循环处理每个工作表
for (sheet_name in xlsx) {
  # 读取当前工作表的数据
  df <- readxl::read_excel(excel_file, sheet = sheet_name)
  addWorksheet(wb, sheet_name)
  # 在这里可以对df进行进一步的处理
  df1 <- pivot_wider(df[-4], names_from = "RNA categories", values_from = Raw_read_count)
  df1[is.na(df1)] <- 0
  
  writeData(wb, sheet_name, "Raw_read_count", startRow = 1, startCol = 2)
  mergeCells(wb, sheet_name, cols = 2:12, rows = 1)
  addStyle(wb, sheet_name,headerstyle, rows = 1, cols = 2)
  writeData(wb, sheet_name, df1, startRow = 2, startCol = 1)
  
  df2 <- pivot_wider(df[-3], names_from = "RNA categories", values_from = Percentage)
  df2[is.na(df2)] <- 0
  writeData(wb, sheet_name, "Percentage", startRow = nrow(df1) + 4, startCol = 2)
  mergeCells(wb, sheet_name, cols = 2:12, rows = nrow(df1) + 4)
  addStyle(wb, sheet_name,headerstyle, rows = nrow(df1) + 4, cols = 2)
  writeData(wb, sheet_name, df2, startRow = nrow(df1) + 5, startCol = 1)
}

# 保存Excel文件
saveWorkbook(wb, paste0("c_",excel_file), overwrite = TRUE)

# S9 -------
excel_file <- "Table S9 Size distribution and genomic classification of 10-35 nt mapped reads of L.er, hen1-2, hen1-1 and hen1-5.xlsx"
xlsx <- readxl::excel_sheets(excel_file)

# 创建一个Excel工作簿
wb <- createWorkbook()

# 循环处理每个工作表
for (sheet_name in xlsx) {
  # 读取当前工作表的数据
  df <- readxl::read_excel(excel_file, sheet = sheet_name)
  addWorksheet(wb, sheet_name)
  # 在这里可以对df进行进一步的处理
  df1 <- pivot_wider(df[-4], names_from = "RNA categories", values_from = Raw_read_count)
  df1[is.na(df1)] <- 0
  
  writeData(wb, sheet_name, "Raw_read_count", startRow = 1, startCol = 2)
  mergeCells(wb, sheet_name, cols = 2:12, rows = 1)
  addStyle(wb, sheet_name,headerstyle, rows = 1, cols = 2)
  writeData(wb, sheet_name, df1, startRow = 2, startCol = 1)
  
  df2 <- pivot_wider(df[-3], names_from = "RNA categories", values_from = Percentage)
  df2[is.na(df2)] <- 0
  writeData(wb, sheet_name, "Percentage", startRow = nrow(df1) + 4, startCol = 2)
  mergeCells(wb, sheet_name, cols = 2:12, rows = nrow(df1) + 4)
  addStyle(wb, sheet_name,headerstyle, rows = nrow(df1) + 4, cols = 2)
  writeData(wb, sheet_name, df2, startRow = nrow(df1) + 5, startCol = 1)
}

# 保存Excel文件
saveWorkbook(wb, paste0("c_",excel_file), overwrite = TRUE)
