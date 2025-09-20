#!/usr/bin/env Rscript
# R语言基础课程 - 环境设置脚本
# 作者：王运生
# 日期：2025-01-20
# 用途：配置R环境并安装必需的包

# 清理环境
cat("清理R环境...\n")
rm(list = ls())

# 设置基本选项
options(stringsAsFactors = FALSE)  # 避免自动转换为因子
options(digits = 4)                # 数字显示精度
options(warn = 1)                  # 显示警告信息

# 设置随机种子
set.seed(123)

# 检查R版本
cat("检查R版本...\n")
r_version <- R.version.string
cat("当前R版本:", r_version, "\n")

if (R.version$major < 4) {
  warning("建议使用R 4.0或更高版本")
}

# 设置CRAN镜像
cat("设置CRAN镜像...\n")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 定义需要安装的包
required_packages <- c(
  "tidyverse",    # 数据科学包集合
  "ggplot2",      # 数据可视化
  "dplyr",        # 数据处理
  "readr",        # 数据读取
  "readxl",       # Excel文件读取
  "writexl",      # Excel文件写入
  "reshape2",     # 数据重塑
  "corrplot",     # 相关性图
  "gridExtra",    # 图形组合
  "RColorBrewer"  # 颜色方案
)

# 安装和加载CRAN包
cat("检查并安装CRAN包...\n")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("安装包:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
    
    # 验证安装
    if (require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("✓ 成功安装并加载:", pkg, "\n")
    } else {
      cat("✗ 安装失败:", pkg, "\n")
    }
  } else {
    cat("✓ 已安装:", pkg, "\n")
  }
}

# 安装Bioconductor包管理器
cat("检查BiocManager...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("安装BiocManager...\n")
  install.packages("BiocManager")
}

# 定义需要安装的Bioconductor包
bioc_packages <- c(
  "Biostrings",      # 生物序列分析
  "GenomicRanges",   # 基因组区间操作
  "IRanges"          # 区间操作基础包
)

# 安装和加载Bioconductor包
cat("检查并安装Bioconductor包...\n")
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("安装Bioconductor包:", pkg, "\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
    
    # 验证安装
    if (require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("✓ 成功安装并加载:", pkg, "\n")
    } else {
      cat("✗ 安装失败:", pkg, "\n")
    }
  } else {
    cat("✓ 已安装:", pkg, "\n")
  }
}

# 创建示例数据
cat("创建示例数据文件...\n")

# 基因信息数据
gene_info <- data.frame(
  gene_name = c("BRCA1", "BRCA2", "TP53", "EGFR", "MYC"),
  chromosome = c(17, 13, 17, 7, 8),
  start_position = c(43044295, 32315086, 7668402, 55019017, 127735434),
  end_position = c(43125483, 32400266, 7687538, 55207337, 127742951),
  gene_type = c("protein_coding", "protein_coding", "protein_coding", 
                "protein_coding", "protein_coding"),
  mutation_status = c("mutated", "wild_type", "mutated", "wild_type", "mutated"),
  expression = c(5.2, 3.8, 7.1, 4.5, 6.3),
  stringsAsFactors = FALSE
)

# 保存基因信息
write.csv(gene_info, "gene_expression.csv", row.names = FALSE)

# 创建表达矩阵
set.seed(123)
gene_names <- c("BRCA1", "BRCA2", "TP53", "EGFR", "MYC")
sample_names <- c("Control1", "Control2", "Treatment1", "Treatment2")

expression_matrix <- matrix(
  round(abs(rnorm(20, mean = 5, sd = 1.5)), 2),
  nrow = 5, ncol = 4,
  dimnames = list(gene_names, sample_names)
)

# 保存表达矩阵
write.table(expression_matrix, "expression_matrix.txt", 
            sep = "\t", quote = FALSE)

# 创建样本信息
sample_info <- data.frame(
  sample_id = sample_names,
  condition = c("control", "control", "treatment", "treatment"),
  batch = c(1, 1, 2, 2),
  stringsAsFactors = FALSE
)

write.table(sample_info, "sample_info.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# 创建FASTA序列文件
fasta_content <- ">gene1|BRCA1
ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>gene2|TP53
ATGCCCGGGAAATTTAAATTTCCCGGGAAATAGATCGATCGATCGATCGATC
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
>gene3|EGFR
ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"

writeLines(fasta_content, "sequences.fasta")

# 显示会话信息
cat("\n" , "="*50, "\n")
cat("R环境配置完成！\n")
cat("="*50, "\n")

cat("工作目录:", getwd(), "\n")
cat("R版本:", r_version, "\n")

# 显示已加载的包
loaded_packages <- search()[grepl("package:", search())]
cat("已加载的包:\n")
for (pkg in loaded_packages) {
  cat(" -", pkg, "\n")
}

# 显示创建的文件
created_files <- c("gene_expression.csv", "expression_matrix.txt", 
                   "sample_info.txt", "sequences.fasta")
cat("创建的示例数据文件:\n")
for (file in created_files) {
  if (file.exists(file)) {
    cat(" ✓", file, "\n")
  } else {
    cat(" ✗", file, "(创建失败)\n")
  }
}

cat("\n环境设置完成，可以开始R语言学习！\n")

# 简单的功能测试
cat("\n执行简单功能测试...\n")

# 测试基本运算
test_result <- 2 + 3
cat("基本运算测试 (2+3):", test_result, "\n")

# 测试数据框操作
test_df <- head(gene_info, 3)
cat("数据框操作测试:\n")
print(test_df)

# 测试绘图功能
if (require(ggplot2, quietly = TRUE)) {
  cat("ggplot2绘图功能正常\n")
} else {
  cat("ggplot2绘图功能异常\n")
}

cat("\n所有测试完成！\n")