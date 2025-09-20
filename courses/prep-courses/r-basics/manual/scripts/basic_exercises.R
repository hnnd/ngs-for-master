#!/usr/bin/env Rscript
# R语言基础课程 - 基础练习脚本
# 作者：王运生
# 日期：2025-01-20
# 用途：R语言基础语法和数据类型练习

# 清理环境
rm(list = ls())

cat("R语言基础练习开始\n")
cat("==================\n\n")

# 练习1: 基本数据类型
cat("练习1: 基本数据类型\n")
cat("------------------\n")

# 创建不同类型的变量
gene_name <- "BRCA1"
chromosome <- 17
gc_content <- 0.42
is_tumor_suppressor <- TRUE
gene_length <- 81188L  # 整数类型

# 显示变量信息
cat("基因名:", gene_name, "- 类型:", class(gene_name), "\n")
cat("染色体:", chromosome, "- 类型:", class(chromosome), "\n")
cat("GC含量:", gc_content, "- 类型:", class(gc_content), "\n")
cat("是否抑癌基因:", is_tumor_suppressor, "- 类型:", class(is_tumor_suppressor), "\n")
cat("基因长度:", gene_length, "- 类型:", class(gene_length), "\n\n")

# 类型检查函数
check_type <- function(var, var_name) {
  cat(var_name, ":\n")
  cat("  值:", var, "\n")
  cat("  类型:", class(var), "\n")
  cat("  模式:", mode(var), "\n")
  cat("  是否数值:", is.numeric(var), "\n")
  cat("  是否字符:", is.character(var), "\n")
  cat("  是否逻辑:", is.logical(var), "\n\n")
}

check_type(gene_name, "基因名")
check_type(gc_content, "GC含量")

# 练习2: 类型转换
cat("练习2: 类型转换\n")
cat("---------------\n")

# 字符串转数值
expression_str <- "5.2"
expression_num <- as.numeric(expression_str)
cat("字符串 '", expression_str, "' 转为数值:", expression_num, "\n")

# 数值转字符串
chromosome_str <- as.character(chromosome)
cat("数值", chromosome, "转为字符串: '", chromosome_str, "'\n")

# 数值转逻辑值
logical_from_num <- as.logical(c(0, 1, 2, -1))
cat("数值转逻辑值:", logical_from_num, "\n")

# 逻辑值转数值
num_from_logical <- as.numeric(c(TRUE, FALSE, TRUE))
cat("逻辑值转数值:", num_from_logical, "\n\n")

# 练习3: 向量操作
cat("练习3: 向量操作\n")
cat("---------------\n")

# 创建向量
genes <- c("BRCA1", "BRCA2", "TP53", "EGFR", "MYC")
expression_levels <- c(5.2, 3.8, 7.1, 4.5, 6.3)
chromosomes <- c(17, 13, 17, 7, 8)

cat("基因向量:", genes, "\n")
cat("表达水平:", expression_levels, "\n")
cat("染色体:", chromosomes, "\n\n")

# 向量基本信息
cat("向量基本信息:\n")
cat("基因数量:", length(genes), "\n")
cat("表达水平范围:", range(expression_levels), "\n")
cat("染色体唯一值:", unique(chromosomes), "\n\n")

# 向量索引
cat("向量索引:\n")
cat("第一个基因:", genes[1], "\n")
cat("最后一个基因:", genes[length(genes)], "\n")
cat("前三个基因:", genes[1:3], "\n")
cat("第2和第4个基因:", genes[c(2, 4)], "\n\n")

# 逻辑索引
high_expression <- expression_levels > 5.0
cat("高表达基因逻辑向量:", high_expression, "\n")
cat("高表达基因名称:", genes[high_expression], "\n")
cat("高表达基因表达水平:", expression_levels[high_expression], "\n\n")

# 命名向量
named_expression <- expression_levels
names(named_expression) <- genes
cat("命名向量:\n")
print(named_expression)
cat("BRCA1的表达水平:", named_expression["BRCA1"], "\n\n")

# 向量运算
cat("向量运算:\n")
log2_expression <- log2(expression_levels)
cat("Log2转换:", round(log2_expression, 2), "\n")

normalized_expression <- (expression_levels - mean(expression_levels)) / sd(expression_levels)
cat("标准化:", round(normalized_expression, 2), "\n")

# 向量统计
cat("统计信息:\n")
cat("平均值:", mean(expression_levels), "\n")
cat("中位数:", median(expression_levels), "\n")
cat("标准差:", sd(expression_levels), "\n")
cat("最小值:", min(expression_levels), "\n")
cat("最大值:", max(expression_levels), "\n\n")

# 练习4: 字符串处理
cat("练习4: 字符串处理\n")
cat("-----------------\n")

# DNA序列处理
dna_sequence <- "ATCGATCGATCGAAATCG"
cat("原始DNA序列:", dna_sequence, "\n")

# 基本字符串操作
cat("序列长度:", nchar(dna_sequence), "\n")
cat("转为大写:", toupper(dna_sequence), "\n")
cat("转为小写:", tolower(dna_sequence), "\n")

# 字符串分割
nucleotides <- strsplit(dna_sequence, "")[[1]]
cat("前10个核苷酸:", nucleotides[1:10], "\n")

# 模式匹配和计数
at_count <- lengths(regmatches(dna_sequence, gregexpr("AT", dna_sequence)))
cg_count <- lengths(regmatches(dna_sequence, gregexpr("CG", dna_sequence)))
cat("AT二核苷酸数量:", at_count, "\n")
cat("CG二核苷酸数量:", cg_count, "\n")

# 字符串替换
rna_sequence <- gsub("T", "U", dna_sequence)
cat("转录为RNA:", rna_sequence, "\n\n")

# 序列分析函数
analyze_dna_sequence <- function(seq) {
  seq <- toupper(seq)
  length <- nchar(seq)
  
  # 计算各碱基数量
  a_count <- lengths(regmatches(seq, gregexpr("A", seq)))
  t_count <- lengths(regmatches(seq, gregexpr("T", seq)))
  c_count <- lengths(regmatches(seq, gregexpr("C", seq)))
  g_count <- lengths(regmatches(seq, gregexpr("G", seq)))
  
  # 计算GC含量
  gc_content <- (g_count + c_count) / length
  
  # 返回结果
  return(list(
    sequence = seq,
    length = length,
    A = a_count,
    T = t_count,
    C = c_count,
    G = g_count,
    GC_content = round(gc_content, 3),
    AT_content = round((a_count + t_count) / length, 3)
  ))
}

# 测试序列分析函数
cat("序列分析结果:\n")
seq_analysis <- analyze_dna_sequence(dna_sequence)
for (name in names(seq_analysis)) {
  cat(name, ":", seq_analysis[[name]], "\n")
}
cat("\n")

# 练习5: 条件语句和循环
cat("练习5: 条件语句和循环\n")
cat("---------------------\n")

# 表达水平分类函数
classify_expression <- function(level) {
  if (level > 6.0) {
    return("高表达")
  } else if (level > 3.0) {
    return("中等表达")
  } else if (level > 1.0) {
    return("低表达")
  } else {
    return("几乎不表达")
  }
}

# 测试分类函数
test_levels <- c(7.5, 4.2, 2.1, 0.5, 6.8)
cat("表达水平分类:\n")
for (i in 1:length(test_levels)) {
  level <- test_levels[i]
  category <- classify_expression(level)
  cat("表达水平", level, ":", category, "\n")
}
cat("\n")

# 批量处理序列
sequences <- c(
  "ATCGATCG",
  "GCTAGCTA", 
  "TTAACCGG",
  "CGATCGAT",
  "AAATTTCCC"
)

cat("批量序列分析:\n")
for (i in 1:length(sequences)) {
  seq <- sequences[i]
  analysis <- analyze_dna_sequence(seq)
  cat("序列", i, ":", seq, 
      "- 长度:", analysis$length, 
      "- GC含量:", analysis$GC_content, "\n")
}
cat("\n")

# 寻找特定模式
pattern <- "CG"
cat("寻找模式 '", pattern, "':\n")
for (i in 1:length(sequences)) {
  seq <- sequences[i]
  count <- lengths(regmatches(seq, gregexpr(pattern, seq)))
  if (count > 0) {
    positions <- gregexpr(pattern, seq)[[1]]
    cat("序列", i, ":", seq, "- 找到", count, "个", pattern, 
        "在位置:", positions, "\n")
  } else {
    cat("序列", i, ":", seq, "- 未找到", pattern, "\n")
  }
}
cat("\n")

# 练习6: 函数定义
cat("练习6: 函数定义\n")
cat("---------------\n")

# 反向互补函数
reverse_complement <- function(sequence) {
  # 定义互补碱基映射
  complement_map <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  
  # 转为大写并分割
  seq_upper <- toupper(sequence)
  bases <- strsplit(seq_upper, "")[[1]]
  
  # 互补
  complement_bases <- complement_map[bases]
  
  # 反向
  reverse_complement_bases <- rev(complement_bases)
  
  # 合并
  result <- paste(reverse_complement_bases, collapse = "")
  
  return(result)
}

# 测试反向互补函数
test_seq <- "ATCGATCG"
rev_comp <- reverse_complement(test_seq)
cat("原序列:", test_seq, "\n")
cat("反向互补:", rev_comp, "\n\n")

# GC含量计算函数
calculate_gc_content <- function(sequence) {
  seq_upper <- toupper(sequence)
  total_length <- nchar(seq_upper)
  
  if (total_length == 0) {
    return(0)
  }
  
  g_count <- lengths(regmatches(seq_upper, gregexpr("G", seq_upper)))
  c_count <- lengths(regmatches(seq_upper, gregexpr("C", seq_upper)))
  
  gc_content <- (g_count + c_count) / total_length
  return(round(gc_content, 3))
}

# 测试GC含量函数
cat("GC含量计算测试:\n")
test_sequences <- c("ATCG", "GCGC", "ATAT", "CGATCG")
for (seq in test_sequences) {
  gc <- calculate_gc_content(seq)
  cat("序列:", seq, "- GC含量:", gc, "\n")
}
cat("\n")

# 综合分析函数
comprehensive_sequence_analysis <- function(sequences) {
  results <- data.frame(
    sequence = character(),
    length = integer(),
    gc_content = numeric(),
    reverse_complement = character(),
    classification = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:length(sequences)) {
    seq <- sequences[i]
    
    # 基本分析
    seq_length <- nchar(seq)
    gc_content <- calculate_gc_content(seq)
    rev_comp <- reverse_complement(seq)
    
    # 分类
    if (gc_content > 0.6) {
      classification <- "高GC"
    } else if (gc_content > 0.4) {
      classification <- "中等GC"
    } else {
      classification <- "低GC"
    }
    
    # 添加到结果
    results <- rbind(results, data.frame(
      sequence = seq,
      length = seq_length,
      gc_content = gc_content,
      reverse_complement = rev_comp,
      classification = classification,
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}

# 测试综合分析函数
cat("综合序列分析:\n")
analysis_results <- comprehensive_sequence_analysis(sequences)
print(analysis_results)

cat("\n基础练习完成！\n")
cat("================\n")