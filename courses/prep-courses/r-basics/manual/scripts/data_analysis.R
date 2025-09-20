#!/usr/bin/env Rscript
# R语言基础课程 - 数据分析脚本
# 作者：王运生
# 日期：2025-01-20
# 用途：数据处理、统计分析和假设检验

# 加载必需的包
library(dplyr)
library(ggplot2)

# 清理环境
rm(list = ls())

cat("R语言数据分析练习\n")
cat("==================\n\n")

# 读取数据
if (file.exists("gene_expression.csv")) {
  gene_info <- read.csv("gene_expression.csv", stringsAsFactors = FALSE)
} else {
  # 如果文件不存在，创建示例数据
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
}

# 计算基因长度
gene_info$gene_length <- gene_info$end_position - gene_info$start_position + 1

cat("数据概览:\n")
print(head(gene_info))
cat("\n")

# 练习1: 使用dplyr进行数据处理
cat("练习1: dplyr数据处理\n")
cat("-------------------\n")

# 基本dplyr操作
processed_data <- gene_info %>%
  # 添加log2表达水平
  mutate(log2_expression = log2(expression + 1)) %>%
  # 添加基因长度分类
  mutate(length_category = case_when(
    gene_length < 50000 ~ "短基因",
    gene_length < 100000 ~ "中等长度",
    TRUE ~ "长基因"
  )) %>%
  # 添加表达水平分类
  mutate(expression_category = case_when(
    expression > 6.0 ~ "高表达",
    expression > 4.0 ~ "中等表达",
    TRUE ~ "低表达"
  )) %>%
  # 添加染色体分组
  mutate(chromosome_group = case_when(
    chromosome <= 10 ~ "低号染色体",
    chromosome <= 20 ~ "中号染色体",
    TRUE ~ "高号染色体"
  )) %>%
  # 按表达水平排序
  arrange(desc(expression))

cat("处理后的数据:\n")
print(processed_data)
cat("\n")

# 数据筛选示例
high_expression_genes <- gene_info %>%
  filter(expression > 5.0) %>%
  select(gene_name, expression, mutation_status)

cat("高表达基因 (>5.0):\n")
print(high_expression_genes)
cat("\n")

# 分组汇总
summary_by_mutation <- gene_info %>%
  group_by(mutation_status) %>%
  summarise(
    基因数量 = n(),
    平均表达水平 = round(mean(expression), 2),
    表达水平标准差 = round(sd(expression), 2),
    中位表达水平 = round(median(expression), 2),
    最小表达水平 = min(expression),
    最大表达水平 = max(expression),
    平均基因长度 = round(mean(gene_length), 0),
    .groups = 'drop'
  )

cat("按突变状态分组的汇总统计:\n")
print(summary_by_mutation)
cat("\n")

# 复杂的数据处理流水线
complex_analysis <- gene_info %>%
  # 添加多个计算列
  mutate(
    log2_expression = log2(expression + 1),
    expression_rank = rank(-expression),
    is_high_expression = expression > median(expression),
    gene_length_kb = round(gene_length / 1000, 1)
  ) %>%
  # 按染色体分组计算组内统计
  group_by(chromosome) %>%
  mutate(
    chr_mean_expression = mean(expression),
    expression_vs_chr_mean = round(expression - chr_mean_expression, 2),
    chr_gene_count = n()
  ) %>%
  ungroup() %>%
  # 最终筛选和排序
  filter(gene_length > 50000) %>%
  arrange(chromosome, start_position) %>%
  select(gene_name, chromosome, expression, log2_expression, 
         gene_length_kb, expression_rank, chr_mean_expression, 
         expression_vs_chr_mean)

cat("复杂数据处理结果:\n")
print(complex_analysis)
cat("\n")

# 练习2: 描述性统计分析
cat("练习2: 描述性统计分析\n")
cat("---------------------\n")

# 基本描述性统计
expression_values <- gene_info$expression

cat("基因表达水平描述性统计:\n")
cat("样本数:", length(expression_values), "\n")
cat("平均值:", round(mean(expression_values), 3), "\n")
cat("中位数:", round(median(expression_values), 3), "\n")
cat("标准差:", round(sd(expression_values), 3), "\n")
cat("方差:", round(var(expression_values), 3), "\n")
cat("最小值:", min(expression_values), "\n")
cat("最大值:", max(expression_values), "\n")
cat("范围:", round(diff(range(expression_values)), 3), "\n")
cat("变异系数:", round(sd(expression_values) / mean(expression_values), 3), "\n")

# 分位数分析
quantiles <- quantile(expression_values, probs = c(0.25, 0.5, 0.75))
cat("四分位数:\n")
cat("  Q1 (25%):", quantiles[1], "\n")
cat("  Q2 (50%):", quantiles[2], "\n")
cat("  Q3 (75%):", quantiles[3], "\n")
cat("  IQR:", quantiles[3] - quantiles[1], "\n\n")

# 详细统计摘要
cat("详细统计摘要:\n")
print(summary(expression_values))
cat("\n")

# 按分组的详细统计
cat("按突变状态分组的详细统计:\n")
by_mutation_detailed <- gene_info %>%
  group_by(mutation_status) %>%
  summarise(
    样本数 = n(),
    平均值 = round(mean(expression), 3),
    中位数 = round(median(expression), 3),
    标准差 = round(sd(expression), 3),
    最小值 = min(expression),
    最大值 = max(expression),
    Q1 = round(quantile(expression, 0.25), 3),
    Q3 = round(quantile(expression, 0.75), 3),
    变异系数 = round(sd(expression) / mean(expression), 3),
    .groups = 'drop'
  )
print(by_mutation_detailed)
cat("\n")

# 数值变量相关性分析
numeric_data <- gene_info %>%
  select(chromosome, start_position, end_position, gene_length, expression)

correlation_matrix <- cor(numeric_data)
cat("数值变量相关性矩阵:\n")
print(round(correlation_matrix, 3))
cat("\n")

# 练习3: 假设检验
cat("练习3: 假设检验\n")
cat("---------------\n")

# 准备测试数据
mutated_expression <- gene_info$expression[gene_info$mutation_status == "mutated"]
wildtype_expression <- gene_info$expression[gene_info$mutation_status == "wild_type"]

cat("突变基因表达水平:", mutated_expression, "\n")
cat("野生型基因表达水平:", wildtype_expression, "\n\n")

# 单样本t检验
cat("1. 单样本t检验\n")
cat("检验突变基因的平均表达水平是否显著不同于5.0\n")
single_t_test <- t.test(mutated_expression, mu = 5.0)
cat("H0: μ = 5.0 vs H1: μ ≠ 5.0\n")
cat("t统计量:", round(single_t_test$statistic, 4), "\n")
cat("自由度:", single_t_test$parameter, "\n")
cat("p值:", round(single_t_test$p.value, 4), "\n")
cat("95%置信区间: [", round(single_t_test$conf.int[1], 3), 
    ", ", round(single_t_test$conf.int[2], 3), "]\n")
cat("结论:", ifelse(single_t_test$p.value < 0.05, 
                "拒绝原假设，平均表达水平显著不等于5.0", 
                "不能拒绝原假设，平均表达水平不显著不同于5.0"), "\n\n")

# 双样本t检验
cat("2. 双样本t检验\n")
cat("检验突变基因和野生型基因的表达水平是否有显著差异\n")
two_sample_t_test <- t.test(mutated_expression, wildtype_expression)
cat("H0: μ1 = μ2 vs H1: μ1 ≠ μ2\n")
cat("t统计量:", round(two_sample_t_test$statistic, 4), "\n")
cat("自由度:", round(two_sample_t_test$parameter, 2), "\n")
cat("p值:", round(two_sample_t_test$p.value, 4), "\n")
cat("95%置信区间: [", round(two_sample_t_test$conf.int[1], 3), 
    ", ", round(two_sample_t_test$conf.int[2], 3), "]\n")
cat("突变组平均值:", round(mean(mutated_expression), 3), "\n")
cat("野生型组平均值:", round(mean(wildtype_expression), 3), "\n")
cat("结论:", ifelse(two_sample_t_test$p.value < 0.05, 
                "拒绝原假设，两组表达水平有显著差异", 
                "不能拒绝原假设，两组表达水平无显著差异"), "\n\n")

# Welch t检验（不等方差）
cat("3. Welch t检验 (不假设等方差)\n")
welch_t_test <- t.test(mutated_expression, wildtype_expression, var.equal = FALSE)
cat("t统计量:", round(welch_t_test$statistic, 4), "\n")
cat("自由度:", round(welch_t_test$parameter, 2), "\n")
cat("p值:", round(welch_t_test$p.value, 4), "\n\n")

# 方差齐性检验
cat("4. 方差齐性检验 (F检验)\n")
var_test <- var.test(mutated_expression, wildtype_expression)
cat("F统计量:", round(var_test$statistic, 4), "\n")
cat("p值:", round(var_test$p.value, 4), "\n")
cat("结论:", ifelse(var_test$p.value < 0.05, 
                "拒绝原假设，两组方差不相等", 
                "不能拒绝原假设，两组方差相等"), "\n\n")

# Wilcoxon秩和检验（非参数检验）
cat("5. Wilcoxon秩和检验 (非参数)\n")
wilcox_test <- wilcox.test(mutated_expression, wildtype_expression)
cat("W统计量:", wilcox_test$statistic, "\n")
cat("p值:", round(wilcox_test$p.value, 4), "\n")
cat("结论:", ifelse(wilcox_test$p.value < 0.05, 
                "拒绝原假设，两组分布有显著差异", 
                "不能拒绝原假设，两组分布无显著差异"), "\n\n")

# 相关性检验
cat("6. 相关性检验\n")
cat("检验基因长度与表达水平的相关性\n")
cor_test <- cor.test(gene_info$gene_length, gene_info$expression)
cat("Pearson相关系数:", round(cor_test$estimate, 4), "\n")
cat("t统计量:", round(cor_test$statistic, 4), "\n")
cat("自由度:", cor_test$parameter, "\n")
cat("p值:", round(cor_test$p.value, 4), "\n")
cat("95%置信区间: [", round(cor_test$conf.int[1], 3), 
    ", ", round(cor_test$conf.int[2], 3), "]\n")
cat("结论:", ifelse(cor_test$p.value < 0.05, 
                "存在显著相关性", 
                "不存在显著相关性"), "\n\n")

# Spearman秩相关检验
spearman_test <- cor.test(gene_info$gene_length, gene_info$expression, method = "spearman")
cat("Spearman秩相关系数:", round(spearman_test$estimate, 4), "\n")
cat("p值:", round(spearman_test$p.value, 4), "\n\n")

# 正态性检验
cat("7. 正态性检验\n")
cat("检验表达水平是否服从正态分布\n")
shapiro_test <- shapiro.test(expression_values)
cat("Shapiro-Wilk检验:\n")
cat("W统计量:", round(shapiro_test$statistic, 4), "\n")
cat("p值:", round(shapiro_test$p.value, 4), "\n")
cat("结论:", ifelse(shapiro_test$p.value < 0.05, 
                "拒绝原假设，数据不服从正态分布", 
                "不能拒绝原假设，数据可能服从正态分布"), "\n\n")

# Kolmogorov-Smirnov检验
ks_test <- ks.test(expression_values, "pnorm", mean(expression_values), sd(expression_values))
cat("Kolmogorov-Smirnov检验:\n")
cat("D统计量:", round(ks_test$statistic, 4), "\n")
cat("p值:", round(ks_test$p.value, 4), "\n\n")

# 练习4: 回归分析
cat("练习4: 回归分析\n")
cat("---------------\n")

# 简单线性回归
cat("1. 简单线性回归: 表达水平 ~ 基因长度\n")
simple_model <- lm(expression ~ gene_length, data = gene_info)
cat("回归方程: expression =", round(coef(simple_model)[1], 6), 
    "+", round(coef(simple_model)[2], 8), "* gene_length\n")

# 模型摘要
model_summary <- summary(simple_model)
cat("R²:", round(model_summary$r.squared, 4), "\n")
cat("调整R²:", round(model_summary$adj.r.squared, 4), "\n")
cat("F统计量:", round(model_summary$fstatistic[1], 4), "\n")
cat("p值:", round(pf(model_summary$fstatistic[1], 
                   model_summary$fstatistic[2], 
                   model_summary$fstatistic[3], 
                   lower.tail = FALSE), 4), "\n")

# 系数检验
coef_table <- model_summary$coefficients
cat("系数检验:\n")
print(round(coef_table, 6))
cat("\n")

# 残差分析
residuals <- residuals(simple_model)
cat("残差统计:\n")
cat("残差平均值:", round(mean(residuals), 6), "\n")
cat("残差标准差:", round(sd(residuals), 4), "\n")
cat("残差范围:", round(range(residuals), 4), "\n\n")

# 预测
new_lengths <- data.frame(gene_length = c(50000, 75000, 100000))
predictions <- predict(simple_model, new_lengths, interval = "confidence")
cat("预测结果:\n")
result_df <- cbind(new_lengths, round(predictions, 3))
print(result_df)
cat("\n")

# 多元回归（如果有更多变量）
if (ncol(numeric_data) > 2) {
  cat("2. 多元线性回归\n")
  # 创建一些额外的预测变量
  gene_info$log_length <- log10(gene_info$gene_length)
  gene_info$chr_numeric <- as.numeric(gene_info$chromosome)
  
  multi_model <- lm(expression ~ gene_length + chr_numeric + log_length, data = gene_info)
  multi_summary <- summary(multi_model)
  
  cat("多元回归R²:", round(multi_summary$r.squared, 4), "\n")
  cat("调整R²:", round(multi_summary$adj.r.squared, 4), "\n")
  cat("F统计量:", round(multi_summary$fstatistic[1], 4), "\n")
  
  cat("系数:\n")
  print(round(multi_summary$coefficients, 6))
  cat("\n")
}

# 练习5: 方差分析
cat("练习5: 方差分析\n")
cat("---------------\n")

# 单因素方差分析
cat("单因素方差分析: 表达水平 ~ 突变状态\n")
anova_model <- aov(expression ~ mutation_status, data = gene_info)
anova_summary <- summary(anova_model)
print(anova_summary)

# 提取F统计量和p值
f_stat <- anova_summary[[1]][["F value"]][1]
p_value <- anova_summary[[1]][["Pr(>F)"]][1]
cat("F统计量:", round(f_stat, 4), "\n")
cat("p值:", round(p_value, 4), "\n")
cat("结论:", ifelse(p_value < 0.05, 
                "拒绝原假设，不同突变状态的表达水平有显著差异", 
                "不能拒绝原假设，不同突变状态的表达水平无显著差异"), "\n\n")

# 效应大小计算
ss_total <- sum(anova_summary[[1]][["Sum Sq"]])
ss_between <- anova_summary[[1]][["Sum Sq"]][1]
eta_squared <- ss_between / ss_total
cat("效应大小 (η²):", round(eta_squared, 4), "\n\n")

cat("数据分析练习完成！\n")
cat("==================\n")