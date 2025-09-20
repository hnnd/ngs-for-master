# R语言基础实践操作手册

## 课程信息
- **课程名称**：高通量测序数据分析 - R语言基础
- **主讲教师**：王运生教授
- **联系邮箱**：wangys@hunau.edu.cn
- **办公室**：16教420室
- **上课地点**：105机房

## 实验目标
通过本次实践课程，学生将能够：
1. 配置R和RStudio开发环境
2. 掌握R语言基础语法和数据结构
3. 熟练使用R进行数据处理和统计分析
4. 创建专业的数据可视化图表
5. 了解并使用生物信息学相关R包
6. 完成基因表达数据的完整分析流程

## 环境要求

### 软件环境
- **R版本**：4.0或更高版本
- **RStudio**：最新版本（推荐）
- **必需R包**：
  ```r
  install.packages(c("tidyverse", "ggplot2", "dplyr", "readr", 
                     "readxl", "writexl", "reshape2", "corrplot"))
  
  # Bioconductor包
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install(c("Biostrings", "GenomicRanges"))
  ```

### 硬件要求
- **内存**：至少4GB RAM
- **存储空间**：至少3GB可用空间
- **网络**：用于下载包和数据

### 数据准备
本次实验将使用以下示例数据文件：
- `gene_expression.csv`：基因表达数据
- `sample_info.txt`：样本信息
- `sequences.fasta`：DNA序列数据

## 操作步骤

### 步骤1：R环境设置和验证

#### 1.1 检查R安装
```r
# 检查R版本
R.version.string
version

# 检查工作目录
getwd()

# 设置工作目录（如果需要）
setwd("/path/to/your/working/directory")

# 查看已安装的包
installed.packages()[1:10, c("Package", "Version")]
```

**预期结果**：显示R 4.0+版本信息和当前工作目录

#### 1.2 安装和加载必需包
```r
# 检查并安装基础包
required_packages <- c("tidyverse", "ggplot2", "dplyr", "readr", 
                       "readxl", "writexl", "reshape2", "corrplot")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# 安装Bioconductor包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

bioc_packages <- c("Biostrings", "GenomicRanges")
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

# 验证安装
sessionInfo()
```

#### 1.3 RStudio环境配置
```r
# 设置RStudio选项
options(stringsAsFactors = FALSE)  # 避免自动转换为因子
options(digits = 4)                # 数字显示精度

# 清理环境
rm(list = ls())  # 清除所有变量

# 设置随机种子（用于可重现的结果）
set.seed(123)

print("R环境配置完成！")
```

**运行脚本**：`scripts/setup.R`

### 步骤2：R基础语法练习

#### 2.1 数据类型和基本操作
```r
# 基本数据类型
gene_name <- "BRCA1"
chromosome <- 17
gc_content <- 0.42
is_tumor_suppressor <- TRUE

# 查看数据类型
cat("基因名:", gene_name, "类型:", class(gene_name), "\n")
cat("染色体:", chromosome, "类型:", class(chromosome), "\n")
cat("GC含量:", gc_content, "类型:", class(gc_content), "\n")
cat("是否抑癌基因:", is_tumor_suppressor, "类型:", class(is_tumor_suppressor), "\n")

# 类型转换练习
numeric_string <- "123.45"
converted_number <- as.numeric(numeric_string)
cat("转换前:", numeric_string, "类型:", class(numeric_string), "\n")
cat("转换后:", converted_number, "类型:", class(converted_number), "\n")

# 逻辑运算
high_expression <- gc_content > 0.4
low_chromosome <- chromosome < 20
result <- high_expression & low_chromosome
cat("高GC含量且低染色体号:", result, "\n")
```

#### 2.2 向量操作练习
```r
# 创建向量
genes <- c("BRCA1", "BRCA2", "TP53", "EGFR", "MYC")
expression_levels <- c(5.2, 3.8, 7.1, 4.5, 6.3)
chromosomes <- c(17, 13, 17, 7, 8)

# 向量信息
cat("基因数量:", length(genes), "\n")
cat("表达水平范围:", range(expression_levels), "\n")
cat("平均表达水平:", mean(expression_levels), "\n")

# 向量索引
cat("第一个基因:", genes[1], "\n")
cat("最后一个基因:", genes[length(genes)], "\n")
cat("前三个基因:", genes[1:3], "\n")

# 逻辑索引
high_expr_genes <- genes[expression_levels > 5.0]
cat("高表达基因:", high_expr_genes, "\n")

# 命名向量
named_expression <- expression_levels
names(named_expression) <- genes
print(named_expression)

# 向量运算
log2_expression <- log2(expression_levels)
normalized_expression <- (expression_levels - mean(expression_levels)) / sd(expression_levels)

cat("Log2转换后的表达水平:\n")
print(log2_expression)
cat("标准化后的表达水平:\n")
print(normalized_expression)
```

#### 2.3 字符串处理练习
```r
# DNA序列处理
dna_sequence <- "ATCGATCGATCGAAATCG"

# 基本字符串操作
cat("序列长度:", nchar(dna_sequence), "\n")
cat("转为大写:", toupper(dna_sequence), "\n")
cat("转为小写:", tolower(dna_sequence), "\n")

# 字符串分割
nucleotides <- strsplit(dna_sequence, "")[[1]]
cat("核苷酸列表:", nucleotides[1:10], "...\n")

# 模式匹配
cg_positions <- gregexpr("CG", dna_sequence)[[1]]
cat("CG二核苷酸位置:", cg_positions, "\n")

# 字符串替换
rna_sequence <- gsub("T", "U", dna_sequence)
cat("转录为RNA:", rna_sequence, "\n")

# 序列分析函数
analyze_sequence <- function(seq) {
  seq <- toupper(seq)
  length <- nchar(seq)
  a_count <- lengths(regmatches(seq, gregexpr("A", seq)))
  t_count <- lengths(regmatches(seq, gregexpr("T", seq)))
  c_count <- lengths(regmatches(seq, gregexpr("C", seq)))
  g_count <- lengths(regmatches(seq, gregexpr("G", seq)))
  gc_content <- (g_count + c_count) / length
  
  return(list(
    length = length,
    A = a_count, T = t_count, C = c_count, G = g_count,
    GC_content = gc_content
  ))
}

# 测试序列分析函数
seq_analysis <- analyze_sequence(dna_sequence)
cat("序列分析结果:\n")
print(seq_analysis)
```

**运行脚本**：`scripts/basic_exercises.R`

### 步骤3：数据结构操作

#### 3.1 矩阵操作练习
```r
# 创建基因表达矩阵
set.seed(123)
gene_names <- c("BRCA1", "BRCA2", "TP53", "EGFR", "MYC")
sample_names <- c("Control1", "Control2", "Treatment1", "Treatment2")

# 生成模拟表达数据
expression_matrix <- matrix(
  round(rnorm(20, mean = 5, sd = 1.5), 2),
  nrow = 5, ncol = 4,
  dimnames = list(gene_names, sample_names)
)

# 确保表达值为正
expression_matrix <- abs(expression_matrix)

cat("基因表达矩阵:\n")
print(expression_matrix)

# 矩阵基本信息
cat("矩阵维度:", dim(expression_matrix), "\n")
cat("行名:", rownames(expression_matrix), "\n")
cat("列名:", colnames(expression_matrix), "\n")

# 矩阵统计
cat("每个基因的平均表达水平:\n")
gene_means <- rowMeans(expression_matrix)
print(gene_means)

cat("每个样本的平均表达水平:\n")
sample_means <- colMeans(expression_matrix)
print(sample_means)

# 矩阵操作
cat("转置矩阵:\n")
print(t(expression_matrix))

# Log2转换
log2_matrix <- log2(expression_matrix + 1)  # 加1避免log(0)
cat("Log2转换后的矩阵:\n")
print(log2_matrix)

# 添加新样本
new_sample <- round(rnorm(5, mean = 4.5, sd = 1.2), 2)
new_sample <- abs(new_sample)
expanded_matrix <- cbind(expression_matrix, Treatment3 = new_sample)
cat("添加新样本后的矩阵:\n")
print(expanded_matrix)
```

#### 3.2 数据框操作练习
```r
# 创建基因信息数据框
gene_info <- data.frame(
  gene_name = c("BRCA1", "BRCA2", "TP53", "EGFR", "MYC"),
  chromosome = c(17, 13, 17, 7, 8),
  start_position = c(43044295, 32315086, 7668402, 55019017, 127735434),
  end_position = c(43125483, 32400266, 7687538, 55207337, 127742951),
  gene_type = c("protein_coding", "protein_coding", "protein_coding", 
                "protein_coding", "protein_coding"),
  mutation_status = c("mutated", "wild_type", "mutated", "wild_type", "mutated"),
  stringsAsFactors = FALSE
)

# 添加表达数据
gene_info$expression <- gene_means[gene_info$gene_name]

# 数据框基本信息
cat("数据框结构:\n")
str(gene_info)

cat("数据框摘要:\n")
summary(gene_info)

cat("前几行数据:\n")
print(head(gene_info))

# 计算基因长度
gene_info$gene_length <- gene_info$end_position - gene_info$start_position + 1

# 添加表达分类
gene_info$expression_category <- ifelse(gene_info$expression > 5.0, "high", "low")

# 数据筛选
high_expr_genes <- subset(gene_info, expression > 5.0)
cat("高表达基因:\n")
print(high_expr_genes)

chr17_genes <- subset(gene_info, chromosome == 17)
cat("17号染色体基因:\n")
print(chr17_genes)

# 数据排序
sorted_by_expression <- gene_info[order(-gene_info$expression), ]
cat("按表达水平降序排列:\n")
print(sorted_by_expression)

# 分组统计
by_mutation <- by(gene_info$expression, gene_info$mutation_status, 
                  function(x) c(mean = mean(x), sd = sd(x), count = length(x)))
cat("按突变状态分组的表达统计:\n")
print(by_mutation)
```

#### 3.3 列表操作练习
```r
# 创建分析结果列表
analysis_results <- list(
  raw_data = expression_matrix,
  gene_info = gene_info,
  statistics = list(
    total_genes = nrow(gene_info),
    total_samples = ncol(expression_matrix),
    mean_expression = mean(gene_info$expression),
    high_expression_count = sum(gene_info$expression > 5.0)
  ),
  parameters = list(
    log_transform = TRUE,
    normalization_method = "z-score",
    p_value_cutoff = 0.05
  )
)

# 访问列表元素
cat("原始数据矩阵:\n")
print(analysis_results$raw_data)

cat("统计信息:\n")
print(analysis_results$statistics)

# 修改列表元素
analysis_results$parameters$fold_change_cutoff <- 2.0
analysis_results$quality_control <- list(
  samples_passed = 4,
  genes_filtered = 0
)

# 列表信息
cat("列表结构:\n")
str(analysis_results, max.level = 2)

# 应用函数到列表
numeric_stats <- analysis_results$statistics[sapply(analysis_results$statistics, is.numeric)]
cat("数值型统计信息:\n")
print(numeric_stats)
```

**运行脚本**：`scripts/data_structures.R`

### 步骤4：数据处理和统计分析

#### 4.1 使用dplyr进行数据处理
```r
library(dplyr)

# 使用管道操作符进行数据处理
processed_data <- gene_info %>%
  # 添加log2表达水平
  mutate(log2_expression = log2(expression + 1)) %>%
  # 添加基因长度分类
  mutate(length_category = case_when(
    gene_length < 50000 ~ "short",
    gene_length < 100000 ~ "medium",
    TRUE ~ "long"
  )) %>%
  # 筛选高表达基因
  filter(expression > 4.0) %>%
  # 按表达水平排序
  arrange(desc(expression)) %>%
  # 选择特定列
  select(gene_name, chromosome, expression, log2_expression, 
         mutation_status, length_category)

cat("处理后的数据:\n")
print(processed_data)

# 分组汇总
summary_by_mutation <- gene_info %>%
  group_by(mutation_status) %>%
  summarise(
    count = n(),
    mean_expression = mean(expression),
    sd_expression = sd(expression),
    median_expression = median(expression),
    min_expression = min(expression),
    max_expression = max(expression),
    .groups = 'drop'
  )

cat("按突变状态分组的汇总统计:\n")
print(summary_by_mutation)

# 复杂的数据处理
complex_analysis <- gene_info %>%
  # 添加多个新列
  mutate(
    log2_expression = log2(expression + 1),
    expression_rank = rank(-expression),
    is_high_expression = expression > median(expression),
    chromosome_group = case_when(
      chromosome <= 10 ~ "low_chr",
      chromosome <= 20 ~ "mid_chr",
      TRUE ~ "high_chr"
    )
  ) %>%
  # 按染色体分组
  group_by(chromosome_group) %>%
  # 计算组内统计
  mutate(
    group_mean_expression = mean(expression),
    expression_vs_group_mean = expression - group_mean_expression
  ) %>%
  ungroup() %>%
  # 最终筛选和排序
  filter(gene_length > 50000) %>%
  arrange(chromosome, start_position)

cat("复杂数据处理结果:\n")
print(complex_analysis)
```

#### 4.2 描述性统计分析
```r
# 基本描述性统计
expression_values <- gene_info$expression

cat("基因表达水平描述性统计:\n")
cat("样本数:", length(expression_values), "\n")
cat("平均值:", mean(expression_values), "\n")
cat("中位数:", median(expression_values), "\n")
cat("标准差:", sd(expression_values), "\n")
cat("方差:", var(expression_values), "\n")
cat("最小值:", min(expression_values), "\n")
cat("最大值:", max(expression_values), "\n")
cat("范围:", diff(range(expression_values)), "\n")

# 分位数
quantiles <- quantile(expression_values, probs = c(0.25, 0.5, 0.75))
cat("分位数 (25%, 50%, 75%):", quantiles, "\n")

# 详细摘要
cat("详细统计摘要:\n")
print(summary(expression_values))

# 按分组的描述性统计
cat("按突变状态分组的描述性统计:\n")
by_mutation_detailed <- gene_info %>%
  group_by(mutation_status) %>%
  summarise(
    n = n(),
    mean = mean(expression),
    median = median(expression),
    sd = sd(expression),
    min = min(expression),
    max = max(expression),
    q25 = quantile(expression, 0.25),
    q75 = quantile(expression, 0.75),
    .groups = 'drop'
  )
print(by_mutation_detailed)

# 相关性分析
numeric_data <- gene_info %>%
  select(chromosome, start_position, end_position, gene_length, expression)

correlation_matrix <- cor(numeric_data)
cat("数值变量相关性矩阵:\n")
print(round(correlation_matrix, 3))
```

#### 4.3 假设检验
```r
# 准备测试数据
mutated_expression <- gene_info$expression[gene_info$mutation_status == "mutated"]
wildtype_expression <- gene_info$expression[gene_info$mutation_status == "wild_type"]

cat("突变基因表达水平:", mutated_expression, "\n")
cat("野生型基因表达水平:", wildtype_expression, "\n")

# 单样本t检验
# 检验突变基因的平均表达水平是否显著不同于5.0
single_t_test <- t.test(mutated_expression, mu = 5.0)
cat("单样本t检验结果 (突变基因 vs μ=5.0):\n")
cat("t统计量:", single_t_test$statistic, "\n")
cat("p值:", single_t_test$p.value, "\n")
cat("95%置信区间:", single_t_test$conf.int, "\n")

# 双样本t检验
# 检验突变基因和野生型基因的表达水平是否有显著差异
two_sample_t_test <- t.test(mutated_expression, wildtype_expression)
cat("双样本t检验结果 (突变 vs 野生型):\n")
cat("t统计量:", two_sample_t_test$statistic, "\n")
cat("p值:", two_sample_t_test$p.value, "\n")
cat("95%置信区间:", two_sample_t_test$conf.int, "\n")

# Wilcoxon秩和检验（非参数检验）
wilcox_test <- wilcox.test(mutated_expression, wildtype_expression)
cat("Wilcoxon秩和检验结果:\n")
cat("W统计量:", wilcox_test$statistic, "\n")
cat("p值:", wilcox_test$p.value, "\n")

# 相关性检验
cor_test <- cor.test(gene_info$gene_length, gene_info$expression)
cat("基因长度与表达水平相关性检验:\n")
cat("相关系数:", cor_test$estimate, "\n")
cat("p值:", cor_test$p.value, "\n")
cat("95%置信区间:", cor_test$conf.int, "\n")

# 正态性检验
shapiro_test <- shapiro.test(expression_values)
cat("表达水平正态性检验 (Shapiro-Wilk):\n")
cat("W统计量:", shapiro_test$statistic, "\n")
cat("p值:", shapiro_test$p.value, "\n")
```

**运行脚本**：`scripts/data_analysis.R`

### 步骤5：数据可视化

#### 5.1 基础绘图
```r
# 设置图形参数
par(mfrow = c(2, 2))  # 2x2的图形布局

# 1. 散点图
plot(gene_info$gene_length, gene_info$expression,
     xlab = "Gene Length (bp)", ylab = "Expression Level",
     main = "Gene Length vs Expression",
     pch = 19, col = "blue")
# 添加回归线
abline(lm(expression ~ gene_length, data = gene_info), col = "red", lwd = 2)

# 2. 柱状图
barplot(gene_info$expression, 
        names.arg = gene_info$gene_name,
        main = "Gene Expression Levels",
        ylab = "Expression Level",
        col = rainbow(nrow(gene_info)),
        las = 2)  # 旋转x轴标签

# 3. 箱线图
boxplot(expression ~ mutation_status, data = gene_info,
        main = "Expression by Mutation Status",
        xlab = "Mutation Status", ylab = "Expression Level",
        col = c("lightblue", "lightcoral"))

# 4. 直方图
hist(gene_info$expression,
     main = "Distribution of Expression Levels",
     xlab = "Expression Level", ylab = "Frequency",
     col = "lightgreen", breaks = 5)

# 重置图形参数
par(mfrow = c(1, 1))

# 保存图形
png("basic_plots.png", width = 800, height = 600)
par(mfrow = c(2, 2))
# 重复上面的绘图代码...
dev.off()
```

#### 5.2 ggplot2高级可视化
```r
library(ggplot2)
library(reshape2)

# 1. 基因表达柱状图
p1 <- ggplot(gene_info, aes(x = reorder(gene_name, -expression), y = expression)) +
  geom_bar(stat = "identity", aes(fill = mutation_status)) +
  scale_fill_manual(values = c("mutated" = "red", "wild_type" = "blue")) +
  labs(title = "Gene Expression Levels by Mutation Status",
       x = "Gene Name", y = "Expression Level",
       fill = "Mutation Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)

# 2. 散点图与回归线
p2 <- ggplot(gene_info, aes(x = gene_length, y = expression)) +
  geom_point(aes(color = mutation_status), size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = c("mutated" = "red", "wild_type" = "blue")) +
  labs(title = "Gene Length vs Expression Level",
       x = "Gene Length (bp)", y = "Expression Level",
       color = "Mutation Status") +
  theme_classic()

print(p2)

# 3. 箱线图
p3 <- ggplot(gene_info, aes(x = mutation_status, y = expression)) +
  geom_boxplot(aes(fill = mutation_status), alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2) +
  scale_fill_manual(values = c("mutated" = "red", "wild_type" = "blue")) +
  labs(title = "Expression Distribution by Mutation Status",
       x = "Mutation Status", y = "Expression Level") +
  theme_minimal() +
  theme(legend.position = "none")

print(p3)

# 4. 表达矩阵热图
# 准备热图数据
heatmap_data <- melt(expression_matrix)
colnames(heatmap_data) <- c("Gene", "Sample", "Expression")

p4 <- ggplot(heatmap_data, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = median(heatmap_data$Expression)) +
  labs(title = "Gene Expression Heatmap",
       x = "Sample", y = "Gene") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p4)

# 5. 多面板图
p5 <- ggplot(gene_info, aes(x = expression)) +
  geom_histogram(bins = 3, fill = "skyblue", alpha = 0.7) +
  facet_wrap(~mutation_status) +
  labs(title = "Expression Distribution by Mutation Status",
       x = "Expression Level", y = "Count") +
  theme_bw()

print(p5)

# 保存ggplot图形
ggsave("gene_expression_barplot.png", p1, width = 10, height = 6)
ggsave("length_vs_expression.png", p2, width = 8, height = 6)
ggsave("expression_boxplot.png", p3, width = 6, height = 6)
ggsave("expression_heatmap.png", p4, width = 8, height = 6)
```

#### 5.3 相关性和网络图
```r
library(corrplot)

# 相关性热图
numeric_data <- gene_info %>%
  select(chromosome, gene_length, expression) %>%
  as.matrix()

correlation_matrix <- cor(numeric_data)

# 使用corrplot绘制相关性矩阵
corrplot(correlation_matrix, method = "color", type = "upper",
         order = "hclust", tl.cex = 0.8, tl.col = "black",
         addCoef.col = "black", number.cex = 0.7)

# 保存相关性图
png("correlation_plot.png", width = 600, height = 600)
corrplot(correlation_matrix, method = "color", type = "upper",
         order = "hclust", tl.cex = 0.8, tl.col = "black",
         addCoef.col = "black", number.cex = 0.7)
dev.off()

# 成对散点图
pairs(numeric_data, main = "Pairwise Scatter Plots",
      pch = 19, col = "blue")
```

**运行脚本**：`scripts/visualization.R`

### 步骤6：生物信息学R包应用

#### 6.1 Biostrings包序列分析
```r
library(Biostrings)

# 创建DNA序列
dna_sequences <- DNAStringSet(c(
  gene1 = "ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCG",
  gene2 = "ATGCCCGGGAAATTTAAATTTCCCGGGAAATAGATCGATCGATCGATCGATC",
  gene3 = "ATGAAATTTCCCGGGAAATAGATCGATCGATCGATCGATCGATCGATCGATCG"
))

cat("DNA序列集合:\n")
print(dna_sequences)

# 序列基本信息
cat("序列长度:\n")
print(width(dna_sequences))

cat("序列名称:\n")
print(names(dna_sequences))

# 序列操作
cat("反向互补序列:\n")
rev_comp <- reverseComplement(dna_sequences)
print(rev_comp)

# 翻译为蛋白质序列
cat("翻译为蛋白质序列:\n")
protein_sequences <- translate(dna_sequences)
print(protein_sequences)

# 碱基频率分析
cat("碱基频率分析:\n")
base_freq <- alphabetFrequency(dna_sequences)
print(base_freq)

# GC含量计算
gc_content <- letterFrequency(dna_sequences, "GC", as.prob = TRUE)
cat("GC含量:\n")
print(gc_content)

# 二核苷酸频率
cat("二核苷酸频率:\n")
dinuc_freq <- dinucleotideFrequency(dna_sequences)
print(dinuc_freq[, 1:8])  # 显示前8列

# 模式匹配
cat("查找ATG起始密码子:\n")
atg_matches <- vmatchPattern("ATG", dna_sequences)
print(atg_matches)

# 计算CG二核苷酸数量
cg_counts <- vcountPattern("CG", dna_sequences)
cat("CG二核苷酸数量:\n")
print(cg_counts)
```

#### 6.2 GenomicRanges包基因组区间分析
```r
library(GenomicRanges)

# 创建基因组区间对象
gene_ranges <- GRanges(
  seqnames = paste0("chr", gene_info$chromosome),
  ranges = IRanges(
    start = gene_info$start_position,
    end = gene_info$end_position
  ),
  strand = "+",
  gene_name = gene_info$gene_name,
  expression = gene_info$expression,
  mutation_status = gene_info$mutation_status
)

cat("基因组区间对象:\n")
print(gene_ranges)

# 区间基本信息
cat("区间宽度 (基因长度):\n")
print(width(gene_ranges))

cat("起始位置:\n")
print(start(gene_ranges))

cat("结束位置:\n")
print(end(gene_ranges))

# 区间操作
# 获取基因上游2kb区域
promoter_regions <- flank(gene_ranges, width = 2000, start = TRUE)
cat("启动子区域 (上游2kb):\n")
print(promoter_regions)

# 调整基因区间大小
resized_genes <- resize(gene_ranges, width = 5000, fix = "center")
cat("调整为5kb的基因区间:\n")
print(resized_genes)

# 平移区间
shifted_genes <- shift(gene_ranges, 1000)
cat("向右平移1kb的基因区间:\n")
print(shifted_genes)

# 查找重叠
overlaps <- findOverlaps(gene_ranges, gene_ranges)
cat("基因间重叠分析:\n")
print(overlaps)

# 按染色体分组
by_chr <- split(gene_ranges, seqnames(gene_ranges))
cat("按染色体分组的基因数量:\n")
print(lengths(by_chr))
```

#### 6.3 综合生物信息学分析
```r
# 整合序列和基因组位置信息进行综合分析

# 创建分析报告函数
generate_gene_report <- function(gene_data, sequences) {
  report <- list()
  
  # 基本统计
  report$basic_stats <- list(
    total_genes = nrow(gene_data),
    chromosomes = unique(gene_data$chromosome),
    mean_expression = mean(gene_data$expression),
    mean_length = mean(gene_data$gene_length)
  )
  
  # 表达分析
  report$expression_analysis <- list(
    high_expression_genes = gene_data$gene_name[gene_data$expression > 5.0],
    low_expression_genes = gene_data$gene_name[gene_data$expression <= 5.0],
    expression_by_mutation = by(gene_data$expression, 
                               gene_data$mutation_status, mean)
  )
  
  # 序列分析（如果提供了序列）
  if (!missing(sequences)) {
    gc_contents <- letterFrequency(sequences, "GC", as.prob = TRUE)
    report$sequence_analysis <- list(
      mean_gc_content = mean(gc_contents),
      gc_range = range(gc_contents),
      sequences_with_high_gc = names(sequences)[gc_contents > 0.5]
    )
  }
  
  return(report)
}

# 生成报告
comprehensive_report <- generate_gene_report(gene_info, dna_sequences)

cat("综合基因分析报告:\n")
cat("===================\n")
cat("基本统计:\n")
print(comprehensive_report$basic_stats)

cat("\n表达分析:\n")
print(comprehensive_report$expression_analysis)

cat("\n序列分析:\n")
print(comprehensive_report$sequence_analysis)

# 创建综合可视化
library(gridExtra)

# 组合多个图形
p_expr <- ggplot(gene_info, aes(x = gene_name, y = expression)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression", x = "Gene", y = "Expression")

p_length <- ggplot(gene_info, aes(x = gene_name, y = gene_length)) +
  geom_bar(stat = "identity", fill = "orange") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Length", x = "Gene", y = "Length (bp)")

# 组合图形
combined_plot <- grid.arrange(p_expr, p_length, ncol = 2)

# 保存组合图形
ggsave("comprehensive_analysis.png", combined_plot, width = 12, height = 6)
```

**运行脚本**：`scripts/bioinformatics_analysis.R`

### 步骤7：文件操作和数据导入导出

#### 7.1 创建示例数据文件
```r
# 创建基因表达CSV文件
write.csv(gene_info, "gene_expression_output.csv", row.names = FALSE)

# 创建表达矩阵文件
write.table(expression_matrix, "expression_matrix.txt", 
            sep = "\t", quote = FALSE)

# 创建样本信息文件
sample_info <- data.frame(
  sample_id = colnames(expression_matrix),
  condition = c("control", "control", "treatment", "treatment"),
  batch = c(1, 1, 2, 2),
  stringsAsFactors = FALSE
)

write.table(sample_info, "sample_info.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("数据文件创建完成:\n")
cat("- gene_expression_output.csv\n")
cat("- expression_matrix.txt\n")
cat("- sample_info.txt\n")
```

#### 7.2 读取和验证数据文件
```r
# 读取CSV文件
gene_data_read <- read.csv("gene_expression_output.csv", stringsAsFactors = FALSE)
cat("读取的基因数据:\n")
print(head(gene_data_read))

# 读取表达矩阵
expr_matrix_read <- read.table("expression_matrix.txt", 
                               header = TRUE, row.names = 1)
cat("读取的表达矩阵:\n")
print(expr_matrix_read)

# 读取样本信息
sample_data_read <- read.table("sample_info.txt", 
                               header = TRUE, stringsAsFactors = FALSE)
cat("读取的样本信息:\n")
print(sample_data_read)

# 数据验证
cat("数据验证:\n")
cat("基因数据行数:", nrow(gene_data_read), "\n")
cat("表达矩阵维度:", dim(expr_matrix_read), "\n")
cat("样本信息行数:", nrow(sample_data_read), "\n")

# 检查数据一致性
genes_match <- all(gene_data_read$gene_name %in% rownames(expr_matrix_read))
samples_match <- all(colnames(expr_matrix_read) %in% sample_data_read$sample_id)

cat("基因名称匹配:", genes_match, "\n")
cat("样本名称匹配:", samples_match, "\n")
```

## 预期结果

### 步骤1预期结果
- R环境成功配置，版本4.0+
- 所有必需包安装完成
- RStudio环境正常工作

### 步骤2预期结果
- 基本数据类型操作正确
- 向量运算和索引功能正常
- 字符串处理和序列分析准确

### 步骤3预期结果
- 矩阵操作和统计计算正确
- 数据框筛选和排序功能正常
- 列表操作和数据结构理解准确

### 步骤4预期结果
- dplyr数据处理流程顺畅
- 描述性统计计算准确
- 假设检验结果合理

### 步骤5预期结果
- 基础图形绘制清晰
- ggplot2图表专业美观
- 相关性可视化准确

### 步骤6预期结果
- Biostrings序列分析正确
- GenomicRanges区间操作准确
- 生物信息学分析流程完整

### 步骤7预期结果
- 数据文件成功创建和读取
- 文件格式正确
- 数据一致性验证通过

## 故障排除

### 常见问题1：包安装失败
**症状**：install.packages()报错
**解决方案**：
```r
# 更换CRAN镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 或者手动选择镜像
chooseCRANmirror()

# 检查网络连接
capabilities()
```

### 常见问题2：Bioconductor包安装问题
**症状**：BiocManager::install()失败
**解决方案**：
```r
# 更新BiocManager
install.packages("BiocManager")

# 检查Bioconductor版本
BiocManager::version()

# 使用特定版本
BiocManager::install(version = "3.14")
```

### 常见问题3：中文字符显示问题
**症状**：图表中中文显示异常
**解决方案**：
```r
# 设置字符编码
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

# 在ggplot中设置字体
theme_set(theme_minimal(base_family = "Arial Unicode MS"))
```

### 常见问题4：内存不足
**症状**：处理大数据时R崩溃
**解决方案**：
```r
# 检查内存使用
memory.limit()
memory.size()

# 清理内存
gc()
rm(list = ls())

# 增加内存限制（Windows）
memory.limit(size = 8000)  # 8GB
```

### 常见问题5：工作目录问题
**症状**：找不到文件
**解决方案**：
```r
# 检查当前工作目录
getwd()

# 设置工作目录
setwd("/path/to/your/directory")

# 列出当前目录文件
list.files()

# 使用相对路径
file.exists("data/gene_expression.csv")
```

## 扩展练习

### 练习1：基因表达数据分析
使用提供的基因表达数据完成以下任务：
1. 计算每个基因在不同样本间的变异系数
2. 识别表达模式相似的基因对
3. 进行层次聚类分析
4. 创建基因表达热图

### 练习2：统计分析项目
设计并执行以下统计分析：
1. 比较不同突变状态基因的表达差异
2. 分析基因长度与表达水平的关系
3. 进行多重比较校正
4. 计算效应大小和统计功效

### 练习3：数据可视化挑战
创建以下高级可视化图表：
1. 交互式散点图（使用plotly）
2. 基因网络图（使用igraph）
3. 基因组位置图（使用karyoploteR）
4. 多维数据降维可视化（PCA/t-SNE）

### 练习4：生物信息学流水线
构建一个完整的分析流水线：
1. 从FASTA文件读取序列
2. 计算序列特征（长度、GC含量、密码子使用）
3. 整合表达数据进行关联分析
4. 生成自动化分析报告

## 参考资料

### 官方文档
- [R官方网站](https://www.r-project.org/)
- [RStudio文档](https://rstudio.com/resources/)
- [CRAN包文档](https://cran.r-project.org/)
- [Bioconductor文档](https://bioconductor.org/help/)

### R包文档
- [tidyverse](https://www.tidyverse.org/)
- [ggplot2](https://ggplot2.tidyverse.org/)
- [dplyr](https://dplyr.tidyverse.org/)
- [Biostrings](https://bioconductor.org/packages/Biostrings/)
- [GenomicRanges](https://bioconductor.org/packages/GenomicRanges/)

### 推荐书籍
- 《R语言实战》
- 《ggplot2数据分析与图形艺术》
- 《R数据科学》
- 《Bioconductor案例研究》

### 在线资源
- [R for Data Science](https://r4ds.had.co.nz/)
- [Advanced R](https://adv-r.hadley.nz/)
- [R Graphics Cookbook](https://r-graphics.org/)

## 课后作业

1. **基础练习**：完成所有步骤中的R代码练习
2. **数据分析项目**：使用提供的基因表达数据进行完整的探索性数据分析
3. **可视化作业**：创建至少5种不同类型的生物信息学相关图表
4. **编程挑战**：编写一个R函数来自动化基因表达数据的质量控制流程

**提交要求**：
- 所有R脚本文件（.R格式）
- 分析结果和图表文件
- 详细的分析报告（包含代码、结果和解释）

**截止时间**：下次课前

---

**注意事项**：
- 保存好所有R脚本，后续课程会继续使用
- 熟悉RStudio的使用，提高编程效率
- 多练习数据处理和可视化技能
- 关注R代码的规范性和可读性
- 学会使用R的帮助系统和在线资源