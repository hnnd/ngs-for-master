---
marp: true
theme: ngs-course
paginate: true
header: '高通量测序数据分析'
footer: '王运生 | 2025'
---

<!-- 
R语言基础课程幻灯片
课程名称：高通量测序数据分析 - 预备课程
主讲教师：王运生教授
联系邮箱：wangys@hunau.edu.cn
办公室：16教420室
上课地点：105机房
-->

<!-- _class: title -->
# R语言基础
## 高通量测序数据分析 - 预备课程

**主讲教师：** 王运生  
**联系邮箱：** wangys@hunau.edu.cn  
**办公室：** 16教420室  
**上课地点：** 105机房  

---

<!-- _class: toc -->
# 本次课程内容

1. R语言简介与环境配置
2. 基础语法与数据类型
3. 数据结构详解
4. 数据处理与操作
5. 统计分析基础
6. 数据可视化
7. 生物信息学R包

---

<!-- _class: content -->
# 学习目标

**本次课程结束后，您将能够：**
- 掌握R语言基础编程技能
- 熟练使用R进行数据处理和统计分析
- 创建专业的数据可视化图表
- 了解R在生物信息学中的应用
- 使用R包进行生物数据分析

---

<!-- _class: content -->
# R语言简介


## 什么是R？

- **统计计算语言**：专为数据分析和统计建模设计
- **开源免费**：GNU GPL许可证，完全免费使用
- **跨平台**：Windows、Linux、macOS都支持
- **丰富的包生态**：CRAN、Bioconductor等包仓库

---

<!-- _class: content -->
# 为什么生物信息学选择R？

## R的优势

- **统计分析强大**：内置丰富的统计函数
- **数据可视化优秀**：ggplot2等强大绘图系统
- **生物信息学生态**：Bioconductor项目提供专业工具
- **活跃的社区**：大量教程、文档和支持

---

<!-- _class: multi-column -->
# R在生物信息学中的应用

<div class="columns">
<div class="column">

## 数据分析
- **基因表达分析**：RNA-seq、微阵列数据
- **基因组学**：变异分析、注释
- **蛋白质组学**：质谱数据分析
- **统计建模**：回归、分类、聚类

</div>
<div class="column">

## 可视化
- **基因表达热图**：ComplexHeatmap
- **基因组浏览**：Gviz、karyoploteR
- **网络分析**：igraph、visNetwork
- **统计图表**：ggplot2生态系统

</div>
</div>

---

<!-- _class: content -->
# R环境配置

## 安装R

### 1. 下载安装R
```bash
# Ubuntu/Debian
sudo apt update
sudo apt install r-base r-base-dev

# CentOS/RHEL
sudo yum install R

# 或从官网下载：https://www.r-project.org/
```

---

### 2. 安装RStudio（推荐）
- 下载地址：https://rstudio.com/products/rstudio/download/
- 提供友好的集成开发环境
- 支持代码编辑、调试、包管理

---

<!-- _class: code -->
# 第一个R程序

```r
# 在R控制台中执行
print("Hello, 生物信息学世界！")

# 基本计算
2 + 3
10 * 5
sqrt(16)

# 赋值操作
x <- 5
y = 10
result <- x + y
print(result)

# 获取帮助
?print
help(mean)
```r

---

<!-- _class: content -->
# R基础语法

## 语法特点

- **大小写敏感**：`Variable`和`variable`是不同的
- **赋值操作**：使用`<-`或`=`（推荐使用`<-`）
- **注释**：使用`#`进行单行注释
- **向量化操作**：大多数函数支持向量输入

---

<!-- _class: code -->
# 基础语法示例

```r
# 这是注释
# R中的赋值
gene_name <- "BRCA1"
expression_level <- 5.2

# 向量化操作
numbers <- c(1, 2, 3, 4, 5)
squared <- numbers^2
print(squared)  # 输出: 1 4 9 16 25

# 逻辑操作
high_expression <- expression_level > 5.0
print(high_expression)  # 输出: TRUE
```

---

<!-- _class: multi-column -->
# 数据类型

<div class="columns">
<div class="column">

## 基本数据类型
```r
# 数值型 (numeric)
age <- 25
gc_content <- 0.42

# 字符型 (character)
gene_name <- "BRCA1"
sequence <- "ATCGATCG"

# 逻辑型 (logical)
is_coding <- TRUE
has_mutation <- FALSE

# 整数型 (integer)
chromosome <- 17L
```

</div>

<div class="column">

## 类型检查和转换
```r
# 查看类型
class(42)           # "numeric"
typeof("DNA")       # "character"
is.numeric(3.14)    # TRUE
is.character("A")   # TRUE

# 类型转换
as.numeric("123")   # 123
as.character(42)    # "42"
as.logical(1)       # TRUE
```

</div>
</div>

---

<!-- _class: content -->
# 向量 (Vector)

## 向量基础操作

```r
# 创建向量
genes <- c("BRCA1", "BRCA2", "TP53", "EGFR")
expression <- c(5.2, 3.8, 7.1, 4.5)
chromosomes <- c(17, 13, 17, 7)

# 访问元素
first_gene <- genes[1]        # "BRCA1" (R索引从1开始)
last_gene <- genes[length(genes)]  # "EGFR"
subset <- genes[2:3]          # "BRCA2" "TP53"

# 向量运算
log_expression <- log2(expression)
normalized <- (expression - mean(expression)) / sd(expression)

# 逻辑索引
high_expr_genes <- genes[expression > 5.0]
print(high_expr_genes)  # "BRCA1" "TP53"
```

---

<!-- _class: code -->
# 向量高级操作

```r
# 命名向量
named_expr <- c(BRCA1=5.2, BRCA2=3.8, TP53=7.1, EGFR=4.5)
print(named_expr["BRCA1"])  # 5.2

# 向量函数
length(expression)    # 向量长度
sum(expression)       # 求和
mean(expression)      # 平均值
median(expression)    # 中位数
sd(expression)        # 标准差
range(expression)     # 范围

# 向量合并
all_data <- c(expression, c(6.3, 4.8))
print(all_data)

# 重复和序列
rep(1, 5)            # 重复: 1 1 1 1 1
seq(1, 10, by=2)     # 序列: 1 3 5 7 9
1:10                 # 简写: 1 2 3 4 5 6 7 8 9 10
```

---

<!-- _class: content -->
# 矩阵 (Matrix)

## 矩阵基础操作

```r
# 创建矩阵
expression_matrix <- matrix(c(5.2, 3.8, 7.1, 4.5, 6.3, 4.8), 
                           nrow=2, ncol=3, byrow=TRUE)

# 添加行名和列名
rownames(expression_matrix) <- c("Sample1", "Sample2")
colnames(expression_matrix) <- c("BRCA1", "BRCA2", "TP53")

print(expression_matrix)
#         BRCA1 BRCA2 TP53
# Sample1   5.2   3.8  7.1
# Sample2   4.5   6.3  4.8

# 访问矩阵元素
expression_matrix[1, 2]      # 第1行第2列
expression_matrix[1, ]       # 第1行所有列
expression_matrix[, "BRCA1"] # BRCA1列
```

---

<!-- _class: code -->
# 矩阵运算

```r
# 矩阵统计
rowMeans(expression_matrix)  # 每行平均值
colMeans(expression_matrix)  # 每列平均值
apply(expression_matrix, 1, sd)  # 每行标准差
apply(expression_matrix, 2, max) # 每列最大值

# 矩阵运算
t(expression_matrix)         # 转置
expression_matrix * 2        # 标量乘法
log2(expression_matrix)      # 对数变换

# 矩阵合并
new_sample <- c(4.9, 5.1, 6.8)
rbind(expression_matrix, Sample3=new_sample)  # 按行合并

new_gene <- c(3.5, 4.2)
cbind(expression_matrix, MYC=new_gene)        # 按列合并
```

---

<!-- _class: content -->
# 数据框 (Data Frame)

## 数据框基础操作

```r
# 创建数据框
gene_data <- data.frame(
  gene_name = c("BRCA1", "BRCA2", "TP53", "EGFR"),
  chromosome = c(17, 13, 17, 7),
  expression = c(5.2, 3.8, 7.1, 4.5),
  mutation_status = c("mutated", "wild_type", "mutated", "wild_type"),
  stringsAsFactors = FALSE  # 避免自动转换为因子
)

# 查看数据框
head(gene_data)      # 前几行
str(gene_data)       # 数据结构
summary(gene_data)   # 统计摘要
dim(gene_data)       # 维度
names(gene_data)     # 列名
```

---

<!-- _class: code -->
# 数据框操作

```r
# 访问列
gene_data$gene_name              # 使用$符号
gene_data[["expression"]]        # 使用[[]]
gene_data[, "chromosome"]        # 使用索引

# 访问行
gene_data[1, ]                   # 第1行
gene_data[gene_data$expression > 5, ]  # 条件筛选

# 添加新列
gene_data$log_expression <- log2(gene_data$expression)
gene_data$high_expression <- gene_data$expression > 5.0

# 排序
gene_data[order(gene_data$expression), ]           # 按表达量升序
gene_data[order(-gene_data$expression), ]          # 按表达量降序

# 子集选择
subset(gene_data, expression > 5.0)               # 高表达基因
subset(gene_data, chromosome == 17)               # 17号染色体基因
```

---

<!-- _class: content -->
# 列表 (List)

## 列表基础操作

```r
# 创建列表
analysis_results <- list(
  genes = c("BRCA1", "BRCA2", "TP53"),
  expression_matrix = matrix(1:9, nrow=3),
  metadata = data.frame(
    sample_id = c("S1", "S2", "S3"),
    condition = c("control", "treatment", "control")
  ),
  parameters = list(
    p_value_cutoff = 0.05,
    fold_change_cutoff = 2.0
  )
)

# 访问列表元素
analysis_results$genes                    # 使用$
analysis_results[["expression_matrix"]]   # 使用[[]]
analysis_results[[1]]                     # 使用数字索引
```

---

<!-- _class: code -->
# 列表高级操作

```r
# 列表信息
length(analysis_results)     # 列表长度
names(analysis_results)      # 元素名称
str(analysis_results)        # 结构信息

# 添加和删除元素
analysis_results$new_data <- c(1, 2, 3)  # 添加
analysis_results$new_data <- NULL         # 删除

# 列表合并
list1 <- list(a=1, b=2)
list2 <- list(c=3, d=4)
combined <- c(list1, list2)

# 应用函数到列表
numeric_list <- list(a=1:5, b=6:10, c=11:15)
lapply(numeric_list, mean)    # 返回列表
sapply(numeric_list, mean)    # 返回向量
```

---

<!-- _class: content -->
# 因子 (Factor)

## 因子基础操作

```r
# 创建因子
treatment <- factor(c("control", "treatment", "control", "treatment"))
grade <- factor(c("low", "medium", "high", "medium"), 
                levels=c("low", "medium", "high"), 
                ordered=TRUE)

print(treatment)
# [1] control   treatment control   treatment
# Levels: control treatment

print(grade)
# [1] low    medium high   medium
# Levels: low < medium < high

# 因子信息
levels(treatment)    # 水平
nlevels(treatment)   # 水平数量
table(treatment)     # 频数统计
```

---

<!-- _class: code -->
# 因子操作

```r
# 重新编码因子
treatment_coded <- factor(treatment, 
                         levels=c("control", "treatment"),
                         labels=c("对照组", "处理组"))

# 添加水平
levels(treatment) <- c(levels(treatment), "new_treatment")

# 合并因子水平
treatment_simple <- factor(treatment, 
                          levels=c("control", "treatment", "new_treatment"),
                          labels=c("control", "treatment", "treatment"))

# 因子在统计中的应用
# 用于分组分析
by(gene_data$expression, gene_data$mutation_status, mean)
```

---

<!-- _class: content -->
# 数据读取和写入

## 读取数据

```r
# 读取CSV文件
gene_expr <- read.csv("gene_expression.csv", header=TRUE)
gene_expr <- read.table("data.txt", sep="\t", header=TRUE)

# 使用readr包（更快，更智能）
library(readr)
gene_expr <- read_csv("gene_expression.csv")
gene_expr <- read_tsv("data.tsv")

# 读取Excel文件
library(readxl)
gene_expr <- read_excel("data.xlsx", sheet=1)

# 读取其他格式
load("data.RData")           # R数据文件
gene_expr <- readRDS("data.rds")  # R对象文件
```

---

<!-- _class: code -->
# 写入数据

```r
# 写入CSV文件
write.csv(gene_data, "output.csv", row.names=FALSE)
write.table(gene_data, "output.txt", sep="\t", row.names=FALSE)

# 使用readr包
library(readr)
write_csv(gene_data, "output.csv")
write_tsv(gene_data, "output.tsv")

# 保存R对象
save(gene_data, file="gene_data.RData")
saveRDS(gene_data, "gene_data.rds")

# 写入Excel文件
library(writexl)
write_xlsx(gene_data, "output.xlsx")
```

---

<!-- _class: content -->
# 数据处理 - dplyr包

## dplyr基础操作

```r
library(dplyr)

# 主要动词函数
# filter()  - 筛选行
# select()  - 选择列
# mutate()  - 创建新列
# arrange() - 排序
# summarise() - 汇总
# group_by() - 分组

# 示例数据
gene_data <- data.frame(
  gene_name = c("BRCA1", "BRCA2", "TP53", "EGFR", "MYC"),
  chromosome = c(17, 13, 17, 7, 8),
  expression = c(5.2, 3.8, 7.1, 4.5, 6.3),
  mutation_status = c("mutated", "wild_type", "mutated", "wild_type", "mutated")
)
```

---

<!-- _class: code -->
# dplyr操作示例

```r
# 筛选高表达基因
high_expr <- gene_data %>%
  filter(expression > 5.0)

# 选择特定列
selected <- gene_data %>%
  select(gene_name, expression)

# 创建新列
gene_data <- gene_data %>%
  mutate(
    log_expression = log2(expression),
    expression_category = ifelse(expression > 5.0, "high", "low")
  )

# 排序
sorted <- gene_data %>%
  arrange(desc(expression))

# 分组汇总
summary_stats <- gene_data %>%
  group_by(mutation_status) %>%
  summarise(
    count = n(),
    mean_expr = mean(expression),
    sd_expr = sd(expression)
  )
```

---

<!-- _class: content -->
# 统计分析基础

## 描述性统计

```r
# 基本统计函数
expression <- c(5.2, 3.8, 7.1, 4.5, 6.3, 4.8, 5.9)

mean(expression)      # 平均值
median(expression)    # 中位数
sd(expression)        # 标准差
var(expression)       # 方差
min(expression)       # 最小值
max(expression)       # 最大值
range(expression)     # 范围
quantile(expression)  # 分位数

# 汇总统计
summary(expression)
```

---

<!-- _class: code -->
# 假设检验

```r
# t检验
control <- c(4.2, 3.8, 4.5, 3.9, 4.1)
treatment <- c(6.1, 5.8, 6.3, 5.9, 6.0)

# 单样本t检验
t.test(control, mu=4.0)

# 双样本t检验
t.test(control, treatment)

# 配对t检验
before <- c(4.2, 3.8, 4.5, 3.9, 4.1)
after <- c(4.8, 4.2, 5.1, 4.3, 4.6)
t.test(before, after, paired=TRUE)

# 方差分析
groups <- factor(rep(c("A", "B", "C"), each=5))
values <- c(rnorm(5, 4), rnorm(5, 5), rnorm(5, 6))
anova_result <- aov(values ~ groups)
summary(anova_result)
```

---

<!-- _class: content -->
# 相关性和回归分析

```r
# 相关性分析
x <- c(1, 2, 3, 4, 5)
y <- c(2.1, 3.9, 6.2, 7.8, 10.1)

cor(x, y)                    # 皮尔逊相关系数
cor(x, y, method="spearman") # 斯皮尔曼相关系数
cor.test(x, y)              # 相关性检验

# 线性回归
model <- lm(y ~ x)
summary(model)

# 预测
new_x <- data.frame(x=c(6, 7))
predict(model, new_x)

# 多元回归
gene_expr <- c(5.2, 3.8, 7.1, 4.5, 6.3)
gene_length <- c(2300, 1800, 2100, 1900, 2500)
gc_content <- c(0.45, 0.38, 0.52, 0.41, 0.48)

multi_model <- lm(gene_expr ~ gene_length + gc_content)
summary(multi_model)
```

---

<!-- _class: content -->
# 数据可视化 - 基础绘图

## 基础图形函数

```r
# 散点图
plot(gene_length, gene_expr, 
     xlab="Gene Length", ylab="Expression Level",
     main="Gene Length vs Expression")

# 柱状图
barplot(gene_expr, names.arg=c("BRCA1", "BRCA2", "TP53", "EGFR", "MYC"),
        main="Gene Expression Levels")

# 直方图
hist(rnorm(1000), main="Normal Distribution", xlab="Value")

# 箱线图
boxplot(gene_expr ~ mutation_status, data=gene_data,
        main="Expression by Mutation Status")

# 线图
plot(1:10, (1:10)^2, type="l", main="Quadratic Function")
```r

---

<!-- _class: content -->
# 数据可视化 - ggplot2

## ggplot2基础

```r
library(ggplot2)

# 基本语法：数据 + 映射 + 几何对象
ggplot(data=gene_data, aes(x=gene_name, y=expression)) +
  geom_bar(stat="identity")

# 散点图
ggplot(gene_data, aes(x=chromosome, y=expression)) +
  geom_point(size=3) +
  labs(title="Expression vs Chromosome",
       x="Chromosome", y="Expression Level")

# 按组着色
ggplot(gene_data, aes(x=gene_name, y=expression, fill=mutation_status)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(title="Gene Expression by Mutation Status")
```r

---

<!-- _class: code -->
# ggplot2高级图形

```r
# 箱线图
ggplot(gene_data, aes(x=mutation_status, y=expression)) +
  geom_boxplot() +
  geom_jitter(width=0.2) +
  theme_classic()

# 热图
library(reshape2)
expr_matrix <- matrix(rnorm(50), nrow=5)
rownames(expr_matrix) <- paste0("Gene", 1:5)
colnames(expr_matrix) <- paste0("Sample", 1:10)

melted <- melt(expr_matrix)
ggplot(melted, aes(x=Var2, y=Var1, fill=value)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red") +
  theme_minimal() +
  labs(title="Gene Expression Heatmap")

# 多面板图
ggplot(gene_data, aes(x=expression)) +
  geom_histogram(bins=10) +
  facet_wrap(~mutation_status) +
  theme_bw()
```

---

<!-- _class: content -->
# 生物信息学R包 - Bioconductor

## Bioconductor简介

- **专业的生物信息学R包集合**
- **高质量的包**：严格的审查和文档要求
- **统一的数据结构**：S4对象系统
- **丰富的功能**：基因组学、转录组学、蛋白质组学

```r
# 安装Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# 安装生物信息学包
BiocManager::install("Biostrings")
BiocManager::install("GenomicRanges")
BiocManager::install("DESeq2")
```

---

<!-- _class: code -->
# Biostrings包 - 序列分析

```r
library(Biostrings)

# 创建DNA序列
dna_seq <- DNAString("ATCGATCGATCG")
dna_set <- DNAStringSet(c(
  seq1="ATCGATCGATCG",
  seq2="GCTAGCTAGCTA",
  seq3="TTAACCGGTTAA"
))

# 序列操作
length(dna_seq)                    # 序列长度
reverseComplement(dna_seq)         # 反向互补
translate(dna_seq)                 # 翻译

# 序列分析
alphabetFrequency(dna_seq)         # 碱基频率
dinucleotideFrequency(dna_seq)     # 二核苷酸频率

# 模式匹配
matchPattern("CG", dna_seq)        # 查找CG
countPattern("AT", dna_seq)        # 计数AT

# GC含量
letterFrequency(dna_seq, "GC", as.prob=TRUE)
```

---

<!-- _class: content -->
# GenomicRanges包 - 基因组区间

```r
library(GenomicRanges)

# 创建基因组区间
gr <- GRanges(
  seqnames = c("chr1", "chr1", "chr2"),
  ranges = IRanges(start=c(100, 200, 300), end=c(150, 250, 350)),
  strand = c("+", "-", "+"),
  gene_id = c("gene1", "gene2", "gene3")
)

# 区间操作
width(gr)                          # 区间宽度
start(gr)                          # 起始位置
end(gr)                            # 结束位置
strand(gr)                         # 链方向

# 区间运算
resize(gr, width=100)              # 调整大小
shift(gr, 50)                      # 平移
flank(gr, width=20)                # 侧翼区域

# 重叠分析
findOverlaps(gr, gr)               # 查找重叠
```

---

<!-- _class: content -->
# 实际应用案例

## 基因表达数据分析流程

```r
# 1. 数据读取和预处理
library(tidyverse)

# 读取数据
expr_data <- read_csv("gene_expression.csv")
sample_info <- read_csv("sample_info.csv")

# 数据清理
expr_clean <- expr_data %>%
  filter(!is.na(expression)) %>%
  mutate(log2_expr = log2(expression + 1))

# 2. 探索性数据分析
summary(expr_clean$expression)
hist(expr_clean$log2_expr)

# 3. 差异表达分析
control_samples <- sample_info %>% filter(condition == "control") %>% pull(sample_id)
treatment_samples <- sample_info %>% filter(condition == "treatment") %>% pull(sample_id)
```

---

<!-- _class: code -->
# 基因表达分析（续）

```r
# 4. 统计检验
perform_ttest <- function(gene_data, control_ids, treatment_ids) {
  control_expr <- gene_data[gene_data$sample_id %in% control_ids, "expression"]
  treatment_expr <- gene_data[gene_data$sample_id %in% treatment_ids, "expression"]
  
  test_result <- t.test(control_expr, treatment_expr)
  
  return(data.frame(
    p_value = test_result$p.value,
    fold_change = mean(treatment_expr) / mean(control_expr),
    control_mean = mean(control_expr),
    treatment_mean = mean(treatment_expr)
  ))
}

# 5. 可视化结果
library(ggplot2)

# 火山图
results %>%
  mutate(
    log2_fc = log2(fold_change),
    neg_log10_p = -log10(p_value),
    significant = p_value < 0.05 & abs(log2_fc) > 1
  ) %>%
  ggplot(aes(x=log2_fc, y=neg_log10_p, color=significant)) +
  geom_point() +
  scale_color_manual(values=c("grey", "red")) +
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-Log10 P-value")
```

---

<!-- _class: content -->
# 包管理和环境

## 包的安装和加载

```r
# 安装包
install.packages("ggplot2")        # 从CRAN安装
BiocManager::install("DESeq2")     # 从Bioconductor安装
devtools::install_github("user/package")  # 从GitHub安装

# 加载包
library(ggplot2)
require(dplyr)

# 检查包
installed.packages()               # 已安装的包
.libPaths()                       # 包安装路径
sessionInfo()                     # 会话信息

# 更新包
update.packages()                 # 更新所有包
BiocManager::install()            # 更新Bioconductor包
```

---

<!-- _class: content -->
# 调试和最佳实践

## 调试技巧

```r
# 调试函数
debug(my_function)        # 进入调试模式
undebug(my_function)      # 退出调试模式
browser()                 # 在代码中设置断点

# 错误处理
tryCatch({
  result <- risky_operation()
}, error = function(e) {
  print(paste("Error:", e$message))
  return(NULL)
})

# 代码性能
system.time({
  # 要测试的代码
})

# 内存使用
object.size(large_object)
```

---

<!-- _class: content -->
# 代码规范

## R编程最佳实践

```r
# 好的命名规范
gene_expression_data <- data.frame()  # 变量名用下划线
calculate_gc_content <- function() {} # 函数名用下划线

# 代码组织
library(tidyverse)  # 在脚本开头加载包
library(ggplot2)

# 函数文档
#' Calculate GC content of DNA sequence
#' 
#' @param sequence A character string representing DNA sequence
#' @return Numeric value representing GC content proportion
#' @examples
#' calculate_gc_content("ATCGATCG")
calculate_gc_content <- function(sequence) {
  # 函数实现
}

# 使用管道操作符提高可读性
result <- data %>%
  filter(condition) %>%
  mutate(new_column = calculation) %>%
  arrange(sort_column)
```

---

<!-- _class: summary -->
# 本次课程总结

## 主要内容回顾
- R语言基础语法和数据类型
- 数据结构：向量、矩阵、数据框、列表
- 数据处理和统计分析
- 数据可视化：基础绘图和ggplot2
- 生物信息学R包介绍

---

## 重要概念
- **向量化操作**：R的核心特性，提高计算效率
- **数据框**：最常用的数据结构，类似Excel表格
- **管道操作**：使用%>%提高代码可读性
- **ggplot2语法**：图层语法，灵活强大的可视化
- **Bioconductor**：生物信息学专业包集合

---

## 下次课程预告
- 高通量测序技术原理
- 测序平台比较
- 数据格式和质量控制

**作业/练习：**
- 完成R基础编程练习
- 使用R分析提供的基因表达数据
- 创建数据可视化图表

---

<!-- _class: end -->
# 谢谢大家！

**有问题请联系：**
- 邮箱：wangys@hunau.edu.cn
- 办公室：16教420室

**推荐学习资源：**
- R官方网站：https://www.r-project.org/
- RStudio教程：https://education.rstudio.com/
- Bioconductor：https://bioconductor.org/
- 《R语言实战》、《ggplot2数据分析与图形艺术》