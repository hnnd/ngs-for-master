#!/usr/bin/env Rscript
# R语言基础课程 - 数据可视化脚本
# 作者：王运生
# 日期：2025-01-20
# 用途：基础绘图和ggplot2高级可视化

# 加载必需的包
library(ggplot2)
library(dplyr)
library(reshape2)
library(corrplot)
library(RColorBrewer)
library(gridExtra)

# 清理环境
rm(list = ls())

cat("R语言数据可视化练习\n")
cat("====================\n\n")

# 读取或创建数据
if (file.exists("gene_expression.csv")) {
  gene_info <- read.csv("gene_expression.csv", stringsAsFactors = FALSE)
} else {
  # 创建示例数据
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

# 创建表达矩阵数据
set.seed(123)
gene_names <- gene_info$gene_name
sample_names <- c("Control1", "Control2", "Treatment1", "Treatment2")

expression_matrix <- matrix(
  round(abs(rnorm(20, mean = 5, sd = 1.5)), 2),
  nrow = 5, ncol = 4,
  dimnames = list(gene_names, sample_names)
)

cat("开始数据可视化练习...\n\n")

# 练习1: 基础绘图系统
cat("练习1: 基础绘图系统\n")
cat("-------------------\n")

# 设置图形输出
png("basic_plots_combined.png", width = 1200, height = 900)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 2))

# 1. 散点图
plot(gene_info$gene_length, gene_info$expression,
     xlab = "基因长度 (bp)", ylab = "表达水平",
     main = "基因长度 vs 表达水平",
     pch = 19, col = "blue", cex = 1.5)
# 添加回归线
abline(lm(expression ~ gene_length, data = gene_info), col = "red", lwd = 2)
# 添加相关系数
cor_coef <- cor(gene_info$gene_length, gene_info$expression)
text(min(gene_info$gene_length), max(gene_info$expression), 
     paste("r =", round(cor_coef, 3)), pos = 4)

# 2. 柱状图
colors <- rainbow(nrow(gene_info))
barplot(gene_info$expression, 
        names.arg = gene_info$gene_name,
        main = "基因表达水平",
        ylab = "表达水平",
        col = colors,
        las = 2,
        cex.names = 0.8)

# 3. 箱线图
mutation_colors <- c("lightblue", "lightcoral")
boxplot(expression ~ mutation_status, data = gene_info,
        main = "按突变状态分组的表达水平",
        xlab = "突变状态", ylab = "表达水平",
        col = mutation_colors,
        names = c("突变", "野生型"))
# 添加数据点
stripchart(expression ~ mutation_status, data = gene_info,
           vertical = TRUE, method = "jitter", add = TRUE, pch = 19, cex = 1.2)

# 4. 直方图
hist(gene_info$expression,
     main = "表达水平分布",
     xlab = "表达水平", ylab = "频数",
     col = "lightgreen", breaks = 5,
     border = "black")
# 添加密度曲线
lines(density(gene_info$expression), col = "red", lwd = 2)

# 5. 饼图
mutation_table <- table(gene_info$mutation_status)
pie(mutation_table,
    main = "突变状态分布",
    col = c("lightblue", "lightcoral"),
    labels = paste(names(mutation_table), "\n", mutation_table, "个"))

# 6. 线图
plot(1:nrow(gene_info), sort(gene_info$expression),
     type = "b", pch = 19, col = "purple",
     main = "表达水平排序",
     xlab = "基因排名", ylab = "表达水平",
     lwd = 2)
grid()

dev.off()
cat("基础图形已保存为 'basic_plots_combined.png'\n\n")

# 练习2: ggplot2基础可视化
cat("练习2: ggplot2基础可视化\n")
cat("-----------------------\n")

# 1. 基因表达柱状图
p1 <- ggplot(gene_info, aes(x = reorder(gene_name, -expression), y = expression)) +
  geom_bar(stat = "identity", aes(fill = mutation_status), alpha = 0.8) +
  scale_fill_manual(values = c("mutated" = "#E74C3C", "wild_type" = "#3498DB"),
                    labels = c("突变", "野生型")) +
  labs(title = "基因表达水平（按突变状态分组）",
       x = "基因名称", y = "表达水平",
       fill = "突变状态") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "top"
  ) +
  geom_text(aes(label = round(expression, 1)), 
            vjust = -0.5, size = 3.5)

ggsave("gene_expression_barplot.png", p1, width = 10, height = 6, dpi = 300)
print(p1)

# 2. 散点图与回归线
p2 <- ggplot(gene_info, aes(x = gene_length/1000, y = expression)) +
  geom_point(aes(color = mutation_status, size = chromosome), alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("mutated" = "#E74C3C", "wild_type" = "#3498DB"),
                     labels = c("突变", "野生型")) +
  scale_size_continuous(range = c(3, 8)) +
  labs(title = "基因长度与表达水平的关系",
       x = "基因长度 (kb)", y = "表达水平",
       color = "突变状态", size = "染色体号") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  ) +
  annotate("text", x = min(gene_info$gene_length/1000), y = max(gene_info$expression),
           label = paste("r =", round(cor(gene_info$gene_length, gene_info$expression), 3)),
           hjust = 0, vjust = 1, size = 4, fontface = "italic")

ggsave("length_vs_expression_scatter.png", p2, width = 10, height = 6, dpi = 300)
print(p2)

# 3. 箱线图
p3 <- ggplot(gene_info, aes(x = mutation_status, y = expression)) +
  geom_boxplot(aes(fill = mutation_status), alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8) +
  scale_fill_manual(values = c("mutated" = "#E74C3C", "wild_type" = "#3498DB")) +
  scale_x_discrete(labels = c("突变", "野生型")) +
  labs(title = "不同突变状态的表达水平分布",
       x = "突变状态", y = "表达水平") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(size = 12)
  ) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
               fill = "yellow", color = "black")

ggsave("expression_boxplot.png", p3, width = 8, height = 6, dpi = 300)
print(p3)

# 练习3: 高级可视化
cat("练习3: 高级可视化\n")
cat("----------------\n")

# 1. 表达矩阵热图
heatmap_data <- melt(expression_matrix)
colnames(heatmap_data) <- c("Gene", "Sample", "Expression")

p4 <- ggplot(heatmap_data, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C",
                       midpoint = median(heatmap_data$Expression),
                       name = "表达水平") +
  labs(title = "基因表达热图",
       x = "样本", y = "基因") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  geom_text(aes(label = round(Expression, 1)), size = 3)

ggsave("expression_heatmap.png", p4, width = 8, height = 6, dpi = 300)
print(p4)

# 2. 多面板图
gene_info$expression_category <- ifelse(gene_info$expression > median(gene_info$expression), 
                                       "高表达", "低表达")

p5 <- ggplot(gene_info, aes(x = gene_length/1000)) +
  geom_histogram(aes(fill = mutation_status), bins = 3, alpha = 0.7, 
                 position = "identity") +
  scale_fill_manual(values = c("mutated" = "#E74C3C", "wild_type" = "#3498DB"),
                    labels = c("突变", "野生型")) +
  facet_wrap(~expression_category, scales = "free_y") +
  labs(title = "基因长度分布（按表达水平分组）",
       x = "基因长度 (kb)", y = "频数",
       fill = "突变状态") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave("length_distribution_faceted.png", p5, width = 10, height = 6, dpi = 300)
print(p5)

# 3. 小提琴图
p6 <- ggplot(gene_info, aes(x = mutation_status, y = expression)) +
  geom_violin(aes(fill = mutation_status), alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("mutated" = "#E74C3C", "wild_type" = "#3498DB")) +
  scale_x_discrete(labels = c("突变", "野生型")) +
  labs(title = "表达水平分布的小提琴图",
       x = "突变状态", y = "表达水平") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  )

ggsave("expression_violin.png", p6, width = 8, height = 6, dpi = 300)
print(p6)

# 练习4: 相关性和网络可视化
cat("练习4: 相关性可视化\n")
cat("------------------\n")

# 准备数值数据
numeric_data <- gene_info %>%
  select(chromosome, gene_length, expression) %>%
  rename(染色体 = chromosome, 基因长度 = gene_length, 表达水平 = expression)

# 计算相关性矩阵
correlation_matrix <- cor(numeric_data)

# 使用corrplot绘制相关性矩阵
png("correlation_matrix.png", width = 800, height = 800)
corrplot(correlation_matrix, 
         method = "color", 
         type = "upper",
         order = "hclust", 
         tl.cex = 1.2, 
         tl.col = "black",
         addCoef.col = "black", 
         number.cex = 1.0,
         col = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(200),
         title = "变量相关性矩阵",
         mar = c(0,0,2,0))
dev.off()

cat("相关性矩阵图已保存为 'correlation_matrix.png'\n")

# 成对散点图
png("pairwise_scatter.png", width = 1000, height = 1000)
pairs(numeric_data, 
      main = "变量间成对散点图",
      pch = 19, 
      col = "blue",
      cex = 1.5)
dev.off()

cat("成对散点图已保存为 'pairwise_scatter.png'\n")

# 练习5: 综合可视化面板
cat("练习5: 综合可视化面板\n")
cat("--------------------\n")

# 创建多个子图
p_bar <- ggplot(gene_info, aes(x = reorder(gene_name, expression), y = expression)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  coord_flip() +
  labs(title = "基因表达水平", x = "基因", y = "表达水平") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"))

p_scatter <- ggplot(gene_info, aes(x = gene_length/1000, y = expression)) +
  geom_point(aes(color = mutation_status), size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_manual(values = c("mutated" = "red", "wild_type" = "blue")) +
  labs(title = "长度 vs 表达", x = "基因长度 (kb)", y = "表达水平", color = "状态") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position = "bottom")

p_box <- ggplot(gene_info, aes(x = mutation_status, y = expression)) +
  geom_boxplot(aes(fill = mutation_status), alpha = 0.7) +
  geom_jitter(width = 0.2) +
  scale_fill_manual(values = c("mutated" = "red", "wild_type" = "blue")) +
  scale_x_discrete(labels = c("突变", "野生型")) +
  labs(title = "表达分布", x = "突变状态", y = "表达水平") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position = "none")

p_hist <- ggplot(gene_info, aes(x = expression)) +
  geom_histogram(bins = 5, fill = "orange", alpha = 0.7, color = "black") +
  geom_density(aes(y = ..density.. * nrow(gene_info) * 0.5), color = "red", size = 1) +
  labs(title = "表达分布直方图", x = "表达水平", y = "频数") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"))

# 组合图形
combined_plot <- grid.arrange(p_bar, p_scatter, p_box, p_hist, 
                             ncol = 2, nrow = 2,
                             top = "基因表达数据综合分析")

# 保存组合图形
ggsave("comprehensive_analysis.png", combined_plot, width = 12, height = 10, dpi = 300)
cat("综合分析图已保存为 'comprehensive_analysis.png'\n")

# 练习6: 交互式可视化准备
cat("练习6: 高级图形定制\n")
cat("-------------------\n")

# 创建专业的发表级图形
publication_plot <- ggplot(gene_info, aes(x = gene_length/1000, y = expression)) +
  geom_point(aes(color = mutation_status, shape = mutation_status), 
             size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", 
              fill = "gray", alpha = 0.3) +
  scale_color_manual(values = c("mutated" = "#D32F2F", "wild_type" = "#1976D2"),
                     labels = c("Mutated", "Wild-type"),
                     name = "Mutation Status") +
  scale_shape_manual(values = c("mutated" = 16, "wild_type" = 17),
                     labels = c("Mutated", "Wild-type"),
                     name = "Mutation Status") +
  labs(title = "Relationship between Gene Length and Expression Level",
       subtitle = "Analysis of 5 cancer-related genes",
       x = "Gene Length (kb)", 
       y = "Expression Level (log2 FPKM)",
       caption = "Data source: Simulated gene expression data") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray50"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = c(0.8, 0.2),
    legend.background = element_rect(fill = "white", color = "black"),
    plot.caption = element_text(size = 8, color = "gray50", hjust = 1),
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_line(color = "gray95", size = 0.25)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3)))

ggsave("publication_quality_plot.png", publication_plot, 
       width = 10, height = 8, dpi = 300)
print(publication_plot)

# 创建数据摘要表格可视化
summary_data <- gene_info %>%
  group_by(mutation_status) %>%
  summarise(
    Count = n(),
    Mean_Expression = round(mean(expression), 2),
    SD_Expression = round(sd(expression), 2),
    Mean_Length = round(mean(gene_length)/1000, 1),
    .groups = 'drop'
  ) %>%
  mutate(mutation_status = ifelse(mutation_status == "mutated", "突变", "野生型"))

# 表格可视化
table_plot <- ggplot(summary_data, aes(x = mutation_status, y = 1)) +
  geom_tile(fill = "lightblue", color = "white", size = 2) +
  geom_text(aes(label = paste("数量:", Count, "\n",
                              "平均表达:", Mean_Expression, "\n",
                              "表达标准差:", SD_Expression, "\n",
                              "平均长度:", Mean_Length, "kb")),
            size = 4, fontface = "bold") +
  labs(title = "基因数据摘要表") +
  theme_void() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 12, face = "bold")
  ) +
  scale_x_discrete(position = "bottom")

ggsave("summary_table_plot.png", table_plot, width = 10, height = 4, dpi = 300)
print(table_plot)

# 输出所有创建的图形文件列表
cat("\n创建的图形文件:\n")
plot_files <- c(
  "basic_plots_combined.png",
  "gene_expression_barplot.png", 
  "length_vs_expression_scatter.png",
  "expression_boxplot.png",
  "expression_heatmap.png",
  "length_distribution_faceted.png",
  "expression_violin.png",
  "correlation_matrix.png",
  "pairwise_scatter.png",
  "comprehensive_analysis.png",
  "publication_quality_plot.png",
  "summary_table_plot.png"
)

for (file in plot_files) {
  if (file.exists(file)) {
    cat("✓", file, "\n")
  } else {
    cat("✗", file, "(未创建)\n")
  }
}

cat("\n数据可视化练习完成！\n")
cat("======================\n")