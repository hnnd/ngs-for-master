#!/usr/bin/env Rscript
# 多组学数据整合分析脚本
# 作者：王运生
# 日期：2025-01-21
# 用法：Rscript data_integration.R

cat("=== 多组学数据整合分析 ===\n")

# 加载必要的包
suppressMessages({
    library(mixOmics)
    library(ggplot2)
    library(pheatmap)
    library(VennDiagram)
    library(RColorBrewer)
})

# 设置工作目录
if (!dir.exists("~/ngs-analysis/lesson-08")) {
    stop("请先运行setup.sh设置环境")
}
setwd("~/ngs-analysis/lesson-08")

# 创建输出目录
dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

# 1. 数据读取和预处理
cat("1. 读取数据...\n")
genomics <- read.csv("data/genomics_data.csv", row.names = 1)
transcriptomics <- read.csv("data/transcriptomics_data.csv", row.names = 1)
proteomics <- read.csv("data/proteomics_data.csv", row.names = 1)
metadata <- read.csv("data/metadata.csv", row.names = 1)

cat("数据维度:\n")
cat("  基因组:", dim(genomics), "\n")
cat("  转录组:", dim(transcriptomics), "\n")
cat("  蛋白质组:", dim(proteomics), "\n")
cat("  样本信息:", dim(metadata), "\n")

# 2. 数据质量评估
cat("\n2. 数据质量评估...\n")
cat("缺失值统计:\n")
cat("  基因组:", sum(is.na(genomics)), "\n")
cat("  转录组:", sum(is.na(transcriptomics)), "\n")
cat("  蛋白质组:", sum(is.na(proteomics)), "\n")

# 3. 数据预处理
cat("\n3. 数据预处理...\n")

# 转录组数据log2转换和标准化
transcriptomics_log <- log2(transcriptomics + 1)
transcriptomics_scaled <- scale(transcriptomics_log)

# 蛋白质组数据标准化
proteomics_scaled <- scale(proteomics)

# 基因组数据（已经是0/1/2编码，不需要标准化）
genomics_processed <- genomics

# 处理缺失值（使用均值填补）
transcriptomics_scaled[is.na(transcriptomics_scaled)] <- 0
proteomics_scaled[is.na(proteomics_scaled)] <- 0

# 保存预处理后的数据
write.csv(genomics_processed, "results/genomics_processed.csv")
write.csv(transcriptomics_scaled, "results/transcriptomics_scaled.csv")
write.csv(proteomics_scaled, "results/proteomics_scaled.csv")

cat("数据预处理完成\n")

# 4. 数据分布可视化
cat("\n4. 生成数据分布图...\n")

# 转录组数据分布（预处理前后对比）
pdf("plots/data_distribution.pdf", width = 12, height = 8)
par(mfrow = c(2, 3))

# 原始转录组数据分布
hist(as.matrix(transcriptomics), main = "转录组原始数据分布", 
     xlab = "表达量", col = "lightblue", breaks = 50)

# Log2转换后分布
hist(as.matrix(transcriptomics_log), main = "转录组Log2转换后分布", 
     xlab = "Log2(表达量+1)", col = "lightgreen", breaks = 50)

# 标准化后分布
hist(as.matrix(transcriptomics_scaled), main = "转录组标准化后分布", 
     xlab = "标准化值", col = "lightcoral", breaks = 50)

# 蛋白质组数据分布
hist(as.matrix(proteomics), main = "蛋白质组原始数据分布", 
     xlab = "丰度", col = "lightyellow", breaks = 50)

hist(as.matrix(proteomics_scaled), main = "蛋白质组标准化后分布", 
     xlab = "标准化值", col = "lightpink", breaks = 50)

# 基因组数据分布
barplot(table(as.matrix(genomics)), main = "基因组数据分布", 
        xlab = "基因型", ylab = "频数", col = "lightgray")

dev.off()

# 5. 样本相关性分析
cat("\n5. 样本相关性分析...\n")

# 计算样本间相关性
cor_genomics <- cor(t(genomics_processed))
cor_transcriptomics <- cor(t(transcriptomics_scaled))
cor_proteomics <- cor(t(proteomics_scaled))

# 绘制相关性热图
pdf("plots/sample_correlation_heatmaps.pdf", width = 15, height = 5)
par(mfrow = c(1, 3))

# 基因组相关性
pheatmap(cor_genomics, 
         main = "基因组样本相关性",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = FALSE, show_colnames = FALSE,
         annotation_row = metadata[, "group", drop = FALSE],
         filename = NA)

# 转录组相关性
pheatmap(cor_transcriptomics, 
         main = "转录组样本相关性",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = FALSE, show_colnames = FALSE,
         annotation_row = metadata[, "group", drop = FALSE],
         filename = NA)

# 蛋白质组相关性
pheatmap(cor_proteomics, 
         main = "蛋白质组样本相关性",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = FALSE, show_colnames = FALSE,
         annotation_row = metadata[, "group", drop = FALSE],
         filename = NA)

dev.off()

# 6. PCA分析
cat("\n6. 主成分分析...\n")

# 对每个组学数据进行PCA
pca_genomics <- prcomp(genomics_processed, scale. = TRUE)
pca_transcriptomics <- prcomp(transcriptomics_scaled)
pca_proteomics <- prcomp(proteomics_scaled)

# 绘制PCA图
pdf("plots/pca_analysis.pdf", width = 15, height = 5)
par(mfrow = c(1, 3))

# 基因组PCA
plot(pca_genomics$x[,1], pca_genomics$x[,2], 
     col = as.factor(metadata$group), pch = 19,
     xlab = paste0("PC1 (", round(summary(pca_genomics)$importance[2,1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca_genomics)$importance[2,2]*100, 1), "%)"),
     main = "基因组PCA")
legend("topright", levels(as.factor(metadata$group)), 
       col = 1:length(levels(as.factor(metadata$group))), pch = 19)

# 转录组PCA
plot(pca_transcriptomics$x[,1], pca_transcriptomics$x[,2], 
     col = as.factor(metadata$group), pch = 19,
     xlab = paste0("PC1 (", round(summary(pca_transcriptomics)$importance[2,1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca_transcriptomics)$importance[2,2]*100, 1), "%)"),
     main = "转录组PCA")
legend("topright", levels(as.factor(metadata$group)), 
       col = 1:length(levels(as.factor(metadata$group))), pch = 19)

# 蛋白质组PCA
plot(pca_proteomics$x[,1], pca_proteomics$x[,2], 
     col = as.factor(metadata$group), pch = 19,
     xlab = paste0("PC1 (", round(summary(pca_proteomics)$importance[2,1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca_proteomics)$importance[2,2]*100, 1), "%)"),
     main = "蛋白质组PCA")
legend("topright", levels(as.factor(metadata$group)), 
       col = 1:length(levels(as.factor(metadata$group))), pch = 19)

dev.off()

# 7. mixOmics多组学整合分析
cat("\n7. mixOmics多组学整合分析...\n")

# 准备数据列表
X <- list(genomics = t(genomics_processed),
          transcriptomics = t(transcriptomics_scaled),
          proteomics = t(proteomics_scaled))

# 准备分组信息
Y <- as.factor(metadata$group)

# 执行多组学PLS-DA分析
cat("执行PLS-DA分析...\n")
result.diablo <- block.plsda(X, Y, ncomp = 2)

# 查看结果摘要
cat("PLS-DA分析完成\n")
print(result.diablo)

# 8. 可视化PLS-DA结果
cat("\n8. 生成PLS-DA结果图...\n")

# 样本得分图
pdf("plots/multiomics_plsda_samples.pdf", width = 10, height = 8)
plotIndiv(result.diablo, ind.names = FALSE, legend = TRUE, 
          title = "多组学PLS-DA样本得分图",
          ellipse = TRUE, ellipse.level = 0.95)
dev.off()

# 变量载荷图
pdf("plots/multiomics_plsda_variables.pdf", width = 12, height = 8)
plotVar(result.diablo, cutoff = 0.5, title = "多组学PLS-DA变量载荷图")
dev.off()

# 9. 网络分析
cat("\n9. 构建多组学网络...\n")

# 构建相关性网络
pdf("plots/multiomics_network.pdf", width = 12, height = 10)
network(result.diablo, blocks = c(1,2,3), 
        color.node = c("darkorchid", "brown1", "lightgreen"),
        cutoff = 0.4,
        save = 'pdf', name.save = 'multiomics_network')
dev.off()

# 10. 提取重要特征
cat("\n10. 提取重要特征...\n")

# 提取每个组学的重要特征
important.features <- selectVar(result.diablo, comp = 1)

# 保存重要特征
sink("results/important_features.txt")
cat("=== 重要特征列表 ===\n\n")

cat("基因组重要特征 (Top 10):\n")
if (length(important.features$genomics$name) > 0) {
    print(head(important.features$genomics$name, 10))
} else {
    cat("无显著特征\n")
}

cat("\n转录组重要特征 (Top 10):\n")
if (length(important.features$transcriptomics$name) > 0) {
    print(head(important.features$transcriptomics$name, 10))
} else {
    cat("无显著特征\n")
}

cat("\n蛋白质组重要特征 (Top 10):\n")
if (length(important.features$proteomics$name) > 0) {
    print(head(important.features$proteomics$name, 10))
} else {
    cat("无显著特征\n")
}
sink()

# 11. 生成综合分析报告
cat("\n11. 生成分析报告...\n")

sink("results/integration_analysis_report.txt")
cat("=== 多组学数据整合分析报告 ===\n")
cat("分析时间:", as.character(Sys.time()), "\n\n")

cat("1. 数据概况:\n")
cat("   样本数:", nrow(metadata), "\n")
cat("   基因组特征数:", ncol(genomics_processed), "\n")
cat("   转录组特征数:", ncol(transcriptomics_scaled), "\n")
cat("   蛋白质组特征数:", ncol(proteomics_scaled), "\n")
cat("   分组:", paste(levels(Y), collapse = ", "), "\n\n")

cat("2. 数据预处理:\n")
cat("   转录组数据: Log2转换 + 标准化\n")
cat("   蛋白质组数据: 标准化\n")
cat("   基因组数据: 无需处理\n")
cat("   缺失值: 均值填补\n\n")

cat("3. PLS-DA分析结果:\n")
cat("   成分数:", result.diablo$ncomp, "\n")
cat("   数据块数:", length(result.diablo$names$blocks), "\n")
cat("   数据块:", paste(result.diablo$names$blocks, collapse = ", "), "\n\n")

cat("4. 重要特征数量:\n")
cat("   基因组:", length(important.features$genomics$name), "\n")
cat("   转录组:", length(important.features$transcriptomics$name), "\n")
cat("   蛋白质组:", length(important.features$proteomics$name), "\n\n")

cat("5. 输出文件:\n")
cat("   预处理数据: results/\n")
cat("   可视化图表: plots/\n")
cat("   重要特征: results/important_features.txt\n")
cat("   分析报告: results/integration_analysis_report.txt\n\n")

cat("分析完成！\n")
sink()

# 12. 清理和总结
cat("\n=== 分析完成 ===\n")
cat("输出文件:\n")
cat("  - 预处理数据: results/\n")
cat("  - 可视化图表: plots/\n")
cat("  - 重要特征: results/important_features.txt\n")
cat("  - 分析报告: results/integration_analysis_report.txt\n")

cat("\n多组学数据整合分析完成！\n")
cat("请查看生成的图表和报告文件。\n")