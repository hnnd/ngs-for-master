#!/usr/bin/env Rscript

# DESeq2差异表达分析脚本
# 课程：高通量测序数据分析 - 第5次课
# 作者：王运生
# 邮箱：wangys@hunau.edu.cn

# 加载必要的包
suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(dplyr)
    library(tibble)
})

# 设置参数
WORK_DIR <- "~/ngs-analysis/lesson-05"
PADJ_CUTOFF <- 0.05
FC_CUTOFF <- 1.0  # log2 fold change cutoff

# 日志函数
log_info <- function(msg) {
    cat("[INFO]", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", msg, "\n")
}

log_error <- function(msg) {
    cat("[ERROR]", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", msg, "\n")
}

# 设置工作目录
setwd(WORK_DIR)

# 创建输出目录
dir.create("results/deseq2", showWarnings = FALSE, recursive = TRUE)
dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)

log_info("开始DESeq2差异表达分析")

# 1. 读取数据
log_info("读取基因计数矩阵...")

count_file <- "results/quantification/gene_counts.txt"
if (!file.exists(count_file)) {
    log_error(paste("计数文件不存在:", count_file))
    quit(status = 1)
}

# 读取featureCounts输出文件
countData <- read.table(count_file, header = TRUE, row.names = 1, skip = 1)

# 只保留计数列（去除注释列：Chr, Start, End, Strand, Length）
countData <- countData[, 6:ncol(countData)]

# 修改列名（去除路径和扩展名）
colnames(countData) <- gsub("results.alignment.|_sorted.bam", "", colnames(countData))
colnames(countData) <- gsub("\\.", "_", colnames(countData))

log_info(paste("读取到", nrow(countData), "个基因，", ncol(countData), "个样本"))
log_info("样本名称:")
print(colnames(countData))

# 2. 创建样本信息表
log_info("创建样本信息表...")

# 根据样本名称自动推断分组
sample_names <- colnames(countData)
conditions <- ifelse(grepl("^ctrl", sample_names), "control", "treatment")

sampleInfo <- data.frame(
    row.names = sample_names,
    condition = factor(conditions, levels = c("control", "treatment")),
    batch = factor(rep(1:3, length.out = length(sample_names)))
)

log_info("样本信息:")
print(sampleInfo)

# 检查样本名称一致性
if (!all(rownames(sampleInfo) == colnames(countData))) {
    log_error("样本名称不匹配")
    quit(status = 1)
}

# 3. 创建DESeq2对象
log_info("创建DESeq2对象...")

dds <- DESeqDataSetFromMatrix(
    countData = countData,
    colData = sampleInfo,
    design = ~ condition
)

# 过滤低表达基因
log_info("过滤低表达基因...")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

log_info(paste("过滤后保留", nrow(dds), "个基因"))

# 4. 运行差异分析
log_info("运行DESeq2差异分析...")

dds <- DESeq(dds)

# 获取结果
res <- results(dds, contrast = c("condition", "treatment", "control"))

# 查看结果摘要
log_info("差异分析结果摘要:")
print(summary(res))

# 5. 保存结果
log_info("保存分析结果...")

# 保存完整结果
res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    arrange(padj)

write.csv(res_df, "results/deseq2/deseq2_results.csv", row.names = FALSE)

# 获取显著差异基因
sig_genes <- res_df %>%
    filter(padj < PADJ_CUTOFF & abs(log2FoldChange) > FC_CUTOFF)

write.csv(sig_genes, "results/deseq2/significant_genes.csv", row.names = FALSE)

log_info(paste("显著差异基因数量:", nrow(sig_genes)))
log_info(paste("上调基因:", sum(sig_genes$log2FoldChange > 0)))
log_info(paste("下调基因:", sum(sig_genes$log2FoldChange < 0)))

# 6. 数据转换用于可视化
log_info("进行数据转换...")

# 方差稳定化转换
vsd <- vst(dds, blind = FALSE)

# rlog转换（如果样本数较少）
if (ncol(countData) <= 30) {
    rld <- rlog(dds, blind = FALSE)
}

# 7. 生成可视化图表
log_info("生成可视化图表...")

# MA图
log_info("生成MA图...")
pdf("results/plots/MA_plot.pdf", width = 8, height = 6)
plotMA(res, main = "MA Plot", ylim = c(-5, 5))
abline(h = c(-FC_CUTOFF, FC_CUTOFF), col = "red", lty = 2)
dev.off()

# 火山图
log_info("生成火山图...")
volcano_data <- res_df %>%
    mutate(
        significance = case_when(
            padj < PADJ_CUTOFF & log2FoldChange > FC_CUTOFF ~ "Up",
            padj < PADJ_CUTOFF & log2FoldChange < -FC_CUTOFF ~ "Down",
            TRUE ~ "NS"
        )
    )

p_volcano <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = significance), alpha = 0.6) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
    geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    labs(
        title = "Volcano Plot",
        x = "log2 Fold Change",
        y = "-log10(p-value)",
        color = "Significance"
    ) +
    theme_minimal() +
    theme(legend.position = "top")

ggsave("results/plots/volcano_plot.pdf", p_volcano, width = 10, height = 8)

# PCA图
log_info("生成PCA图...")
pdf("results/plots/PCA_plot.pdf", width = 8, height = 6)
plotPCA(vsd, intgroup = "condition") +
    ggtitle("PCA Plot") +
    theme_minimal()
dev.off()

# 样本距离热图
log_info("生成样本距离热图...")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, rownames(sampleInfo), sep = "-")
colnames(sampleDistMatrix) <- NULL

pdf("results/plots/sample_distance_heatmap.pdf", width = 8, height = 6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample Distance Heatmap")
dev.off()

# 差异基因热图
if (nrow(sig_genes) > 0) {
    log_info("生成差异基因热图...")
    
    # 选择前50个最显著的差异基因
    top_genes <- head(sig_genes$gene_id, 50)
    
    if (length(top_genes) > 0) {
        pdf("results/plots/heatmap_top50.pdf", width = 10, height = 12)
        pheatmap(assay(vsd)[top_genes, ],
                 cluster_rows = TRUE,
                 show_rownames = TRUE,
                 cluster_cols = TRUE,
                 annotation_col = sampleInfo,
                 scale = "row",
                 main = "Top 50 Differentially Expressed Genes",
                 fontsize_row = 8)
        dev.off()
    }
}

# 基因表达分布图
log_info("生成基因表达分布图...")
pdf("results/plots/expression_distribution.pdf", width = 12, height = 8)
par(mfrow = c(2, 2))

# 原始计数分布
plotCounts(dds, gene = which.min(res$padj), intgroup = "condition", main = "Most significant gene")

# 标准化计数分布
counts_norm <- counts(dds, normalized = TRUE)
hist(log10(counts_norm + 1), breaks = 50, main = "Normalized counts distribution", 
     xlab = "log10(normalized counts + 1)")

# 离散度估计图
plotDispEsts(dds, main = "Dispersion estimates")

# p值分布
hist(res$pvalue, breaks = 50, main = "P-value distribution", xlab = "p-value")

dev.off()

# 8. 生成分析摘要
log_info("生成分析摘要...")

summary_stats <- list(
    "Total genes analyzed" = nrow(res_df),
    "Significant genes (padj < 0.05, |FC| > 1)" = nrow(sig_genes),
    "Upregulated genes" = sum(sig_genes$log2FoldChange > 0, na.rm = TRUE),
    "Downregulated genes" = sum(sig_genes$log2FoldChange < 0, na.rm = TRUE),
    "Mean |log2FC| of significant genes" = round(mean(abs(sig_genes$log2FoldChange), na.rm = TRUE), 2),
    "Max upregulation (log2FC)" = round(max(sig_genes$log2FoldChange, na.rm = TRUE), 2),
    "Max downregulation (log2FC)" = round(min(sig_genes$log2FoldChange, na.rm = TRUE), 2)
)

# 保存摘要到文件
summary_file <- "results/deseq2/analysis_summary.txt"
cat("RNA-seq差异表达分析摘要\n", file = summary_file)
cat("========================\n\n", file = summary_file, append = TRUE)
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", file = summary_file, append = TRUE)

for (i in 1:length(summary_stats)) {
    cat(names(summary_stats)[i], ":", summary_stats[[i]], "\n", 
        file = summary_file, append = TRUE)
}

cat("\n参数设置:\n", file = summary_file, append = TRUE)
cat("- padj cutoff:", PADJ_CUTOFF, "\n", file = summary_file, append = TRUE)
cat("- log2FC cutoff:", FC_CUTOFF, "\n", file = summary_file, append = TRUE)

# 打印摘要到控制台
log_info("=== 分析摘要 ===")
for (i in 1:length(summary_stats)) {
    log_info(paste(names(summary_stats)[i], ":", summary_stats[[i]]))
}

# 9. 保存R对象
log_info("保存R对象...")
save(dds, res, vsd, sig_genes, file = "results/deseq2/deseq2_objects.RData")

log_info("DESeq2差异表达分析完成！")
log_info("结果文件保存在 results/deseq2/ 目录")
log_info("图表文件保存在 results/plots/ 目录")

# 输出文件列表
log_info("生成的文件:")
output_files <- c(
    "results/deseq2/deseq2_results.csv",
    "results/deseq2/significant_genes.csv",
    "results/deseq2/analysis_summary.txt",
    "results/deseq2/deseq2_objects.RData",
    "results/plots/MA_plot.pdf",
    "results/plots/volcano_plot.pdf",
    "results/plots/PCA_plot.pdf",
    "results/plots/sample_distance_heatmap.pdf",
    "results/plots/expression_distribution.pdf"
)

if (nrow(sig_genes) > 0) {
    output_files <- c(output_files, "results/plots/heatmap_top50.pdf")
}

for (file in output_files) {
    if (file.exists(file)) {
        log_info(paste("✓", file))
    } else {
        log_info(paste("✗", file, "(未生成)"))
    }
}