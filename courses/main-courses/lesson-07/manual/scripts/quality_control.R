# 单细胞数据质量控制专用脚本
# 课程：高通量测序数据分析 - 第7次课
# 作者：王运生
# 日期：2025-01-21
# 用法：source("scripts/quality_control.R")

cat("=== 单细胞数据质量控制分析 ===\n")

# 加载必要的包
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# 质量控制函数
perform_qc_analysis <- function(seurat_obj, 
                               min_features = 200,
                               max_features = 2500, 
                               max_mt_percent = 20,
                               max_rb_percent = 50,
                               plot_prefix = "qc") {
  
  cat("开始质量控制分析...\n")
  
  # 计算质控指标
  cat("计算质控指标...\n")
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  
  # 计算复杂度指标（基因数/UMI数比值）
  seurat_obj[["log10GenesPerUMI"]] <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
  
  # 基本统计
  cat("\n=== 质控指标统计 ===\n")
  
  qc_stats <- data.frame(
    Metric = c("细胞数量", "基因数量", "检测基因数(中位数)", "UMI数量(中位数)", 
               "线粒体基因%(中位数)", "核糖体基因%(中位数)"),
    Value = c(
      ncol(seurat_obj),
      nrow(seurat_obj),
      median(seurat_obj$nFeature_RNA),
      median(seurat_obj$nCount_RNA),
      round(median(seurat_obj$percent.mt), 2),
      round(median(seurat_obj$percent.rb), 2)
    )
  )
  
  print(qc_stats)
  
  # 详细统计
  cat("\n=== 详细统计信息 ===\n")
  cat("检测基因数 - 范围:", range(seurat_obj$nFeature_RNA), "\n")
  cat("UMI数量 - 范围:", range(seurat_obj$nCount_RNA), "\n")
  cat("线粒体基因% - 范围:", round(range(seurat_obj$percent.mt), 2), "\n")
  cat("核糖体基因% - 范围:", round(range(seurat_obj$percent.rb), 2), "\n")
  
  # 创建质控可视化
  cat("\n创建质控图表...\n")
  
  # 1. 小提琴图
  p1 <- VlnPlot(seurat_obj, 
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), 
                ncol = 2, pt.size = 0.1) +
    plot_annotation(title = "质量控制指标分布")
  
  ggsave(paste0("plots/", plot_prefix, "_violin_plots.png"), p1, 
         width = 12, height = 8, dpi = 300)
  
  # 2. 散点图矩阵
  p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    geom_hline(yintercept = max_mt_percent, color = "red", linetype = "dashed") +
    ggtitle("UMI数量 vs 线粒体基因%")
  
  p3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_hline(yintercept = c(min_features, max_features), color = "red", linetype = "dashed") +
    ggtitle("UMI数量 vs 检测基因数")
  
  p4 <- FeatureScatter(seurat_obj, feature1 = "percent.mt", feature2 = "percent.rb") +
    geom_vline(xintercept = max_mt_percent, color = "red", linetype = "dashed") +
    ggtitle("线粒体基因% vs 核糖体基因%")
  
  scatter_plots <- (p2 + p3) / p4
  ggsave(paste0("plots/", plot_prefix, "_scatter_plots.png"), scatter_plots, 
         width = 12, height = 8, dpi = 300)
  
  # 3. 直方图
  p5 <- ggplot(seurat_obj@meta.data, aes(x = nFeature_RNA)) +
    geom_histogram(bins = 50, fill = "skyblue", alpha = 0.7) +
    geom_vline(xintercept = c(min_features, max_features), color = "red", linetype = "dashed") +
    labs(title = "检测基因数分布", x = "检测基因数", y = "细胞数量") +
    theme_bw()
  
  p6 <- ggplot(seurat_obj@meta.data, aes(x = percent.mt)) +
    geom_histogram(bins = 50, fill = "lightcoral", alpha = 0.7) +
    geom_vline(xintercept = max_mt_percent, color = "red", linetype = "dashed") +
    labs(title = "线粒体基因%分布", x = "线粒体基因%", y = "细胞数量") +
    theme_bw()
  
  hist_plots <- p5 + p6
  ggsave(paste0("plots/", plot_prefix, "_histograms.png"), hist_plots, 
         width = 12, height = 4, dpi = 300)
  
  # 4. 复杂度分析
  p7 <- ggplot(seurat_obj@meta.data, aes(x = log10GenesPerUMI)) +
    geom_histogram(bins = 50, fill = "lightgreen", alpha = 0.7) +
    geom_vline(xintercept = 0.8, color = "red", linetype = "dashed") +
    labs(title = "细胞复杂度分布", 
         x = "log10(基因数/UMI数)", 
         y = "细胞数量",
         subtitle = "低复杂度可能表示细胞质量差或双细胞") +
    theme_bw()
  
  ggsave(paste0("plots/", plot_prefix, "_complexity.png"), p7, 
         width = 8, height = 6, dpi = 300)
  
  # 识别异常细胞
  cat("\n=== 异常细胞识别 ===\n")
  
  # 计算异常细胞数量
  low_feature_cells <- sum(seurat_obj$nFeature_RNA < min_features)
  high_feature_cells <- sum(seurat_obj$nFeature_RNA > max_features)
  high_mt_cells <- sum(seurat_obj$percent.mt > max_mt_percent)
  high_rb_cells <- sum(seurat_obj$percent.rb > max_rb_percent)
  low_complexity_cells <- sum(seurat_obj$log10GenesPerUMI < 0.8)
  
  anomaly_stats <- data.frame(
    异常类型 = c("低基因数细胞", "高基因数细胞", "高线粒体%细胞", 
                "高核糖体%细胞", "低复杂度细胞"),
    细胞数量 = c(low_feature_cells, high_feature_cells, high_mt_cells, 
                high_rb_cells, low_complexity_cells),
    百分比 = round(c(low_feature_cells, high_feature_cells, high_mt_cells, 
                    high_rb_cells, low_complexity_cells) / ncol(seurat_obj) * 100, 2)
  )
  
  print(anomaly_stats)
  
  # 应用过滤条件
  cat("\n应用过滤条件...\n")
  cat("过滤前细胞数量:", ncol(seurat_obj), "\n")
  
  # 标记要保留的细胞
  keep_cells <- seurat_obj$nFeature_RNA >= min_features & 
                seurat_obj$nFeature_RNA <= max_features & 
                seurat_obj$percent.mt <= max_mt_percent
  
  filtered_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[keep_cells])
  
  cat("过滤后细胞数量:", ncol(filtered_obj), "\n")
  cat("过滤掉的细胞数量:", ncol(seurat_obj) - ncol(filtered_obj), "\n")
  cat("过滤比例:", round((ncol(seurat_obj) - ncol(filtered_obj)) / ncol(seurat_obj) * 100, 2), "%\n")
  
  # 创建过滤前后对比图
  before_after_data <- data.frame(
    Stage = rep(c("过滤前", "过滤后"), each = ncol(seurat_obj)),
    nFeature_RNA = c(seurat_obj$nFeature_RNA, 
                     ifelse(keep_cells, seurat_obj$nFeature_RNA, NA)),
    percent.mt = c(seurat_obj$percent.mt, 
                   ifelse(keep_cells, seurat_obj$percent.mt, NA)),
    Status = c(rep("所有细胞", ncol(seurat_obj)),
               ifelse(keep_cells, "保留", "过滤"))
  )
  
  before_after_data <- before_after_data[!is.na(before_after_data$nFeature_RNA), ]
  
  p8 <- ggplot(before_after_data, aes(x = Stage, y = nFeature_RNA, fill = Status)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, alpha = 0.8) +
    labs(title = "过滤前后基因数分布对比", y = "检测基因数") +
    theme_bw()
  
  p9 <- ggplot(before_after_data, aes(x = Stage, y = percent.mt, fill = Status)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, alpha = 0.8) +
    labs(title = "过滤前后线粒体基因%对比", y = "线粒体基因%") +
    theme_bw()
  
  comparison_plots <- p8 + p9
  ggsave(paste0("plots/", plot_prefix, "_before_after_filtering.png"), comparison_plots, 
         width = 12, height = 6, dpi = 300)
  
  # 保存质控报告
  qc_report <- list(
    analysis_date = Sys.Date(),
    original_cells = ncol(seurat_obj),
    filtered_cells = ncol(filtered_obj),
    cells_removed = ncol(seurat_obj) - ncol(filtered_obj),
    filter_rate = round((ncol(seurat_obj) - ncol(filtered_obj)) / ncol(seurat_obj) * 100, 2),
    parameters = list(
      min_features = min_features,
      max_features = max_features,
      max_mt_percent = max_mt_percent,
      max_rb_percent = max_rb_percent
    ),
    statistics = qc_stats,
    anomalies = anomaly_stats
  )
  
  saveRDS(qc_report, paste0("results/", plot_prefix, "_report.rds"))
  write.csv(qc_stats, paste0("results/", plot_prefix, "_statistics.csv"), row.names = FALSE)
  write.csv(anomaly_stats, paste0("results/", plot_prefix, "_anomalies.csv"), row.names = FALSE)
  
  cat("\n质控分析完成！\n")
  cat("报告文件已保存到 results/ 目录\n")
  cat("图表文件已保存到 plots/ 目录\n")
  
  return(filtered_obj)
}

# 基因过滤函数
filter_genes <- function(seurat_obj, min_cells = 3) {
  cat("\n=== 基因过滤分析 ===\n")
  
  genes_before <- nrow(seurat_obj)
  
  # 计算每个基因在多少个细胞中表达
  gene_expression_freq <- Matrix::rowSums(seurat_obj@assays$RNA@counts > 0)
  
  cat("基因表达频率统计:\n")
  cat("表达细胞数 < 3:", sum(gene_expression_freq < 3), "个基因\n")
  cat("表达细胞数 >= 3:", sum(gene_expression_freq >= 3), "个基因\n")
  cat("表达细胞数 >= 10:", sum(gene_expression_freq >= 10), "个基因\n")
  
  # 过滤低表达基因
  keep_genes <- gene_expression_freq >= min_cells
  filtered_obj <- subset(seurat_obj, features = rownames(seurat_obj)[keep_genes])
  
  genes_after <- nrow(filtered_obj)
  
  cat("\n基因过滤结果:\n")
  cat("过滤前基因数量:", genes_before, "\n")
  cat("过滤后基因数量:", genes_after, "\n")
  cat("过滤掉的基因数量:", genes_before - genes_after, "\n")
  cat("过滤比例:", round((genes_before - genes_after) / genes_before * 100, 2), "%\n")
  
  # 可视化基因表达频率
  freq_data <- data.frame(
    gene_id = names(gene_expression_freq),
    expression_frequency = gene_expression_freq,
    keep = keep_genes
  )
  
  p_gene_freq <- ggplot(freq_data, aes(x = expression_frequency, fill = keep)) +
    geom_histogram(bins = 50, alpha = 0.7) +
    geom_vline(xintercept = min_cells, color = "red", linetype = "dashed") +
    scale_x_log10() +
    labs(title = "基因表达频率分布", 
         x = "表达该基因的细胞数量 (log10)", 
         y = "基因数量",
         fill = "是否保留") +
    theme_bw()
  
  ggsave("plots/gene_expression_frequency.png", p_gene_freq, 
         width = 10, height = 6, dpi = 300)
  
  return(filtered_obj)
}

# 双细胞检测函数（简化版）
detect_doublets <- function(seurat_obj, expected_doublet_rate = 0.08) {
  cat("\n=== 双细胞检测（简化版）===\n")
  
  # 基于基因数和UMI数的简单双细胞检测
  # 真正的双细胞检测需要使用专门的工具如DoubletFinder
  
  # 计算阈值（基于分位数）
  feature_threshold <- quantile(seurat_obj$nFeature_RNA, 0.95)
  count_threshold <- quantile(seurat_obj$nCount_RNA, 0.95)
  
  # 识别可能的双细胞
  potential_doublets <- seurat_obj$nFeature_RNA > feature_threshold & 
                       seurat_obj$nCount_RNA > count_threshold
  
  n_doublets <- sum(potential_doublets)
  doublet_rate <- n_doublets / ncol(seurat_obj)
  
  cat("潜在双细胞数量:", n_doublets, "\n")
  cat("双细胞比例:", round(doublet_rate * 100, 2), "%\n")
  cat("预期双细胞比例:", round(expected_doublet_rate * 100, 2), "%\n")
  
  if (doublet_rate > expected_doublet_rate * 1.5) {
    cat("警告：双细胞比例可能偏高，建议使用DoubletFinder进行详细检测\n")
  }
  
  # 可视化潜在双细胞
  doublet_data <- data.frame(
    nFeature_RNA = seurat_obj$nFeature_RNA,
    nCount_RNA = seurat_obj$nCount_RNA,
    potential_doublet = potential_doublets
  )
  
  p_doublet <- ggplot(doublet_data, aes(x = nCount_RNA, y = nFeature_RNA, 
                                       color = potential_doublet)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = feature_threshold, color = "red", linetype = "dashed") +
    geom_vline(xintercept = count_threshold, color = "red", linetype = "dashed") +
    labs(title = "潜在双细胞检测", 
         x = "UMI数量", 
         y = "检测基因数",
         color = "潜在双细胞") +
    theme_bw()
  
  ggsave("plots/potential_doublets.png", p_doublet, 
         width = 10, height = 8, dpi = 300)
  
  return(potential_doublets)
}

# 主函数：如果直接运行此脚本
if (sys.nframe() == 0) {
  cat("质量控制脚本已加载\n")
  cat("使用方法:\n")
  cat("1. filtered_obj <- perform_qc_analysis(seurat_obj)\n")
  cat("2. filtered_obj <- filter_genes(filtered_obj)\n")
  cat("3. doublets <- detect_doublets(filtered_obj)\n")
}

cat("质量控制函数加载完成\n")