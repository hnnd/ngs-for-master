# 单细胞轨迹分析脚本（扩展练习）
# 课程：高通量测序数据分析 - 第7次课
# 作者：王运生
# 日期：2025-01-21
# 用法：source("scripts/trajectory_analysis.R")

cat("=== 单细胞轨迹分析（扩展练习）===\n")

# 检查是否已完成基本分析
if (!exists("pbmc") || !file.exists("results/pbmc_analysis_complete.rds")) {
  cat("请先完成基本的单细胞分析\n")
  cat("加载已保存的分析结果...\n")
  
  if (file.exists("results/pbmc_analysis_complete.rds")) {
    pbmc <- readRDS("results/pbmc_analysis_complete.rds")
    cat("成功加载Seurat对象\n")
  } else {
    stop("未找到分析结果文件，请先运行complete_analysis.R")
  }
}

# 检查并安装Monocle3
if (!requireNamespace("monocle3", quietly = TRUE)) {
  cat("安装Monocle3包...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("monocle3")
}

# 加载必要的包
suppressPackageStartupMessages({
  library(monocle3)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

cat("已加载轨迹分析相关包\n")

# 步骤1: 准备轨迹分析数据
cat("\n步骤1: 准备轨迹分析数据\n")

# 选择感兴趣的细胞类型进行轨迹分析
# 这里我们选择T细胞亚群，因为它们可能存在分化轨迹
t_cell_types <- c("Naive CD4 T", "Memory CD4 T", "CD8 T")
pbmc_t <- subset(pbmc, idents = t_cell_types)

cat("选择的细胞类型:", paste(t_cell_types, collapse = ", "), "\n")
cat("用于轨迹分析的细胞数量:", ncol(pbmc_t), "\n")

# 步骤2: 转换为Monocle对象
cat("\n步骤2: 转换为Monocle对象\n")

# 提取表达数据
expression_matrix <- GetAssayData(pbmc_t, assay = "RNA", slot = "counts")

# 提取细胞元数据
cell_metadata <- pbmc_t@meta.data
cell_metadata$cell_type <- Idents(pbmc_t)

# 提取基因元数据
gene_metadata <- data.frame(
  gene_id = rownames(expression_matrix),
  gene_short_name = rownames(expression_matrix)
)
rownames(gene_metadata) <- rownames(expression_matrix)

# 创建Monocle对象
cds <- new_cell_data_set(expression_matrix,
                        cell_metadata = cell_metadata,
                        gene_metadata = gene_metadata)

cat("Monocle对象创建完成\n")
print(cds)

# 步骤3: 预处理数据
cat("\n步骤3: 预处理数据\n")

# 预处理
cds <- preprocess_cds(cds, num_dim = 50)

# 可视化预处理结果
p1 <- plot_pc_variance_explained(cds)
ggsave("plots/15_monocle_pca_variance.png", p1, width = 8, height = 6, dpi = 300)

# 步骤4: 降维分析
cat("\n步骤4: 降维分析\n")

# 使用UMAP降维
cds <- reduce_dimension(cds, reduction_method = "UMAP")

# 可视化降维结果
p2 <- plot_cells(cds, color_cells_by = "cell_type", label_cell_groups = TRUE)
ggsave("plots/16_monocle_umap_celltypes.png", p2, width = 10, height = 8, dpi = 300)

# 步骤5: 聚类分析
cat("\n步骤5: 聚类分析\n")

# 聚类
cds <- cluster_cells(cds)

# 可视化聚类结果
p3 <- plot_cells(cds, color_cells_by = "cluster", label_cell_groups = TRUE)
ggsave("plots/17_monocle_clusters.png", p3, width = 10, height = 8, dpi = 300)

# 比较Monocle聚类和Seurat细胞类型
p4 <- plot_cells(cds, color_cells_by = "cell_type") + 
      ggtitle("Seurat细胞类型")
p5 <- plot_cells(cds, color_cells_by = "cluster") + 
      ggtitle("Monocle聚类")

ggsave("plots/18_clustering_comparison.png", p4 + p5, width = 16, height = 6, dpi = 300)

# 步骤6: 学习轨迹图
cat("\n步骤6: 学习轨迹图\n")

# 学习轨迹
cds <- learn_graph(cds)

# 可视化轨迹图
p6 <- plot_cells(cds, color_cells_by = "cell_type", 
                label_cell_groups = TRUE,
                label_leaves = TRUE,
                label_branch_points = TRUE,
                graph_label_size = 3)

ggsave("plots/19_trajectory_graph.png", p6, width = 12, height = 10, dpi = 300)

# 步骤7: 排序细胞（伪时间分析）
cat("\n步骤7: 排序细胞（伪时间分析）\n")

# 选择根节点（起始点）
# 通常选择最原始的细胞类型作为起始点
# 这里我们选择Naive CD4 T细胞作为起始点

# 获取Naive CD4 T细胞
naive_cd4_cells <- colnames(cds)[colData(cds)$cell_type == "Naive CD4 T"]

# 自动选择根节点
cds <- order_cells(cds, root_cells = naive_cd4_cells[1:10])

# 可视化伪时间
p7 <- plot_cells(cds, color_cells_by = "pseudotime",
                label_cell_groups = TRUE,
                label_leaves = TRUE,
                label_branch_points = TRUE,
                graph_label_size = 3)

ggsave("plots/20_pseudotime.png", p7, width = 12, height = 10, dpi = 300)

# 步骤8: 分析轨迹相关基因
cat("\n步骤8: 分析轨迹相关基因\n")

# 寻找轨迹相关基因
tryCatch({
  # 这个分析可能需要较长时间
  cat("正在寻找轨迹相关基因（这可能需要几分钟）...\n")
  
  # 选择一些已知的T细胞相关基因进行分析
  t_cell_genes <- c("CD3D", "CD3E", "CD4", "CD8A", "IL7R", "CCR7", 
                   "CD44", "CD62L", "GZMB", "PRF1", "IFNG", "IL2")
  
  # 过滤存在的基因
  available_genes <- intersect(t_cell_genes, rownames(cds))
  
  if (length(available_genes) > 0) {
    # 可视化基因表达沿轨迹的变化
    p8 <- plot_cells(cds, genes = available_genes[1:4],
                    show_trajectory_graph = TRUE,
                    label_cell_groups = TRUE,
                    label_leaves = TRUE)
    
    ggsave("plots/21_genes_along_trajectory.png", p8, width = 12, height = 10, dpi = 300)
    
    # 绘制基因表达随伪时间的变化
    p9 <- plot_genes_in_pseudotime(cds[available_genes[1:6],],
                                  color_cells_by = "cell_type",
                                  min_expr = 0.5)
    
    ggsave("plots/22_genes_pseudotime_trend.png", p9, width = 12, height = 8, dpi = 300)
    
    cat("轨迹相关基因分析完成\n")
  } else {
    cat("未找到指定的T细胞标记基因\n")
  }
  
}, error = function(e) {
  cat("轨迹基因分析出错:", e$message, "\n")
  cat("跳过此步骤\n")
})

# 步骤9: 分支分析
cat("\n步骤9: 分支分析\n")

# 分析分支点
tryCatch({
  # 获取分支信息
  branch_info <- cds@principal_graph_aux$UMAP$branch_points
  
  if (length(branch_info) > 0) {
    cat("发现", length(branch_info), "个分支点\n")
    
    # 可视化分支
    p10 <- plot_cells(cds, color_cells_by = "cell_type",
                     label_branch_points = TRUE,
                     label_leaves = TRUE,
                     label_roots = TRUE,
                     graph_label_size = 4)
    
    ggsave("plots/23_branch_analysis.png", p10, width = 12, height = 10, dpi = 300)
  } else {
    cat("未发现明显的分支点\n")
  }
  
}, error = function(e) {
  cat("分支分析出错:", e$message, "\n")
})

# 步骤10: 保存轨迹分析结果
cat("\n步骤10: 保存轨迹分析结果\n")

# 保存Monocle对象
saveRDS(cds, "results/monocle_trajectory_analysis.rds")

# 提取伪时间信息并添加到原始Seurat对象
pseudotime_df <- data.frame(
  cell_id = colnames(cds),
  pseudotime = pseudotime(cds),
  cell_type = colData(cds)$cell_type
)

# 保存伪时间数据
write.csv(pseudotime_df, "results/pseudotime_data.csv", row.names = FALSE)

# 创建轨迹分析总结
trajectory_summary <- list(
  analysis_date = Sys.Date(),
  n_cells_analyzed = ncol(cds),
  cell_types = unique(colData(cds)$cell_type),
  n_clusters = length(unique(clusters(cds))),
  pseudotime_range = range(pseudotime(cds), na.rm = TRUE),
  method = "Monocle3"
)

saveRDS(trajectory_summary, "results/trajectory_summary.rds")

# 创建综合轨迹图
if (exists("p6") && exists("p7")) {
  trajectory_combined <- p6 + p7 + 
    plot_layout(ncol = 2) +
    plot_annotation(title = "T细胞分化轨迹分析",
                   subtitle = "左：轨迹图，右：伪时间")
  
  ggsave("plots/24_trajectory_combined.png", trajectory_combined, 
         width = 20, height = 8, dpi = 300)
}

# 打印轨迹分析总结
cat("\n=== 轨迹分析完成总结 ===\n")
cat("分析日期:", as.character(trajectory_summary$analysis_date), "\n")
cat("分析细胞数:", trajectory_summary$n_cells_analyzed, "\n")
cat("细胞类型:", paste(trajectory_summary$cell_types, collapse = ", "), "\n")
cat("Monocle聚类数:", trajectory_summary$n_clusters, "\n")
cat("伪时间范围:", round(trajectory_summary$pseudotime_range, 3), "\n")

cat("\n输出文件:\n")
cat("  - Monocle对象: results/monocle_trajectory_analysis.rds\n")
cat("  - 伪时间数据: results/pseudotime_data.csv\n")
cat("  - 轨迹总结: results/trajectory_summary.rds\n")
cat("  - 轨迹图表: plots/15-24_*.png\n")

cat("\n=== 轨迹分析完成！ ===\n")
cat("注意：轨迹分析结果需要结合生物学知识进行解释\n")
cat("建议进一步验证关键基因的表达模式和功能\n")

# 清理临时变量
rm(t_cell_types, expression_matrix, cell_metadata, gene_metadata)
rm(naive_cd4_cells, branch_info)
if (exists("available_genes")) rm(available_genes)
if (exists("t_cell_genes")) rm(t_cell_genes)

cat("轨迹分析脚本执行完成！\n")