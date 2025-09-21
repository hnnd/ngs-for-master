# 单细胞RNA测序数据完整分析脚本
# 课程：高通量测序数据分析 - 第7次课
# 作者：王运生
# 日期：2025-01-21
# 用法：source("scripts/complete_analysis.R")

# 加载环境设置
if (file.exists("scripts/setup.R")) {
  source("scripts/setup.R")
} else {
  cat("警告：未找到setup.R文件，请先运行环境设置\n")
}

cat("\n=== 开始单细胞RNA测序数据分析 ===\n")

# 步骤1: 数据导入
cat("\n步骤1: 数据导入\n")

# 检查数据文件是否存在
data_path <- "data/filtered_gene_bc_matrices/hg19/"
if (!dir.exists(data_path)) {
  cat("数据目录不存在，尝试其他路径...\n")
  # 尝试其他可能的路径
  alternative_paths <- c(
    "data/pbmc3k/filtered_gene_bc_matrices/hg19/",
    "data/filtered_gene_bc_matrices/",
    "data/"
  )
  
  for (path in alternative_paths) {
    if (dir.exists(path) && length(list.files(path, pattern = "matrix.mtx|barcodes.tsv|features.tsv")) > 0) {
      data_path <- path
      cat("找到数据路径:", data_path, "\n")
      break
    }
  }
}

# 读取数据
tryCatch({
  pbmc.data <- Read10X(data.dir = data_path)
  cat("成功读取数据，维度:", dim(pbmc.data), "\n")
}, error = function(e) {
  cat("数据读取失败:", e$message, "\n")
  cat("请检查数据文件是否存在于正确路径\n")
  stop("数据读取失败")
})

# 创建Seurat对象
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                          project = "pbmc3k", 
                          min.cells = 3, 
                          min.features = 200)

cat("Seurat对象创建完成\n")
print(pbmc)

# 步骤2: 质量控制
cat("\n步骤2: 质量控制分析\n")

# 计算质控指标
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc[["percent.rb"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")

# 质控统计
cat("质控指标统计:\n")
cat("基因数量 - 中位数:", median(pbmc$nFeature_RNA), "范围:", range(pbmc$nFeature_RNA), "\n")
cat("UMI数量 - 中位数:", median(pbmc$nCount_RNA), "范围:", range(pbmc$nCount_RNA), "\n")
cat("线粒体% - 中位数:", round(median(pbmc$percent.mt), 2), "范围:", round(range(pbmc$percent.mt), 2), "\n")

# 质控可视化
p1 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, pt.size = 0.1)
ggsave("plots/01_qc_violin_plots.png", p1, width = 12, height = 4, dpi = 300)

p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
p3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("plots/02_qc_scatter_plots.png", p2 + p3, width = 10, height = 4, dpi = 300)

# 过滤细胞
cells_before <- ncol(pbmc)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & 
                             nFeature_RNA < 2500 & 
                             percent.mt < 20)
cells_after <- ncol(pbmc)

cat("细胞过滤结果:\n")
cat("过滤前:", cells_before, "细胞\n")
cat("过滤后:", cells_after, "细胞\n")
cat("过滤掉:", cells_before - cells_after, "细胞 (", 
    round((cells_before - cells_after)/cells_before * 100, 1), "%)\n")

# 步骤3: 数据标准化
cat("\n步骤3: 数据标准化\n")

# 标准化
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", 
                     scale.factor = 10000)

# 寻找高变基因
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 可视化高变基因
top10 <- head(VariableFeatures(pbmc), 10)
cat("前10个高变基因:", paste(top10, collapse = ", "), "\n")

p4 <- VariableFeaturePlot(pbmc)
p5 <- LabelPoints(plot = p4, points = top10, repel = TRUE)
ggsave("plots/03_variable_features.png", p5, width = 10, height = 6, dpi = 300)

# 数据缩放
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# 步骤4: 主成分分析
cat("\n步骤4: 主成分分析\n")

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# PCA可视化
p6 <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
ggsave("plots/04_pca_loadings.png", p6, width = 10, height = 6, dpi = 300)

p7 <- DimPlot(pbmc, reduction = "pca")
ggsave("plots/05_pca_plot.png", p7, width = 8, height = 6, dpi = 300)

# 肘部图
p8 <- ElbowPlot(pbmc, ndims = 50)
ggsave("plots/06_elbow_plot.png", p8, width = 8, height = 6, dpi = 300)

# 选择主成分数量
n_pcs <- 15
cat("选择的主成分数量:", n_pcs, "\n")

# 步骤5: 细胞聚类
cat("\n步骤5: 细胞聚类\n")

# 构建邻接图
pbmc <- FindNeighbors(pbmc, dims = 1:n_pcs)

# 聚类
pbmc <- FindClusters(pbmc, resolution = 0.5)

# 聚类统计
cluster_counts <- table(Idents(pbmc))
cat("聚类结果:\n")
print(cluster_counts)

# 步骤6: 非线性降维
cat("\n步骤6: 非线性降维\n")

# UMAP
pbmc <- RunUMAP(pbmc, dims = 1:n_pcs)

p9 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
      ggtitle("UMAP聚类结果")
ggsave("plots/07_umap_clusters.png", p9, width = 10, height = 8, dpi = 300)

# t-SNE（可选）
pbmc <- RunTSNE(pbmc, dims = 1:n_pcs)

p10 <- DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + 
       ggtitle("t-SNE聚类结果")

# 比较UMAP和t-SNE
ggsave("plots/08_umap_tsne_comparison.png", p9 + p10, width = 16, height = 6, dpi = 300)

# 步骤7: 寻找标记基因
cat("\n步骤7: 寻找标记基因\n")

# 寻找所有聚类的标记基因
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25)

cat("找到", nrow(pbmc.markers), "个标记基因\n")

# 获取每个聚类的top基因
top_markers <- pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

cat("每个聚类的top标记基因:\n")
print(top_markers)

# 可视化标记基因
marker_genes_to_plot <- c("MS4A1", "GNLY", "CD3E", "CD14", 
                         "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")

p11 <- VlnPlot(pbmc, features = marker_genes_to_plot)
ggsave("plots/09_marker_genes_violin.png", p11, width = 15, height = 10, dpi = 300)

p12 <- FeaturePlot(pbmc, features = marker_genes_to_plot)
ggsave("plots/10_marker_genes_feature.png", p12, width = 15, height = 12, dpi = 300)

# 热图
top10_markers <- pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)

p13 <- DoHeatmap(pbmc, features = top10_markers$gene) + NoLegend()
ggsave("plots/11_marker_genes_heatmap.png", p13, width = 12, height = 15, dpi = 300)

# 步骤8: 细胞类型注释
cat("\n步骤8: 细胞类型注释\n")

# 基于标记基因的注释
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", 
                    "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# 可视化注释结果
p14 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
       ggtitle("细胞类型注释结果")
ggsave("plots/12_cell_type_annotation.png", p14, width = 12, height = 8, dpi = 300)

# 验证标记基因
canonical_markers <- c("IL7R", "CCR7",    # Naive CD4+ T
                      "CD14", "LYZ",      # CD14+ Mono
                      "IL7R", "S100A4",   # Memory CD4+ T
                      "MS4A1",            # B
                      "CD8A",             # CD8+ T
                      "FCGR3A", "MS4A7",  # FCGR3A+ Mono
                      "GNLY", "NKG7",     # NK
                      "FCER1A", "CST3",   # DC
                      "PPBP")             # Platelet

p15 <- DotPlot(pbmc, features = canonical_markers) + RotatedAxis()
ggsave("plots/13_canonical_markers_dotplot.png", p15, width = 12, height = 6, dpi = 300)

# 细胞类型统计
cell_type_counts <- table(Idents(pbmc))
cell_type_props <- prop.table(cell_type_counts) * 100

cell_type_stats <- data.frame(
  CellType = names(cell_type_counts),
  Count = as.numeric(cell_type_counts),
  Percentage = round(as.numeric(cell_type_props), 2)
)

cat("细胞类型统计:\n")
print(cell_type_stats)

# 步骤9: 保存结果
cat("\n步骤9: 保存分析结果\n")

# 保存Seurat对象
saveRDS(pbmc, file = "results/pbmc_analysis_complete.rds")

# 保存标记基因
write.csv(pbmc.markers, "results/all_markers.csv", row.names = FALSE)
write.csv(top_markers, "results/top_markers.csv", row.names = FALSE)
write.csv(cell_type_stats, "results/cell_type_statistics.csv", row.names = FALSE)

# 创建分析总结
analysis_summary <- list(
  analysis_date = Sys.Date(),
  total_cells = ncol(pbmc),
  total_genes = nrow(pbmc),
  cells_after_qc = cells_after,
  cells_filtered = cells_before - cells_after,
  n_clusters = length(levels(Idents(pbmc))),
  cell_types = levels(Idents(pbmc)),
  n_markers = nrow(pbmc.markers),
  top_variable_genes = head(VariableFeatures(pbmc), 10),
  parameters = list(
    min_features = 200,
    max_features = 2500,
    max_mt_percent = 20,
    n_variable_features = 2000,
    n_pcs = n_pcs,
    clustering_resolution = 0.5
  )
)

saveRDS(analysis_summary, "results/analysis_summary.rds")

# 生成最终报告图
final_plots <- list(
  qc = p1,
  umap_clusters = p9,
  cell_types = p14,
  markers = p15
)

# 创建综合图表
summary_plot <- (p1 / (p9 + p14)) | p15
ggsave("plots/14_analysis_summary.png", summary_plot, width = 20, height = 12, dpi = 300)

# 打印最终总结
cat("\n=== 分析完成总结 ===\n")
cat("分析日期:", as.character(analysis_summary$analysis_date), "\n")
cat("细胞总数:", analysis_summary$total_cells, "\n")
cat("基因总数:", analysis_summary$total_genes, "\n")
cat("质控后细胞数:", analysis_summary$cells_after_qc, "\n")
cat("聚类数量:", analysis_summary$n_clusters, "\n")
cat("标记基因数量:", analysis_summary$n_markers, "\n")

cat("\n识别的细胞类型:\n")
for(i in 1:nrow(cell_type_stats)) {
  cat(sprintf("  %s: %d cells (%.1f%%)\n", 
              cell_type_stats$CellType[i], 
              cell_type_stats$Count[i], 
              cell_type_stats$Percentage[i]))
}

cat("\n输出文件:\n")
cat("  - Seurat对象: results/pbmc_analysis_complete.rds\n")
cat("  - 标记基因: results/all_markers.csv\n")
cat("  - 细胞统计: results/cell_type_statistics.csv\n")
cat("  - 分析总结: results/analysis_summary.rds\n")
cat("  - 图表文件: plots/目录下的PNG文件\n")

cat("\n=== 单细胞RNA测序数据分析完成！ ===\n")

# 清理环境
rm(pbmc.data, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15)
rm(top10, marker_genes_to_plot, canonical_markers, new.cluster.ids)
rm(cells_before, cells_after, n_pcs, cluster_counts, cell_type_counts, cell_type_props)
rm(top_markers, top10_markers, final_plots, summary_plot)

cat("环境已清理，主要对象pbmc和分析结果已保留\n")