# 第7次课实践操作手册

## 课程信息
- **课程名称**：高通量测序数据分析
- **主讲教师**：王运生
- **联系邮箱**：wangys@hunau.edu.cn
- **办公室**：16教420室
- **上课地点**：105机房
- **课程时间**：第7次课（2学时实践）

## 实验目标

### 主要目标
- 掌握单细胞RNA测序数据的基本分析流程
- 学会使用Seurat进行单细胞数据质量控制和预处理
- 掌握细胞聚类、降维可视化和标记基因识别方法
- 学会进行细胞类型注释和结果解读
- 了解轨迹分析的基本概念和方法

### 预期成果
- 完成一个完整的单细胞数据分析项目
- 生成高质量的UMAP可视化图表
- 识别不同细胞类型的标记基因
- 掌握单细胞数据分析的标准流程
- 能够独立进行基本的单细胞数据分析

## 环境要求

### 软件环境
| 软件名称 | 版本要求 | 安装方式 | 说明 |
|---------|---------|---------|------|
| R | ≥4.0.0 | 系统安装 | 统计计算环境 |
| RStudio | ≥1.4.0 | 官网下载 | R集成开发环境 |
| Seurat | ≥4.3.0 | `install.packages("Seurat")` | 单细胞分析主包 |
| dplyr | ≥1.0.0 | `install.packages("dplyr")` | 数据处理 |
| ggplot2 | ≥3.3.0 | `install.packages("ggplot2")` | 数据可视化 |
| patchwork | ≥1.1.0 | `install.packages("patchwork")` | 图表组合 |

### 硬件要求
- **内存**：至少 8 GB RAM（推荐16 GB）
- **存储空间**：至少 5 GB 可用空间
- **CPU**：多核处理器（推荐4核以上）
- **网络**：稳定的互联网连接（用于下载数据和包）

### 数据准备
| 数据文件 | 大小 | 下载链接/位置 | 说明 |
|---------|------|-------------|------|
| pbmc3k_filtered_gene_bc_matrices.tar.gz | ~20MB | 10X Genomics官网 | PBMC 3K数据集 |
| pbmc3k_final_seurat.rds | ~50MB | 课程共享文件夹 | 预处理完成的Seurat对象 |

## 操作步骤

### 步骤1：环境设置和准备工作

#### 1.1 创建工作目录
```bash
# 创建本次实验的工作目录
mkdir -p ~/ngs-analysis/lesson-07
cd ~/ngs-analysis/lesson-07

# 创建子目录结构
mkdir -p {data,scripts,results,plots}
```

#### 1.2 启动R环境并加载包
```r
# 启动RStudio或R控制台
# 设置工作目录
setwd("~/ngs-analysis/lesson-07")

# 加载必要的R包
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# 检查包版本
packageVersion("Seurat")
packageVersion("dplyr")
packageVersion("ggplot2")
```

**预期输出：**
```
[1] '4.3.0'
[1] '1.1.0'
[1] '3.4.0'
```

#### 1.3 下载和准备数据
```r
# 下载PBMC 3K数据集（如果尚未下载）
# 这是10X Genomics提供的经典示例数据集
# 包含约2700个外周血单核细胞(PBMC)

# 方法1：从10X官网下载
# url <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
# download.file(url, destfile = "data/pbmc3k_filtered_gene_bc_matrices.tar.gz")

# 方法2：使用课程提供的数据
# 假设数据已经解压到data目录
```

**检查点：** 确认数据文件已正确下载并位于 `data/` 目录中。

---

### 步骤2：数据导入和初步探索

#### 2.1 读取10X数据

**操作说明：**
使用Seurat的Read10X函数读取10X Genomics格式的数据，这种格式包含三个文件：matrix.mtx（表达矩阵）、features.tsv（基因信息）和barcodes.tsv（细胞条码）。

**执行命令：**
```r
# 读取10X数据
pbmc.data <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19/")

# 查看数据基本信息
dim(pbmc.data)
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# 创建Seurat对象
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                          project = "pbmc3k", 
                          min.cells = 3, 
                          min.features = 200)

# 查看Seurat对象信息
pbmc
```

**参数解释：**
- `min.cells = 3`：过滤在少于3个细胞中表达的基因
- `min.features = 200`：过滤检测到少于200个基因的细胞
- `project = "pbmc3k"`：项目名称标识

**预期输出：**
```
An object of class Seurat 
13714 features across 2700 samples within 1 assay 
Active assay: RNA (13714 features, 0 variable features)
```

#### 2.2 数据基本统计

**执行命令：**
```r
# 查看基因和细胞数量
cat("基因数量:", nrow(pbmc), "\n")
cat("细胞数量:", ncol(pbmc), "\n")

# 查看稀疏性
sparsity <- (sum(pbmc@assays$RNA@counts == 0) / length(pbmc@assays$RNA@counts)) * 100
cat("数据稀疏性:", round(sparsity, 2), "%\n")

# 查看表达量分布
summary(Matrix::rowSums(pbmc@assays$RNA@counts))
summary(Matrix::colSums(pbmc@assays$RNA@counts))
```

**结果验证：**
```r
# 验证数据读取是否正确
head(rownames(pbmc))  # 查看基因名
head(colnames(pbmc))  # 查看细胞条码
```

**检查点：** 确认数据已正确读取，基因和细胞数量符合预期。

---

### 步骤3：质量控制分析

#### 3.1 计算质量控制指标

**操作说明：**
计算每个细胞的质量控制指标，包括检测到的基因数量、总UMI数量、线粒体基因比例等。这些指标用于识别低质量细胞。

**执行命令：**
```r
# 计算线粒体基因比例
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 计算核糖体基因比例
pbmc[["percent.rb"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")

# 查看质控指标
head(pbmc@meta.data, 5)

# 基本统计
summary(pbmc@meta.data$nFeature_RNA)
summary(pbmc@meta.data$nCount_RNA)
summary(pbmc@meta.data$percent.mt)
```

#### 3.2 质量控制可视化

**执行命令：**
```r
# 创建质控指标的小提琴图
p1 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, pt.size = 0.1)
print(p1)

# 保存图片
ggsave("plots/qc_violin_plots.png", p1, width = 12, height = 4, dpi = 300)

# 创建散点图查看指标间关系
p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
p3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2 + p3

# 保存散点图
ggsave("plots/qc_scatter_plots.png", p2 + p3, width = 10, height = 4, dpi = 300)
```

#### 3.3 过滤低质量细胞

**执行命令：**
```r
# 设置过滤阈值
# 这些阈值需要根据数据特点调整
min_features <- 200   # 最少检测基因数
max_features <- 2500  # 最多检测基因数（去除可能的双细胞）
max_mt_percent <- 20  # 最大线粒体基因比例

# 应用过滤条件
pbmc <- subset(pbmc, subset = nFeature_RNA > min_features & 
                             nFeature_RNA < max_features & 
                             percent.mt < max_mt_percent)

# 查看过滤后的细胞数量
cat("过滤前细胞数量:", ncol(pbmc.data), "\n")
cat("过滤后细胞数量:", ncol(pbmc), "\n")
cat("过滤掉的细胞数量:", ncol(pbmc.data) - ncol(pbmc), "\n")
```

**检查点：** 确认过滤参数合理，保留了足够数量的高质量细胞。

---

### 步骤4：数据标准化和特征选择

#### 4.1 数据标准化

**操作说明：**
对基因表达数据进行标准化，消除细胞间测序深度的差异，使不同细胞的表达水平可以进行比较。

**执行命令：**
```r
# 标准化数据
# 使用LogNormalize方法，将每个细胞的表达量标准化到10000，然后取对数
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", 
                     scale.factor = 10000)

# 查看标准化后的数据
pbmc@assays$RNA@data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
```

#### 4.2 寻找高变基因

**执行命令：**
```r
# 识别高变基因
# 这些基因在细胞间表达变异最大，最能反映细胞异质性
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 获取前10个高变基因
top10 <- head(VariableFeatures(pbmc), 10)
print(top10)

# 可视化高变基因
p4 <- VariableFeaturePlot(pbmc)
p5 <- LabelPoints(plot = p4, points = top10, repel = TRUE)
print(p5)

# 保存高变基因图
ggsave("plots/variable_features.png", p5, width = 10, height = 6, dpi = 300)
```

#### 4.3 数据缩放

**执行命令：**
```r
# 缩放数据
# 对所有基因进行Z-score标准化，使每个基因的均值为0，方差为1
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# 查看缩放后的数据
pbmc@assays$RNA@scale.data[c("CD3D", "TCL1A", "MS4A1"), 1:5]
```

**检查点：** 确认数据标准化和缩放步骤完成，高变基因已正确识别。

---

### 步骤5：主成分分析(PCA)

#### 5.1 执行PCA分析

**操作说明：**
使用主成分分析对高变基因进行降维，识别数据中的主要变异源。PCA结果将用于后续的聚类和非线性降维分析。

**执行命令：**
```r
# 执行PCA分析
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# 查看PCA结果
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# 可视化PCA结果
p6 <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
print(p6)

# PCA散点图
p7 <- DimPlot(pbmc, reduction = "pca")
print(p7)

# 保存PCA图
ggsave("plots/pca_plot.png", p7, width = 8, height = 6, dpi = 300)
```

#### 5.2 确定主成分数量

**执行命令：**
```r
# 绘制肘部图确定主成分数量
p8 <- ElbowPlot(pbmc, ndims = 50)
print(p8)

# 保存肘部图
ggsave("plots/elbow_plot.png", p8, width = 8, height = 6, dpi = 300)

# 基于肘部图，选择前15个主成分用于后续分析
n_pcs <- 15
cat("选择的主成分数量:", n_pcs, "\n")
```

**检查点：** 确认PCA分析完成，选择了合适数量的主成分。

---

### 步骤6：细胞聚类分析

#### 6.1 构建邻接图

**操作说明：**
基于PCA结果构建细胞间的邻接图，为聚类分析做准备。使用共享最近邻(SNN)图方法。

**执行命令：**
```r
# 构建SNN图
pbmc <- FindNeighbors(pbmc, dims = 1:n_pcs)

# 查看邻接图信息
pbmc@graphs$RNA_snn[1:5, 1:5]
```

#### 6.2 执行聚类

**执行命令：**
```r
# 执行Louvain聚类
# resolution参数控制聚类的粒度，值越大聚类越细
pbmc <- FindClusters(pbmc, resolution = 0.5)

# 查看聚类结果
head(Idents(pbmc), 5)
table(Idents(pbmc))

# 查看聚类数量
n_clusters <- length(levels(Idents(pbmc)))
cat("聚类数量:", n_clusters, "\n")
```

**参数解释：**
- `resolution = 0.5`：聚类分辨率，可以尝试0.1-2.0之间的值
- 较小的值产生较少的聚类，较大的值产生较多的聚类

**检查点：** 确认聚类分析完成，聚类数量合理（通常5-20个聚类）。

---

### 步骤7：非线性降维可视化

#### 7.1 UMAP降维

**操作说明：**
使用UMAP（Uniform Manifold Approximation and Projection）进行非线性降维，将高维数据投影到二维空间进行可视化。

**执行命令：**
```r
# 执行UMAP降维
pbmc <- RunUMAP(pbmc, dims = 1:n_pcs)

# 可视化UMAP结果
p9 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
      ggtitle("UMAP聚类结果")
print(p9)

# 保存UMAP图
ggsave("plots/umap_clusters.png", p9, width = 10, height = 8, dpi = 300)
```

#### 7.2 t-SNE降维（可选）

**执行命令：**
```r
# 执行t-SNE降维进行比较
pbmc <- RunTSNE(pbmc, dims = 1:n_pcs)

# 可视化t-SNE结果
p10 <- DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + 
       ggtitle("t-SNE聚类结果")
print(p10)

# 比较UMAP和t-SNE
p9 + p10

# 保存比较图
ggsave("plots/umap_tsne_comparison.png", p9 + p10, width = 16, height = 6, dpi = 300)
```

**检查点：** 确认降维可视化完成，聚类在二维空间中分离良好。

---

### 步骤8：寻找标记基因

#### 8.1 差异表达分析

**操作说明：**
为每个聚类寻找标记基因，这些基因在特定聚类中高表达，用于细胞类型注释。

**执行命令：**
```r
# 寻找所有聚类的标记基因
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25)

# 查看标记基因结果
head(pbmc.markers, 10)

# 获取每个聚类的top基因
top_markers <- pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

print(top_markers)
```

**参数解释：**
- `only.pos = TRUE`：只寻找正向标记基因（高表达）
- `min.pct = 0.25`：基因至少在25%的细胞中表达
- `logfc.threshold = 0.25`：log2倍数变化阈值

#### 8.2 可视化标记基因

**执行命令：**
```r
# 可视化top标记基因
p11 <- VlnPlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", 
                                 "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                                 "CD8A"))
print(p11)

# 保存小提琴图
ggsave("plots/marker_genes_violin.png", p11, width = 15, height = 10, dpi = 300)

# 在UMAP上显示标记基因表达
p12 <- FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", 
                                     "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                                     "CD8A"))
print(p12)

# 保存特征图
ggsave("plots/marker_genes_feature.png", p12, width = 15, height = 12, dpi = 300)
```

#### 8.3 热图可视化

**执行命令：**
```r
# 创建标记基因热图
top10_markers <- pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)

p13 <- DoHeatmap(pbmc, features = top10_markers$gene) + NoLegend()
print(p13)

# 保存热图
ggsave("plots/marker_genes_heatmap.png", p13, width = 12, height = 15, dpi = 300)
```

**检查点：** 确认标记基因分析完成，每个聚类都有特异性表达的基因。

---

### 步骤9：细胞类型注释

#### 9.1 基于标记基因的手动注释

**操作说明：**
根据已知的细胞类型标记基因，对聚类进行生物学注释。这需要结合文献知识和数据库信息。

**执行命令：**
```r
# 基于标记基因进行细胞类型注释
# 这里使用PBMC数据的经典注释
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", 
                    "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# 可视化注释结果
p14 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
       ggtitle("细胞类型注释结果")
print(p14)

# 保存注释结果
ggsave("plots/cell_type_annotation.png", p14, width = 12, height = 8, dpi = 300)
```

#### 9.2 验证注释结果

**执行命令：**
```r
# 验证关键标记基因的表达
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
print(p15)

# 保存点图
ggsave("plots/canonical_markers_dotplot.png", p15, width = 12, height = 6, dpi = 300)
```

#### 9.3 细胞类型比例统计

**执行命令：**
```r
# 统计各细胞类型的数量和比例
cell_type_counts <- table(Idents(pbmc))
cell_type_props <- prop.table(cell_type_counts) * 100

# 创建统计表
cell_type_stats <- data.frame(
  CellType = names(cell_type_counts),
  Count = as.numeric(cell_type_counts),
  Percentage = round(as.numeric(cell_type_props), 2)
)

print(cell_type_stats)

# 保存统计结果
write.csv(cell_type_stats, "results/cell_type_statistics.csv", row.names = FALSE)
```

**检查点：** 确认细胞类型注释合理，标记基因表达模式与预期一致。

---

### 步骤10：保存分析结果

#### 10.1 保存Seurat对象

**执行命令：**
```r
# 保存完整的Seurat对象
saveRDS(pbmc, file = "results/pbmc_analysis_complete.rds")

# 保存标记基因结果
write.csv(pbmc.markers, "results/all_markers.csv", row.names = FALSE)
write.csv(top_markers, "results/top_markers.csv", row.names = FALSE)

cat("分析结果已保存到results目录\n")
```

#### 10.2 生成分析报告

**执行命令：**
```r
# 创建分析总结
analysis_summary <- list(
  total_cells = ncol(pbmc),
  total_genes = nrow(pbmc),
  n_clusters = length(levels(Idents(pbmc))),
  cell_types = levels(Idents(pbmc)),
  top_variable_genes = head(VariableFeatures(pbmc), 10)
)

# 保存分析总结
saveRDS(analysis_summary, "results/analysis_summary.rds")

# 打印总结信息
cat("=== 单细胞数据分析总结 ===\n")
cat("细胞总数:", analysis_summary$total_cells, "\n")
cat("基因总数:", analysis_summary$total_genes, "\n")
cat("聚类数量:", analysis_summary$n_clusters, "\n")
cat("识别的细胞类型:\n")
for(ct in analysis_summary$cell_types) {
  cat("  -", ct, "\n")
}
```

**检查点：** 确认所有结果文件已正确保存。

---

## 预期结果

### 主要输出文件
1. **Seurat对象**：`results/pbmc_analysis_complete.rds`
   - 内容：完整的分析结果，包含所有计算数据
   - 用途：后续分析和结果重现

2. **标记基因列表**：`results/all_markers.csv`
   - 内容：所有聚类的差异表达基因
   - 用途：细胞类型注释和功能分析

3. **可视化图表**：`plots/`目录下的PNG文件
   - 内容：质控图、UMAP图、标记基因图等
   - 用途：结果展示和报告制作

### 关键结果指标
- **细胞数量**：应该在2000-2700之间（过滤后）
- **聚类数量**：预期8-10个主要聚类
- **细胞类型**：识别出T细胞、B细胞、单核细胞、NK细胞等
- **数据质量**：线粒体基因比例<20%，基因检测数200-2500

### 成功标准
- [ ] 数据成功导入并通过质量控制
- [ ] 聚类结果在UMAP图上分离良好
- [ ] 每个聚类都有特异性标记基因
- [ ] 细胞类型注释与已知生物学知识一致
- [ ] 所有图表清晰美观，适合展示

## 故障排除

### 常见问题1：内存不足错误
**症状：** R报告内存不足，无法完成计算
**原因：** 数据量大，系统内存不够
**解决方案：**
```r
# 增加内存限制（Windows）
memory.limit(size = 16000)  # 设置为16GB

# 或者减少分析的基因数量
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))  # 只缩放高变基因
```r

### 常见问题2：包加载失败
**症状：** library()命令报错，提示包不存在
**原因：** R包未正确安装
**解决方案：**
```r
# 重新安装问题包
install.packages("Seurat")
install.packages("dplyr")

# 如果仍有问题，尝试从Bioconductor安装
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Seurat")
```

### 常见问题3：图片保存失败
**症状：** ggsave()命令报错
**原因：** 目录不存在或权限问题
**解决方案：**
```r
# 确保目录存在
dir.create("plots", showWarnings = FALSE)

# 检查当前工作目录
getwd()

# 如果需要，更改工作目录
setwd("~/ngs-analysis/lesson-07")
```

### 常见问题4：聚类结果不理想
**症状：** 聚类过多或过少，细胞类型混杂
**原因：** 参数设置不当
**解决方案：**
```r
# 调整聚类分辨率
pbmc <- FindClusters(pbmc, resolution = 0.3)  # 减少聚类数
# 或
pbmc <- FindClusters(pbmc, resolution = 0.8)  # 增加聚类数

# 调整PCA维度
pbmc <- FindNeighbors(pbmc, dims = 1:10)  # 使用更少的PC
```r

### 获取帮助
如果遇到其他问题：
1. 查看R控制台的错误信息
2. 检查Seurat官方文档：https://satijalab.org/seurat/
3. 搜索相关错误信息和解决方案
4. 联系助教或老师：wangys@hunau.edu.cn

## 扩展练习

### 练习1：参数优化
**目标：** 尝试不同的参数设置，观察对结果的影响
**任务：** 
- 尝试不同的聚类分辨率（0.1, 0.3, 0.5, 0.8, 1.0）
- 比较使用不同数量主成分的效果
- 调整质量控制的阈值参数
**提示：** 记录不同参数下的聚类数量和细胞类型注释结果

### 练习2：深入的标记基因分析
**目标：** 进行更详细的标记基因分析和功能注释
**任务：**
- 为每个细胞类型寻找更多的标记基因
- 使用在线数据库验证标记基因的生物学功能
- 创建更详细的标记基因热图和点图
**提示：** 可以使用CellMarker、PanglaoDB等数据库

### 练习3：轨迹分析入门
**目标：** 尝试使用Monocle3进行简单的轨迹分析
**任务：**
- 安装并加载Monocle3包
- 将Seurat对象转换为Monocle对象
- 进行简单的轨迹分析
**提示：** 重点关注T细胞或单核细胞的分化轨迹

### 思考问题
1. 为什么单细胞数据需要特殊的标准化方法？
2. 如何判断聚类结果的好坏？
3. 细胞类型注释的准确性如何验证？
4. 单细胞分析中的批次效应如何处理？
5. 如何选择合适的降维方法和参数？

## 参考资料

### 相关文献
1. Stuart T, et al. Comprehensive Integration of Single-Cell Data. Cell, 2019.
2. Hao Y, et al. Integrated analysis of multimodal single-cell data. Cell, 2021.
3. Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis. Mol Syst Biol, 2019.

### 在线资源
- Seurat官方教程：https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
- 单细胞分析最佳实践：https://www.sc-best-practices.org/
- 10X Genomics数据集：https://www.10xgenomics.com/resources/datasets

### 软件文档
- Seurat文档：https://satijalab.org/seurat/
- R语言官方文档：https://www.r-project.org/
- ggplot2文档：https://ggplot2.tidyverse.org/

### 数据库资源
- CellMarker：http://biocc.hrbmu.edu.cn/CellMarker/
- PanglaoDB：https://panglaodb.se/
- Human Cell Atlas：https://www.humancellatlas.org/

## 附录

### 附录A：完整分析脚本
参见：`scripts/`目录中的`complete_analysis.R`文件

### 附录B：质量控制参数建议
```r
# 不同组织类型的QC参数建议
qc_params <- list(
  PBMC = list(min_features = 200, max_features = 2500, max_mt = 20),
  Brain = list(min_features = 200, max_features = 4000, max_mt = 5),
  Kidney = list(min_features = 200, max_features = 3000, max_mt = 25),
  Liver = list(min_features = 200, max_features = 3500, max_mt = 30)
)
```

### 附录C：常用标记基因列表
```r
# PBMC细胞类型标记基因
marker_genes <- list(
  "T_cells" = c("CD3D", "CD3E", "CD3G"),
  "CD4_T" = c("CD4", "IL7R", "CCR7"),
  "CD8_T" = c("CD8A", "CD8B"),
  "B_cells" = c("CD19", "MS4A1", "CD79A"),
  "Monocytes" = c("CD14", "LYZ", "S100A9"),
  "NK_cells" = c("GNLY", "NKG7", "KLRB1"),
  "Dendritic" = c("FCER1A", "CST3"),
  "Platelets" = c("PPBP", "PF4")
)
```

---

**实验完成时间：** 预计 2 小时  
**难度等级：** 中级  
**最后更新：** 2025年1月