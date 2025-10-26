# 第5次课实践操作手册

## 课程信息
- **课程名称**：高通量测序数据分析
- **主讲教师**：王运生
- **联系邮箱**：wangys@hunau.edu.cn
- **办公室**：16教420室
- **上课地点**：105机房
- **课程时间**：第5次课

## 实验目标

### 主要目标
- 掌握RNA-seq数据分析的完整流程
- 学会使用HISAT2进行转录组序列比对
- 掌握featureCounts进行基因表达定量
- 学会使用DESeq2进行差异表达分析
- 掌握结果可视化和功能富集分析方法

### 预期成果
- 获得高质量的转录组比对结果
- 生成基因表达定量矩阵
- 识别差异表达基因列表
- 创建MA图、火山图等可视化结果
- 完成GO功能富集分析

## 环境要求

### 软件环境
| 软件名称 | 版本要求 | 安装方式 | 说明 |
|---------|---------|---------|------|
| FastQC | ≥0.11.9 | `conda install fastqc` | 质量控制 |
| MultiQC | ≥1.9 | `conda install multiqc` | 报告汇总 |
| Trimmomatic | ≥0.39 | `conda install trimmomatic` | 序列清洗 |
| HISAT2 | ≥2.2.1 | `conda install hisat2` | 序列比对 |
| SAMtools | ≥1.12 | `conda install samtools` | SAM/BAM处理 |
| featureCounts | ≥2.0.1 | `conda install subread` | 基因定量 |
| R | ≥4.1.0 | `conda install r-base` | 统计分析 |
| DESeq2 | ≥1.32.0 | R包安装 | 差异表达分析 |

### 硬件要求
- **内存**：至少 8 GB RAM
- **存储空间**：至少 20 GB 可用空间
- **CPU**：4核心以上推荐
- **网络**：稳定的互联网连接（用于下载数据）

### 数据准备

#### 📥 推荐方式：使用自动化脚本

```bash
# 进入manual目录
cd lesson-05/manual

# 自动准备所有数据
bash scripts/download_data.sh

# 该脚本将自动完成：
# 1. 检查软件依赖
# 2. 下载/生成参考基因组
# 3. 准备测序数据（可选择下载或生成）
# 4. 验证数据完整性
```

#### 📊 数据选项

| 数据类型 | 大小 | 准备时间 | 说明 |
|---------|------|---------|------|
| 演示数据 | ~300MB | 5-10分钟 | 快速学习，推荐初学者 |
| 模拟数据 | ~2GB | 20-30分钟 | 完整分析，使用wgsim生成 |
| 真实数据 | ~10-20GB | 1-2小时 | NCBI SRA数据集 |

#### ℹ️ 更多详细信息

- **自动化脚本选项**: 查看脚本帮助
  ```bash
  bash scripts/download_data.sh --help
  ```

- **完整数据准备指南**: 参考 [../DATA_SOURCES.md](../DATA_SOURCES.md)

## 操作步骤

### 步骤1：环境设置和准备工作

#### 1.1 创建工作目录
```bash
# 创建本次实验的工作目录
mkdir -p ~/ngs-analysis/lesson-05
cd ~/ngs-analysis/lesson-05

# 创建子目录结构
mkdir -p {data,reference,scripts,results,logs,qc}
mkdir -p results/{alignment,quantification,deseq2,plots}
```

#### 1.2 检查软件环境
```bash
# 检查必要软件是否已安装
fastqc --version
multiqc --version
trimmomatic -version
hisat2 --version
samtools --version
featureCounts -v
R --version
```

**预期输出：**
```
FastQC v0.11.9
multiqc, version 1.9
0.39
version 2.2.1
samtools 1.12
featureCounts v2.0.1
R version 4.1.0
```

#### 1.3 使用脚本自动准备数据（推荐）

```bash
# 进入manual目录
cd ~/ngs-analysis/lesson-05

# 运行自动化数据准备脚本
bash manual/scripts/download_data.sh --all

# 脚本将自动：
# ✓ 检查必需的软件工具
# ✓ 下载/生成参考基因组（多个选项可选）
# ✓ 生成/下载测序数据
# ✓ 验证数据的完整性
# ✓ 生成样本信息表
```

**检查点：** 查看日志文件确认数据准备成功

```bash
cat ../logs/data_summary.txt
ls -lh manual/reference/
ls -lh manual/data/
```

#### 1.4 手动准备数据（可选）

如果自动脚本无法运行或需要特定数据，可参考详细指南：

- **完整操作步骤**: 见 [../DATA_SOURCES.md](../DATA_SOURCES.md)
- **参考基因组下载**: 选择Ensembl或UCSC源
- **测序数据获取**:
  - 直接下载NCBI SRA数据
  - 使用wgsim生成模拟数据
  - 使用课程服务器的预置数据

---

### 步骤2：数据质量控制

#### 2.1 原始数据质量评估

**操作说明：**
使用FastQC对原始测序数据进行质量评估，检查序列质量、GC含量、重复序列等指标。

**执行命令：**
```bash
# 回到工作目录
cd ~/ngs-analysis/lesson-05

# 对所有fastq文件进行质量评估
fastqc data/*.fastq.gz -o qc/ -t 4

# 生成汇总报告
multiqc qc/ -o qc/multiqc_report/
```

**参数解释：**
- `-o`：指定输出目录
- `-t 4`：使用4个线程并行处理
- `multiqc`：整合多个FastQC报告

**预期输出：**
```
Analysis complete for sample1_R1.fastq.gz
Analysis complete for sample1_R2.fastq.gz
...
[INFO] multiqc: Report saved to qc/multiqc_report/multiqc_report.html
```

**结果验证：**
```bash
# 查看生成的报告文件
ls qc/
firefox qc/multiqc_report/multiqc_report.html &
```

#### 2.2 数据清洗（如需要）

**操作说明：**
根据质量评估结果，使用Trimmomatic去除低质量序列和接头序列。

**执行命令：**
```bash
# 创建清洗后数据目录
mkdir -p data/clean

# 对每个样本进行数据清洗
for sample in sample1 sample2 sample3 ctrl1 ctrl2 ctrl3; do
    trimmomatic PE -threads 4 \
        data/${sample}_R1.fastq.gz data/${sample}_R2.fastq.gz \
        data/clean/${sample}_R1_clean.fastq.gz data/clean/${sample}_R1_unpaired.fastq.gz \
        data/clean/${sample}_R2_clean.fastq.gz data/clean/${sample}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
```

**参数解释：**
- `PE`：双端测序模式
- `ILLUMINACLIP`：去除Illumina接头序列
- `LEADING:3`：去除开头质量值低于3的碱基
- `TRAILING:3`：去除结尾质量值低于3的碱基
- `SLIDINGWINDOW:4:15`：滑动窗口质量过滤
- `MINLEN:36`：保留长度至少36bp的序列

**检查点：** 确认清洗后的数据质量有明显改善。

---

### 步骤3：建立参考基因组索引

#### 3.1 构建HISAT2索引

**操作说明：**
为参考基因组构建HISAT2比对索引，这是进行序列比对的前提步骤。

**执行命令：**
```bash
# 进入参考基因组目录
cd reference

# 构建HISAT2索引（这个过程需要较长时间）
hisat2-build -p 4 Homo_sapiens.GRCh38.dna.primary_assembly.fa grch38_index

# 检查索引文件
ls grch38_index.*
```

**参数解释：**
- `-p 4`：使用4个线程
- `grch38_index`：索引文件前缀

**预期输出：**
```
Building a SMALL index
Total time for call to driver() for forward index: 00:45:23
Total time for call to driver() for mirror index: 00:42:18
```

**结果验证：**
```bash
# 验证索引文件完整性
ls -lh grch38_index.*
# 应该看到8个索引文件（.1.ht2, .2.ht2, 等）
```

**检查点：** 确认所有索引文件已成功生成。

---

### 步骤4：序列比对

#### 4.1 HISAT2序列比对

**操作说明：**
使用HISAT2将清洗后的RNA-seq reads比对到参考基因组，生成SAM格式的比对文件。

**执行命令：**
```bash
# 回到工作目录
cd ~/ngs-analysis/lesson-05

# 对每个样本进行序列比对
for sample in sample1 sample2 sample3 ctrl1 ctrl2 ctrl3; do
    echo "Processing $sample..."
    hisat2 -x reference/grch38_index \
        -1 data/clean/${sample}_R1_clean.fastq.gz \
        -2 data/clean/${sample}_R2_clean.fastq.gz \
        -S results/alignment/${sample}.sam \
        --threads 4 \
        --rna-strandness RF \
        --summary-file logs/${sample}_alignment_summary.txt
done
```

**参数解释：**
- `-x`：指定索引文件前缀
- `-1/-2`：双端测序的两个文件
- `-S`：输出SAM文件
- `--threads 4`：使用4个线程
- `--rna-strandness RF`：链特异性参数
- `--summary-file`：保存比对统计信息

**预期输出：**
```
1000000 reads; of these:
  1000000 (100.00%) were paired; of these:
    45123 (4.51%) aligned concordantly 0 times
    901234 (90.12%) aligned concordantly exactly 1 time
    53643 (5.36%) aligned concordantly >1 times
95.49% overall alignment rate
```

#### 4.2 SAM到BAM转换和排序

**执行命令：**
```bash
# 将SAM文件转换为BAM格式并排序
for sample in sample1 sample2 sample3 ctrl1 ctrl2 ctrl3; do
    echo "Converting $sample to BAM..."
    samtools view -bS results/alignment/${sample}.sam | \
    samtools sort -o results/alignment/${sample}_sorted.bam
    
    # 建立索引
    samtools index results/alignment/${sample}_sorted.bam
    
    # 删除SAM文件节省空间
    rm results/alignment/${sample}.sam
done
```

**结果验证：**
```bash
# 检查BAM文件
ls -lh results/alignment/*.bam
samtools flagstat results/alignment/sample1_sorted.bam
```

**检查点：** 确认所有样本的BAM文件已成功生成并建立索引。

---

### 步骤5：基因表达定量

#### 5.1 使用featureCounts进行基因计数

**操作说明：**
使用featureCounts根据基因注释文件对比对结果进行计数，生成基因表达矩阵。

**执行命令：**
```bash
# 使用featureCounts对所有样本进行基因计数
featureCounts -a reference/Homo_sapiens.GRCh38.104.gtf \
    -o results/quantification/gene_counts.txt \
    -T 4 \
    -p \
    -B \
    -C \
    -s 2 \
    results/alignment/*_sorted.bam
```

**参数解释：**
- `-a`：基因注释文件（GTF格式）
- `-o`：输出文件名
- `-T 4`：使用4个线程
- `-p`：双端测序模式
- `-B`：只计算正确配对的reads
- `-C`：不计算嵌合reads
- `-s 2`：链特异性设置（2表示反向链特异性）

**预期输出：**
```
        ==========     _____ _    _ ____  _____  ______          _____  
        =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       Rsubread 2.0.1

//========================== featureCounts setting ===========================\\
||                                                                            ||
||             Input files : 6 BAM files                                     ||
||                           sample1_sorted.bam                              ||
||                           sample2_sorted.bam                              ||
||                           ...                                              ||
||                                                                            ||
||             Output file : gene_counts.txt                                 ||
||                 Summary : gene_counts.txt.summary                         ||
||              Annotation : Homo_sapiens.GRCh38.104.gtf (GTF)              ||
||      Dir for temp files : results/quantification                          ||
||                                                                            ||
||                 Threads : 4                                               ||
||                   Level : meta-feature level                              ||
||              Paired-end : yes                                             ||
||      Multimapping reads : not counted                                     ||
|| Multi-overlapping reads : not counted                                     ||
||   Min overlapping bases : 1                                               ||
||                                                                            ||
\\============================================================================//

Process BAM file sample1_sorted.bam...
Process BAM file sample2_sorted.bam...
...

Summary of counting results can be found in file "results/quantification/gene_counts.txt.summary"
```

#### 5.2 检查计数结果

**执行命令：**
```bash
# 查看计数结果文件结构
head results/quantification/gene_counts.txt

# 查看计数统计摘要
cat results/quantification/gene_counts.txt.summary

# 统计检测到的基因数量
tail -n +2 results/quantification/gene_counts.txt | wc -l
```

**预期输出：**
```
# 基因计数矩阵前几行
Geneid  Chr     Start   End     Strand  Length  sample1_sorted.bam  sample2_sorted.bam  ...
ENSG00000000003 X       100627108       100639991       -       4535    156     189     ...
ENSG00000000005 X       100584802       100599885       +       1610    0       2       ...

# 统计摘要
Status  sample1_sorted.bam  sample2_sorted.bam  ...
Assigned        15234567        16789012        ...
Unassigned_Unmapped     1234567 1456789 ...
```

**检查点：** 确认基因计数矩阵已成功生成，大部分reads被成功分配到基因。

---

### 步骤6：差异表达分析

#### 6.1 准备R分析环境

**操作说明：**
创建R脚本进行DESeq2差异表达分析。

**执行命令：**
```bash
# 创建R分析脚本
cat > scripts/deseq2_analysis.R << 'EOF'
# 加载必要的包
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

# 设置工作目录
setwd("~/ngs-analysis/lesson-05")

# 读取计数矩阵
countData <- read.table("results/quantification/gene_counts.txt", 
                        header=TRUE, row.names=1, skip=1)

# 只保留计数列（去除注释列）
countData <- countData[,6:ncol(countData)]

# 修改列名（去除路径和扩展名）
colnames(countData) <- gsub("results.alignment.|_sorted.bam", "", colnames(countData))

# 创建样本信息表
sampleInfo <- data.frame(
    row.names = colnames(countData),
    condition = factor(c(rep("treatment", 3), rep("control", 3))),
    batch = factor(c(1,2,3,1,2,3))
)

print("Sample information:")
print(sampleInfo)

# 创建DESeq2对象
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = sampleInfo,
                              design = ~ condition)

# 过滤低表达基因
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

print(paste("Genes after filtering:", nrow(dds)))

# 运行差异分析
dds <- DESeq(dds)

# 获取结果
res <- results(dds, contrast=c("condition", "treatment", "control"))

# 查看结果摘要
summary(res)

# 保存结果
write.csv(as.data.frame(res), "results/deseq2/deseq2_results.csv")

# 获取显著差异基因
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(sig_genes), "results/deseq2/significant_genes.csv")

print(paste("Significant genes:", nrow(sig_genes)))

EOF
```

#### 6.2 运行DESeq2分析

**执行命令：**
```bash
# 创建DESeq2结果目录
mkdir -p results/deseq2

# 运行R脚本
Rscript scripts/deseq2_analysis.R
```

**预期输出：**
```
Sample information:
         condition batch
sample1  treatment     1
sample2  treatment     2
sample3  treatment     3
ctrl1    control       1
ctrl2    control       2
ctrl3    control       3

[1] "Genes after filtering: 18456"

out of 18456 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 1234, 6.7%
LFC < 0 (down)     : 987, 5.3%
outliers [1]       : 23, 0.12%
low counts [2]     : 0, 0%

[1] "Significant genes: 1456"
```

**检查点：** 确认差异分析成功完成，生成了结果文件。

---

### 步骤7：结果可视化

#### 7.1 创建可视化脚本

**执行命令：**
```bash
# 创建可视化脚本
cat > scripts/visualization.R << 'EOF'
# 加载包和数据
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

setwd("~/ngs-analysis/lesson-05")

# 重新加载DESeq2对象（或从上一步继续）
# 这里假设从上一步继续，如果重新运行需要重新加载数据

# 1. MA图
pdf("results/plots/MA_plot.pdf", width=8, height=6)
plotMA(res, main="MA Plot", ylim=c(-5,5))
dev.off()

# 2. 火山图
pdf("results/plots/volcano_plot.pdf", width=10, height=8)
EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue',
    title = 'Treatment vs Control',
    pCutoff = 0.05,
    FCcutoff = 1.0,
    pointSize = 2.0,
    labSize = 4.0)
dev.off()

# 3. 热图 - 显著差异基因
# 获取标准化后的表达数据
vsd <- vst(dds, blind=FALSE)

# 选择显著差异基因（前50个）
sig_genes_ordered <- sig_genes[order(sig_genes$padj),]
top_genes <- head(rownames(sig_genes_ordered), 50)

# 创建热图
pdf("results/plots/heatmap_top50.pdf", width=10, height=12)
pheatmap(assay(vsd)[top_genes,], 
         cluster_rows=TRUE, 
         show_rownames=TRUE, 
         cluster_cols=TRUE, 
         annotation_col=sampleInfo,
         scale="row",
         main="Top 50 Differentially Expressed Genes")
dev.off()

# 4. PCA图
pdf("results/plots/PCA_plot.pdf", width=8, height=6)
plotPCA(vsd, intgroup="condition") + 
    ggtitle("PCA Plot") +
    theme_minimal()
dev.off()

# 5. 样本距离热图
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, rownames(sampleInfo), sep="-")
colnames(sampleDistMatrix) <- NULL

pdf("results/plots/sample_distance_heatmap.pdf", width=8, height=6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         main="Sample Distance Heatmap")
dev.off()

print("All plots have been generated and saved to results/plots/")

EOF
```

#### 7.2 生成可视化图表

**执行命令：**
```bash
# 创建图表输出目录
mkdir -p results/plots

# 运行可视化脚本
Rscript scripts/visualization.R
```

**结果验证：**
```bash
# 检查生成的图表文件
ls -lh results/plots/
```

**检查点：** 确认所有可视化图表已成功生成。

---

### 步骤8：功能富集分析

#### 8.1 GO富集分析

**执行命令：**
```bash
# 创建功能分析脚本
cat > scripts/functional_analysis.R << 'EOF'
# 加载必要的包
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

setwd("~/ngs-analysis/lesson-05")

# 读取显著差异基因
sig_genes <- read.csv("results/deseq2/significant_genes.csv", row.names=1)

# 获取上调和下调基因
up_genes <- rownames(sig_genes[sig_genes$log2FoldChange > 0,])
down_genes <- rownames(sig_genes[sig_genes$log2FoldChange < 0,])

# 转换基因ID（从Ensembl ID到Gene Symbol）
up_symbols <- bitr(up_genes, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)
down_symbols <- bitr(down_genes, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)

# GO富集分析 - 上调基因
ego_up <- enrichGO(gene = up_symbols$SYMBOL,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

# GO富集分析 - 下调基因
ego_down <- enrichGO(gene = down_symbols$SYMBOL,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05,
                     readable = TRUE)

# 保存结果
write.csv(as.data.frame(ego_up), "results/deseq2/GO_enrichment_up.csv")
write.csv(as.data.frame(ego_down), "results/deseq2/GO_enrichment_down.csv")

# 可视化GO富集结果
pdf("results/plots/GO_enrichment_up.pdf", width=12, height=8)
dotplot(ego_up, showCategory=20) + ggtitle("GO Enrichment - Upregulated Genes")
dev.off()

pdf("results/plots/GO_enrichment_down.pdf", width=12, height=8)
dotplot(ego_down, showCategory=20) + ggtitle("GO Enrichment - Downregulated Genes")
dev.off()

# KEGG通路分析
kk_up <- enrichKEGG(gene = up_symbols$ENTREZID,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)

kk_down <- enrichKEGG(gene = down_symbols$ENTREZID,
                      organism = 'hsa',
                      pvalueCutoff = 0.05)

# 保存KEGG结果
write.csv(as.data.frame(kk_up), "results/deseq2/KEGG_enrichment_up.csv")
write.csv(as.data.frame(kk_down), "results/deseq2/KEGG_enrichment_down.csv")

print("Functional enrichment analysis completed!")
print(paste("Upregulated genes GO terms:", nrow(ego_up)))
print(paste("Downregulated genes GO terms:", nrow(ego_down)))

EOF
```

#### 8.2 运行功能分析

**执行命令：**
```bash
# 运行功能富集分析
Rscript scripts/functional_analysis.R
```

**预期输出：**
```
Functional enrichment analysis completed!
[1] "Upregulated genes GO terms: 45"
[1] "Downregulated genes GO terms: 38"
```

**检查点：** 确认功能富集分析完成，生成了GO和KEGG富集结果。

---

### 步骤9：结果整理和报告生成

#### 9.1 生成分析摘要

**执行命令：**
```bash
# 创建结果摘要脚本
cat > scripts/generate_summary.R << 'EOF'
# 生成分析结果摘要
setwd("~/ngs-analysis/lesson-05")

# 读取各种结果文件
deseq2_results <- read.csv("results/deseq2/deseq2_results.csv", row.names=1)
sig_genes <- read.csv("results/deseq2/significant_genes.csv", row.names=1)

# 创建摘要报告
summary_report <- list(
    "Total genes analyzed" = nrow(deseq2_results),
    "Significant genes (padj < 0.05, |FC| > 2)" = nrow(sig_genes),
    "Upregulated genes" = sum(sig_genes$log2FoldChange > 0),
    "Downregulated genes" = sum(sig_genes$log2FoldChange < 0),
    "Mean log2 fold change (significant genes)" = round(mean(abs(sig_genes$log2FoldChange)), 2),
    "Max upregulation (log2FC)" = round(max(sig_genes$log2FoldChange), 2),
    "Max downregulation (log2FC)" = round(min(sig_genes$log2FoldChange), 2)
)

# 保存摘要
capture.output(summary_report, file="results/analysis_summary.txt")

# 打印摘要
print("=== RNA-seq Analysis Summary ===")
for(i in 1:length(summary_report)) {
    cat(names(summary_report)[i], ":", summary_report[[i]], "\n")
}

EOF

# 运行摘要脚本
Rscript scripts/generate_summary.R
```

#### 9.2 整理最终结果

**执行命令：**
```bash
# 创建最终结果目录结构
mkdir -p final_results/{tables,figures,scripts}

# 复制重要结果文件
cp results/deseq2/deseq2_results.csv final_results/tables/
cp results/deseq2/significant_genes.csv final_results/tables/
cp results/deseq2/GO_enrichment_*.csv final_results/tables/
cp results/plots/*.pdf final_results/figures/
cp scripts/*.R final_results/scripts/
cp results/analysis_summary.txt final_results/

# 生成文件清单
ls -la final_results/ > final_results/file_list.txt
ls -la final_results/*/ >> final_results/file_list.txt

echo "Analysis completed! Results are organized in final_results/ directory."
```

## 预期结果

### 主要输出文件
1. **基因计数矩阵**：`results/quantification/gene_counts.txt`
   - 内容：每个基因在每个样本中的read计数
   - 用途：差异表达分析的输入数据

2. **差异表达结果**：`results/deseq2/deseq2_results.csv`
   - 内容：所有基因的差异表达统计结果
   - 用途：识别显著差异表达基因

3. **显著差异基因**：`results/deseq2/significant_genes.csv`
   - 内容：满足显著性阈值的差异表达基因
   - 用途：后续功能分析

4. **可视化图表**：`results/plots/`目录下的PDF文件
   - 内容：MA图、火山图、热图、PCA图等
   - 用途：结果展示和解释

5. **功能富集结果**：`results/deseq2/GO_enrichment_*.csv`
   - 内容：GO和KEGG富集分析结果
   - 用途：理解差异基因的生物学意义

### 关键结果指标
- **比对率**：应该在85-95%之间
- **基因检出数量**：预期检出15,000-20,000个基因
- **显著差异基因数量**：预期1,000-3,000个（取决于实验条件）
- **功能富集显著性**：p值 < 0.05的GO terms

### 成功标准
- [ ] 所有样本比对率 > 85%
- [ ] 基因计数矩阵成功生成
- [ ] DESeq2分析无错误完成
- [ ] 生成了所有预期的可视化图表
- [ ] 功能富集分析识别出生物学相关的通路

## 故障排除

### 常见问题1：内存不足错误
**症状：** 程序运行时出现"out of memory"错误
**原因：** 系统内存不足以处理大型基因组数据
**解决方案：**
```bash
# 减少并行线程数
hisat2 -x reference/grch38_index --threads 2  # 改为2线程
# 或者使用更小的测试数据集
```

### 常见问题2：R包安装失败
**症状：** 无法加载DESeq2或其他R包
**原因：** R包未正确安装或版本不兼容
**解决方案：**
```r
# 在R中重新安装包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
```

### 常见问题3：基因ID转换失败
**症状：** 功能富集分析时基因ID无法识别
**原因：** 基因ID格式不匹配或版本不一致
**解决方案：**
```r
# 检查基因ID格式
head(rownames(sig_genes))
# 使用biomaRt进行ID转换
library(biomaRt)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
```

### 获取帮助
如果遇到其他问题：
1. 检查错误日志：`cat logs/error.log`
2. 查看软件帮助：`hisat2 --help`, `featureCounts -h`
3. 联系助教或老师：wangys@hunau.edu.cn

## 扩展练习

### 练习1：参数优化比较
**目标：** 比较不同比对参数对结果的影响
**任务：** 使用不同的HISAT2参数重新比对一个样本，比较结果差异
**提示：** 尝试修改`--mp`、`--sp`等参数

### 练习2：批次效应分析
**目标：** 检测和校正批次效应
**任务：** 在DESeq2分析中加入批次因子，比较校正前后的结果
**提示：** 修改design公式为`~ batch + condition`

### 练习3：时间序列分析
**目标：** 分析基因表达的时间动态变化
**任务：** 如果有时间序列数据，使用maSigPro包进行分析
**提示：** 需要安装maSigPro包并准备时间点信息

### 思考问题
1. 为什么RNA-seq数据需要使用负二项分布模型而不是正态分布？
2. 多重检验校正的必要性是什么？FDR和Bonferroni校正有什么区别？
3. 如何解释MA图和火山图中的数据分布模式？

## 参考资料

### 相关文献
1. Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550 (2014).
2. Kim, D., Langmead, B. & Salzberg, S.L. HISAT: a fast spliced aligner with low memory requirements. Nature Methods 12, 357–360 (2015).
3. Liao, Y., Smyth, G.K. and Shi, W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics 30, 923-930 (2014).

### 在线资源
- DESeq2官方教程：https://bioconductor.org/packages/DESeq2/
- HISAT2用户手册：http://daehwankimlab.github.io/hisat2/
- clusterProfiler教程：https://yulab-smu.top/biomedical-knowledge-mining-book/

### 软件文档
- Bioconductor工作流：https://www.bioconductor.org/help/workflows/
- Galaxy RNA-seq教程：https://training.galaxyproject.org/training-material/topics/transcriptomics/

## 附录

### 附录A：完整分析脚本
参见：`final_results/scripts/` 目录中的R脚本文件

### 附录B：配置文件模板
```bash
# HISAT2配置示例
HISAT2_INDEX="reference/grch38_index"
THREADS=4
STRANDNESS="RF"

# DESeq2参数
PADJ_CUTOFF=0.05
FC_CUTOFF=1.0
```

### 附录C：数据格式说明
- **FASTQ格式**：原始测序数据格式，包含序列和质量信息
- **SAM/BAM格式**：序列比对结果格式，BAM是SAM的二进制压缩版本
- **GTF格式**：基因注释格式，包含基因结构信息
- **CSV格式**：差异分析结果格式，可用Excel打开

---

**实验完成时间：** 预计 4-6 小时  
**难度等级：** 中级  
**最后更新：** 2024年课程版本