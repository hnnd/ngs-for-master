# 第6次课实践操作手册

## 课程信息
- **课程名称**：高通量测序数据分析
- **主讲教师**：王运生
- **联系邮箱**：wangys@hunau.edu.cn
- **办公室**：16教420室
- **上课地点**：105机房
- **课程时间**：第6次课（4学时）

## 实验目标

### 主要目标
- 掌握ChIP-seq数据分析的完整流程
- 学会使用MACS2进行Peak calling
- 了解ChIP-seq数据质量控制方法
- 掌握Peak注释和功能富集分析

### 预期成果
- 完成H3K4me3 ChIP-seq数据的完整分析
- 生成高质量的Peak calling结果
- 掌握ChIP-seq结果的可视化方法
- 理解表观遗传修饰的生物学意义

## 环境要求

### 软件环境
| 软件名称 | 版本要求 | 安装方式 | 说明 |
|---------|---------|---------|------|
| MACS2 | >= 2.2.0 | `pip install MACS2` | Peak calling工具 |
| deepTools | >= 3.5.0 | `pip install deepTools` | 可视化工具 |
| SAMtools | >= 1.10 | `conda install samtools` | BAM文件处理 |
| BWA | >= 0.7.17 | `conda install bwa` | 序列比对 |
| FastQC | >= 0.11.9 | `conda install fastqc` | 质量控制 |
| R | >= 4.0 | `conda install r-base` | 统计分析 |
| ChIPseeker | latest | R包安装 | Peak注释 |
| IGV | >= 2.8 | 官网下载 | 基因组浏览器 |

### 硬件要求
- **内存**：至少 8 GB RAM
- **存储空间**：至少 20 GB 可用空间
- **CPU**：多核处理器推荐
- **网络**：稳定的互联网连接（用于下载数据）

### 数据准备

#### 📥 推荐方式：使用自动化脚本

```bash
# 进入manual目录
cd lesson-06/manual

# 自动准备所有数据
bash scripts/download_data.sh

# 该脚本将自动完成：
# ✓ 检查必需的软件
# ✓ 生成参考基因组（演示版本）
# ✓ 生成ChIP-seq模拟数据
# ✓ 创建样本信息表
# ✓ 验证数据完整性
```

#### 📊 数据选项

| 数据类型 | 大小 | 准备时间 | 说明 |
|---------|------|---------|------|
| 演示数据 | ~100MB | 5-10分钟 | 快速学习，推荐初学者 |
| 模拟数据 | ~500MB | 15-20分钟 | 完整分析，使用wgsim生成 |
| 真实数据 | ~5-10GB | 30分钟-2小时 | ENCODE/GEO数据集 |

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
mkdir -p ~/ngs-analysis/lesson-06
cd ~/ngs-analysis/lesson-06

# 创建子目录结构
mkdir -p {data,scripts,results,logs,figures}
mkdir -p results/{fastqc,alignment,peaks,annotation}
```

#### 1.2 检查软件环境
```bash
# 检查必要软件是否已安装
macs2 --version
samtools --version
bwa
fastqc --version
deeptools --version
```

**预期输出：**
```
macs2 2.2.7.1
samtools 1.15.1
bwa 0.7.17-r1188
FastQC v0.11.9
deepTools 3.5.1
```

#### 1.3 使用脚本自动准备数据（推荐）

```bash
# 进入manual目录
cd ~/ngs-analysis/lesson-06/manual

# 运行自动化数据准备脚本
bash scripts/download_data.sh --all

# 脚本将自动：
# ✓ 检查必需的软件工具
# ✓ 生成参考基因组（5Mb演示版本）
# ✓ 生成ChIP-seq模拟数据
# ✓ 创建样本信息表
# ✓ 验证数据完整性
```

**检查点：** 查看日志文件确认数据准备成功

```bash
cat ../logs/data_summary.txt
ls -lh reference/
ls -lh data/
```

#### 1.4 手动准备数据（可选）

如果自动脚本无法运行或需要特定数据，可参考详细指南：

- **完整操作步骤**: 见 [../DATA_SOURCES.md](../DATA_SOURCES.md)
- **参考基因组下载**: 选择Ensembl或UCSC源
- **测序数据获取**:
  - 直接下载ENCODE/GEO数据
  - 使用wgsim生成模拟数据
  - 使用课程服务器的预置数据

---

### 步骤2：数据质量控制

#### 2.1 FastQC质量评估

**操作说明：**
使用FastQC评估原始测序数据的质量，检查序列质量分布、GC含量、重复序列等指标。

**执行命令：**
```bash
# 对ChIP和Input数据进行质量评估
cd ~/ngs-analysis/lesson-06
fastqc data/H3K4me3_ChIP.fastq.gz -o results/fastqc/
fastqc data/Input_control.fastq.gz -o results/fastqc/

# 查看结果
ls results/fastqc/
```

**参数解释：**
- `-o`：指定输出目录
- `fastq.gz`：可以直接处理压缩文件

**预期输出：**
```
H3K4me3_ChIP_fastqc.html
H3K4me3_ChIP_fastqc.zip
Input_control_fastqc.html
Input_control_fastqc.zip
```

**结果验证：**
```bash
# 在浏览器中查看HTML报告
firefox results/fastqc/H3K4me3_ChIP_fastqc.html &
```

#### 2.2 数据预处理（可选）

如果数据质量较差，可以进行质量过滤和接头去除：

```bash
# 使用Trimmomatic进行数据清洗（如果需要）
trimmomatic SE -phred33 \
    data/H3K4me3_ChIP.fastq.gz \
    results/H3K4me3_ChIP_trimmed.fastq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

**检查点：** FastQC报告显示数据质量良好，Per base sequence quality大部分区域在绿色范围内。

---

### 步骤3：序列比对

#### 3.1 构建参考基因组索引

**操作说明：**
使用BWA为参考基因组构建索引，这是进行序列比对的前提步骤。

**执行命令：**
```bash
# 构建BWA索引（如果尚未构建）
cd data
bwa index hg38.fa

# 同时构建samtools索引
samtools faidx hg38.fa
```

**参数解释：**
- `index`：BWA索引构建命令
- `faidx`：SAMtools FASTA索引

**预期输出：**
```
hg38.fa.amb
hg38.fa.ann
hg38.fa.bwt
hg38.fa.pac
hg38.fa.sa
hg38.fa.fai
```

#### 3.2 序列比对

**执行命令：**
```bash
# 比对ChIP-seq数据
cd ~/ngs-analysis/lesson-06
bwa mem -t 4 data/hg38.fa data/H3K4me3_ChIP.fastq.gz | \
    samtools sort -@ 4 -o results/alignment/H3K4me3_ChIP.bam

# 比对Input对照数据
bwa mem -t 4 data/hg38.fa data/Input_control.fastq.gz | \
    samtools sort -@ 4 -o results/alignment/Input_control.bam

# 建立BAM索引
samtools index results/alignment/H3K4me3_ChIP.bam
samtools index results/alignment/Input_control.bam
```

**参数解释：**
- `-t 4`：使用4个线程
- `-@ 4`：SAMtools使用4个线程
- `sort`：按坐标排序
- `index`：建立BAM索引

**结果验证：**
```bash
# 检查比对统计
samtools flagstat results/alignment/H3K4me3_ChIP.bam
samtools flagstat results/alignment/Input_control.bam
```

**检查点：** 比对率应该在80%以上，重复率在合理范围内（<30%）。

---

### 步骤4：Peak Calling

#### 4.1 使用MACS2进行Peak calling

**操作说明：**
MACS2是最常用的ChIP-seq Peak calling工具，它使用泊松分布模型识别显著富集的区域。

**执行命令：**
```bash
# 使用MACS2进行Peak calling
cd ~/ngs-analysis/lesson-06
macs2 callpeak \
    -t results/alignment/H3K4me3_ChIP.bam \
    -c results/alignment/Input_control.bam \
    -f BAM \
    -g hs \
    -n H3K4me3 \
    -q 0.05 \
    --outdir results/peaks/ \
    --call-summits
```

**参数解释：**
- `-t`：处理组（ChIP）BAM文件
- `-c`：对照组（Input）BAM文件
- `-f BAM`：输入文件格式
- `-g hs`：基因组大小（hs=human，2.7e9）
- `-n`：输出文件前缀
- `-q 0.05`：FDR阈值
- `--call-summits`：检测峰顶位置

**预期输出：**
```
H3K4me3_peaks.narrowPeak
H3K4me3_peaks.xls
H3K4me3_summits.bed
H3K4me3_treat_pileup.bdg
H3K4me3_control_lambda.bdg
```

#### 4.2 Peak结果统计

**执行命令：**
```bash
# 统计Peak数量
wc -l results/peaks/H3K4me3_peaks.narrowPeak

# 查看Peak质量分布
head -20 results/peaks/H3K4me3_peaks.narrowPeak

# 计算FRiP score
total_reads=$(samtools view -c results/alignment/H3K4me3_ChIP.bam)
peak_reads=$(bedtools intersect -a results/alignment/H3K4me3_ChIP.bam \
    -b results/peaks/H3K4me3_peaks.narrowPeak -u | samtools view -c)
frip=$(echo "scale=4; $peak_reads / $total_reads" | bc)
echo "FRiP Score: $frip"
```

**检查点：** Peak数量应该在10,000-50,000之间，FRiP score应该>5%。

---

### 步骤5：数据可视化

#### 5.1 生成BigWig文件

**操作说明：**
BigWig文件用于在基因组浏览器中可视化ChIP-seq信号。

**执行命令：**
```bash
# 使用deepTools生成标准化的BigWig文件
bamCoverage -b results/alignment/H3K4me3_ChIP.bam \
    -o results/H3K4me3_ChIP.bw \
    --normalizeUsing RPKM \
    --binSize 10 \
    --numberOfProcessors 4

bamCoverage -b results/alignment/Input_control.bam \
    -o results/Input_control.bw \
    --normalizeUsing RPKM \
    --binSize 10 \
    --numberOfProcessors 4
```

**参数解释：**
- `--normalizeUsing RPKM`：使用RPKM标准化
- `--binSize 10`：10bp窗口大小
- `--numberOfProcessors 4`：使用4个处理器

#### 5.2 生成Peak区域热图

**执行命令：**
```bash
# 计算Peak区域的信号矩阵
computeMatrix reference-point \
    -S results/H3K4me3_ChIP.bw \
    -R results/peaks/H3K4me3_summits.bed \
    --referencePoint center \
    -b 2000 -a 2000 \
    --skipZeros \
    -o results/H3K4me3_matrix.gz \
    --numberOfProcessors 4

# 生成热图
plotHeatmap -m results/H3K4me3_matrix.gz \
    -out figures/H3K4me3_heatmap.png \
    --colorMap Blues \
    --whatToShow 'heatmap and colorbar'
```

#### 5.3 生成Profile图

**执行命令：**
```bash
# 生成Profile图
plotProfile -m results/H3K4me3_matrix.gz \
    -out figures/H3K4me3_profile.png \
    --numPlotsPerRow 1 \
    --plotTitle "H3K4me3 Signal Profile"
```

**检查点：** 热图和Profile图显示H3K4me3信号在Peak中心区域富集。

---

### 步骤6：Peak注释和功能分析

#### 6.1 安装和加载R包

**执行命令：**
```bash
# 启动R环境
R
```

在R中执行：
```r
# 安装必要的R包（如果尚未安装）
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", 
                       "clusterProfiler", "org.Hs.eg.db"))

# 加载包
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
```

#### 6.2 Peak注释

```r
# 读取Peak文件
peaks <- readPeakFile("results/peaks/H3K4me3_peaks.narrowPeak")

# 获取基因组注释
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Peak注释
peakAnno <- annotatePeak(peaks, tssRegion=c(-2000, 2000), TxDb=txdb)

# 查看注释结果
print(peakAnno)

# 可视化Peak基因组分布
pdf("figures/peak_annotation.pdf", width=10, height=6)
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
plotDistToTSS(peakAnno)
dev.off()
```

#### 6.3 功能富集分析

```r
# 获取Peak相关基因
genes <- as.data.frame(peakAnno)$geneId
genes <- genes[!is.na(genes)]

# GO富集分析
ego <- enrichGO(gene = genes,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable = TRUE)

# 可视化GO结果
pdf("figures/GO_enrichment.pdf", width=12, height=8)
dotplot(ego, showCategory=20)
barplot(ego, showCategory=20)
dev.off()

# KEGG通路分析
kk <- enrichKEGG(gene = genes,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)

pdf("figures/KEGG_enrichment.pdf", width=12, height=8)
dotplot(kk, showCategory=20)
dev.off()

# 保存结果
write.csv(as.data.frame(ego), "results/annotation/GO_enrichment.csv")
write.csv(as.data.frame(kk), "results/annotation/KEGG_enrichment.csv")

# 退出R
quit(save="no")
```

**检查点：** GO分析应该显示与转录调控、染色质修饰相关的生物过程富集。

---

### 步骤7：IGV可视化

#### 7.1 准备IGV文件

**执行命令：**
```bash
# 转换Peak文件为IGV兼容格式
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"Peak_"NR,$5}' \
    results/peaks/H3K4me3_peaks.narrowPeak > results/H3K4me3_peaks.bed

# 创建IGV session文件
cat > results/H3K4me3_session.xml << 'EOF'
<?xml version="1.0" encoding="UTF-8"?>
<Session genome="hg38" version="8">
    <Resources>
        <Resource path="results/H3K4me3_ChIP.bw"/>
        <Resource path="results/Input_control.bw"/>
        <Resource path="results/H3K4me3_peaks.bed"/>
    </Resources>
</Session>
EOF
```

#### 7.2 IGV可视化操作

1. 启动IGV：`igv`
2. 加载人类基因组hg38
3. 加载BigWig文件和Peak文件
4. 导航到感兴趣的基因区域（如：chr1:1,000,000-1,100,000）
5. 调整轨道高度和颜色
6. 截图保存重要区域

**检查点：** 在IGV中可以看到H3K4me3信号在启动子区域的富集，与Peak区域一致。

---

### 步骤8：结果整理和报告

#### 8.1 生成分析报告

**执行命令：**
```bash
# 创建分析总结
cat > results/analysis_summary.txt << 'EOF'
ChIP-seq Analysis Summary
========================

Sample: H3K4me3 ChIP-seq
Date: $(date)

Data Quality:
- Total reads: $(samtools view -c results/alignment/H3K4me3_ChIP.bam)
- Mapped reads: $(samtools view -c -F 4 results/alignment/H3K4me3_ChIP.bam)
- Mapping rate: $(samtools flagstat results/alignment/H3K4me3_ChIP.bam | grep "mapped (" | cut -d'(' -f2 | cut -d' ' -f1)

Peak Calling Results:
- Total peaks: $(wc -l < results/peaks/H3K4me3_peaks.narrowPeak)
- FRiP Score: $frip

Peak Annotation:
- See results/annotation/ for detailed annotation results
- GO and KEGG enrichment analysis completed

Files Generated:
- BAM files: results/alignment/
- Peak files: results/peaks/
- Visualization: figures/
- Annotation: results/annotation/
EOF
```

#### 8.2 文件整理

```bash
# 整理重要结果文件
mkdir -p final_results
cp results/peaks/H3K4me3_peaks.narrowPeak final_results/
cp results/H3K4me3_ChIP.bw final_results/
cp figures/*.png final_results/
cp results/annotation/*.csv final_results/

# 压缩结果
tar -czf ChIP-seq_analysis_results.tar.gz final_results/
```

## 预期结果

### 主要输出文件
1. **Peak文件**：`results/peaks/H3K4me3_peaks.narrowPeak`
   - 内容：Peak区域坐标和统计信息
   - 用途：下游分析和可视化

2. **BigWig文件**：`results/H3K4me3_ChIP.bw`
   - 内容：标准化的信号轨迹
   - 用途：基因组浏览器可视化

3. **注释结果**：`results/annotation/GO_enrichment.csv`
   - 内容：功能富集分析结果
   - 用途：生物学解释

### 关键结果指标
- **Peak数量**：应该在10,000-50,000之间
- **FRiP Score**：预期值>5%（H3K4me3通常>10%）
- **Peak分布**：主要位于启动子区域（~60%）

### 成功标准
- [ ] 所有命令执行无错误
- [ ] 生成了预期的输出文件
- [ ] FRiP score在合理范围内
- [ ] Peak主要富集在启动子区域
- [ ] GO分析显示转录相关功能富集

## 故障排除

### 常见问题1：内存不足错误
**症状：** MACS2或其他工具报告内存不足
**原因：** 数据文件过大，系统内存不够
**解决方案：**
```bash
# 减少并行线程数
macs2 callpeak ... --buffer-size 10000
# 或者使用更严格的过滤参数
```

### 常见问题2：Peak数量异常
**症状：** Peak数量过多（>100,000）或过少（<1,000）
**原因：** 参数设置不当或数据质量问题
**解决方案：**
```bash
# 调整q-value阈值
macs2 callpeak ... -q 0.01  # 更严格
macs2 callpeak ... -q 0.1   # 更宽松

# 检查数据质量
samtools flagstat results/alignment/H3K4me3_ChIP.bam
```

### 常见问题3：R包安装失败
**症状：** BiocManager安装包时出错
**原因：** 网络问题或依赖缺失
**解决方案：**
```r
# 更换镜像源
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# 或者手动安装依赖
```

### 获取帮助
如果遇到其他问题：
1. 检查错误日志：`cat logs/error.log`
2. 查看软件帮助：`macs2 callpeak --help`
3. 联系助教或老师：wangys@hunau.edu.cn

## 扩展练习

### 练习1：不同参数比较
**目标：** 比较不同MACS2参数对结果的影响
**任务：** 使用不同的q-value（0.01, 0.05, 0.1）重复Peak calling，比较结果差异
**提示：** 使用不同的输出前缀区分结果

### 练习2：多样本比较
**目标：** 分析不同组蛋白修饰的差异
**任务：** 如果有H3K27ac数据，比较其与H3K4me3的Peak分布差异
**提示：** 使用bedtools intersect分析重叠情况

### 思考问题
1. 为什么H3K4me3主要富集在启动子区域？
2. FRiP score的高低说明了什么？
3. 如何解释GO富集分析的结果？

## 参考资料

### 相关文献
1. Zhang Y, et al. Model-based analysis of ChIP-Seq (MACS). Genome Biol. 2008.
2. Landt SG, et al. ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia. Genome Res. 2012.

### 在线资源
- MACS2文档：https://github.com/macs3-project/MACS
- deepTools文档：https://deeptools.readthedocs.io/
- ChIPseeker教程：https://bioconductor.org/packages/ChIPseeker/

### 软件文档
- IGV用户指南：https://software.broadinstitute.org/software/igv/
- ENCODE ChIP-seq标准：https://www.encodeproject.org/chip-seq/

## 附录

### 附录A：完整脚本文件
参见：`scripts/` 目录中的相关脚本文件

### 附录B：配置文件模板
```bash
# MACS2配置文件示例
[DEFAULT]
gsize = hs
qvalue = 0.05
format = BAM
```

### 附录C：数据格式说明

**narrowPeak格式说明：**
1. chrom - 染色体名称
2. chromStart - Peak起始位置
3. chromEnd - Peak结束位置
4. name - Peak名称
5. score - Peak得分
6. strand - 链方向（通常为"."）
7. signalValue - 信号值
8. pValue - P值
9. qValue - Q值（FDR校正后）
10. peak - Peak顶点相对位置

---

**实验完成时间：** 预计 4 小时  
**难度等级：** 中级  
**最后更新：** 2025年