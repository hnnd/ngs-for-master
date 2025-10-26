# 第6次课 - ChIP-seq数据来源与准备指南

## 概述

本课程提供了多种获取ChIP-seq数据的方式，满足不同的学习需求：
- **演示数据** - 快速学习，了解分析流程
- **模拟数据** - 完整的生物学模式，适合参数优化学习
- **真实数据** - 来自公共数据库的实际研究数据

## 数据准备方式

### 方式1：自动化脚本准备（推荐）

最简单的方式是使用提供的自动化脚本：

```bash
# 进入manual目录
cd lesson-06/manual

# 使用脚本准备所有数据
bash scripts/download_data.sh

# 或指定准备特定数据
bash scripts/download_data.sh --reference    # 仅下载参考基因组
bash scripts/download_data.sh --sequencing   # 仅准备ChIP-seq数据
bash scripts/download_data.sh --validate     # 仅验证数据完整性
```

### 方式2：手动准备

#### 步骤1：准备参考基因组

**选项A：从Ensembl下载完整基因组**

```bash
cd reference

# 下载完整人类基因组
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

**选项B：仅下载演示数据（22号染色体片段）**

```bash
cd reference

# 从Ensembl下载22号染色体
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
gunzip Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz

# 提取前5Mb作为演示数据
samtools faidx Homo_sapiens.GRCh38.dna.chromosome.22.fa
samtools faidx Homo_sapiens.GRCh38.dna.chromosome.22.fa 22:1-5000000 > hg38_chr22_demo.fa
```

#### 步骤2：下载基因注释文件

```bash
cd reference

# 下载GTF注释文件
wget ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip Homo_sapiens.GRCh38.110.gtf.gz
```

#### 步骤3：准备ChIP-seq测序数据

**选项A：从NCBI SRA下载真实数据**

```bash
cd data

# 安装sra-toolkit（如果未安装）
conda install -c bioconda sra-tools

# 下载ChIP-seq数据集 (以H3K4me3为例)
# 数据来源: ENCODE项目或GEO数据库

# 示例：下载来自ENCODE的H3K4me3数据
fastq-dump --split-files --gzip -X 500000 SRR037648 -O .
```

**选项B：生成模拟数据（推荐用于快速学习）**

```bash
cd data

# 安装wgsim工具
conda install -c bioconda wgsim

# 生成ChIP样本数据
wgsim -N 500000 -1 50 -2 50 -r 0.001 -R 0.15 -X 0.3 \
    ../reference/hg38_chr22_demo.fa \
    H3K4me3_ChIP_R1.fastq \
    H3K4me3_ChIP_R2.fastq

# 生成Input对照数据
wgsim -N 500000 -1 50 -2 50 -r 0.001 -R 0.15 -X 0.3 \
    ../reference/hg38_chr22_demo.fa \
    Input_control_R1.fastq \
    Input_control_R2.fastq

# 压缩文件
gzip *.fastq
```

## 数据详细说明

### 参考基因组

| 特征 | 详情 |
|------|------|
| **物种** | 人类 (*Homo sapiens*) |
| **基因组版本** | GRCh38 (hg38) |
| **完整基因组大小** | ~3.1GB |
| **演示数据** | 22号染色体的5Mb片段（~5MB） |
| **来源** | Ensembl / UCSC基因组浏览器 |
| **格式** | FASTA |

### ChIP-seq数据

#### 演示数据特征

```
样本设计:
  - ChIP样本: H3K4me3 (三甲基化组蛋白H3K4)
  - 对照样本: Input control (全基因组DNA)

序列特征:
  - 测序类型: Paired-end (PE) 双端测序
  - 读长: 50bp × 50bp (R1 × R2)
  - 测序深度: ~500k reads/样本
  - 预期Peak数: 1000-5000个

数据质量:
  - 平均Phred质量分数: Q30 (~99.9%准确率)
  - GC含量: ~40-45% (人类基因组典型值)
```

### 样本信息

```
ChIP样本 - H3K4me3:
  - 标记: 三甲基化H3K4
  - 位置: 活跃启动子和转录开始位点
  - 预期特点: 尖锐、定位准确的peaks

Input对照:
  - 类型: 未经ChIP的全基因组DNA
  - 用途: 用作背景校正和质量评估
  - 关键: 必须与ChIP样本来自同一生物样本
```

## 数据验证

### 完整性检查

```bash
# 1. 检查文件大小
ls -lh reference/
ls -lh data/*.fastq.gz

# 2. 验证FASTA格式
head -10 reference/*.fa

# 3. 验证FASTQ格式
zcat data/H3K4me3_ChIP_R1.fastq.gz | head -8

# 4. 检查PE数据配对
r1_reads=$(zcat data/H3K4me3_ChIP_R1.fastq.gz | wc -l | awk '{print $1/4}')
r2_reads=$(zcat data/H3K4me3_ChIP_R2.fastq.gz | wc -l | awk '{print $1/4}')
if [ $r1_reads -eq $r2_reads ]; then
    echo "✓ PE数据配对正确"
fi
```

## 磁盘空间需求

| 数据类型 | 大小 | 说明 |
|---------|------|------|
| 参考基因组完整版 | ~3GB | 可选 |
| 参考基因组演示版 | ~5MB | 推荐 |
| 基因注释 | ~50MB | 必需 |
| ChIP-seq数据（演示） | ~100MB | 快速学习 |
| ChIP-seq数据（完整） | ~5-10GB | 完整分析 |
| 分析结果 | ~500MB | BAM和中间文件 |
| **总计（推荐）** | **~10GB** | 演示+完整分析 |

## 常见问题

### Q1: 如何快速开始？

使用演示数据：
```bash
bash scripts/download_data.sh --all
bash ../scripts/chipseq_pipeline.sh --all
```

### Q2: 可以使用自己的数据吗？

完全可以！只需确保：
1. FASTQ格式的ChIP-seq配对末端数据
2. 相应物种的参考基因组FASTA
3. 修改脚本中的文件路径和参数

### Q3: MACS2的参数如何调整？

主要参数说明：
```bash
-q          # q值阈值（默认0.05）
-m          # fold enrichment最小值（默认10）
--keep-dup  # 重复序列处理方式
--nomodel   # 是否模拟片段长度
```

## 获取更多帮助

- ChIP-seq最佳实践：https://www.encodeproject.org/
- 基因组浏览器：https://genome.ucsc.edu
- GEO数据库：https://www.ncbi.nlm.nih.gov/geo

---

**最后更新**: 2024年版本
**推荐数据版本**: GRCh38 Release 110
