# 第6次课：ChIP-seq与表观基因组分析

## 课程概述

本次课程将深入介绍ChIP-seq技术的原理、实验设计、数据分析流程以及表观基因组学的相关概念。学生将学习如何使用MACS2进行Peak calling，掌握ChIP-seq数据的质量控制方法，并了解表观遗传修饰的功能注释。

## 学习目标

- 理解ChIP-seq技术的基本原理和实验流程
- 掌握ChIP-seq数据分析的标准流程
- 学会使用MACS2进行Peak calling
- 了解表观基因组学的基本概念和分析方法
- 掌握ChIP-seq结果的可视化和功能注释

## 课程内容

### 理论部分（2学时）
- ChIP-seq技术原理
- 实验设计和质量控制
- Peak calling算法和工具
- 表观遗传修饰类型和功能
- 数据可视化和注释方法

### 实践部分（2学时）
- ChIP-seq数据预处理
- 使用MACS2进行Peak calling
- 结果质量评估和可视化
- Peak注释和功能富集分析

## 文件结构

```
lesson-06/
├── README.md                           # 本文件
├── DATA_SOURCES.md                     # 数据来源和准备指南（新）
├── slides/                             # 理论课幻灯片
│   ├── slides.md                      # Marp源文件
│   └── slides.html                    # 生成的HTML文件
├── manual/                             # 实践操作手册
│   ├── manual.md                      # 手册主文件
│   ├── data/                          # ChIP-seq数据目录
│   │   ├── README.md                  # 数据说明
│   │   ├── sample_info.txt            # 样本信息表
│   │   └── *.fastq.gz                 # ChIP-seq测序数据
│   ├── reference/                     # 参考基因组目录
│   │   ├── *.fa                       # 参考基因组FASTA
│   │   ├── *.gtf                      # 基因注释GTF
│   │   └── *.bwt                      # BWA索引文件
│   └── scripts/                       # 分析脚本（新）
│       ├── download_data.sh           # 自动化数据准备脚本
│       ├── chipseq_pipeline.sh        # 完整分析流水线脚本
│       └── peak_analysis.R            # Peak注释和可视化脚本
├── results/                            # 分析结果目录（运行时生成）
│   ├── fastqc/                        # 质量控制结果
│   ├── alignment/                     # 序列比对结果
│   ├── peaks/                         # Peak calling结果
│   ├── annotation/                    # Peak注释结果
│   └── plots/                         # 可视化图表
└── images/                             # 课程专用图片
    ├── chip_principle.svg             # ChIP-seq原理图
    ├── chipseq_workflow.svg           # ChIP-seq工作流程
    ├── macs2_algorithm.svg            # MACS2算法图
    ├── peak_calling.svg               # Peak calling流程
    └── epigenome_marks.svg            # 表观遗传标记
```

## 快速开始指南

### 1️⃣ 准备工作环境

```bash
# 进入实验目录
cd lesson-06/manual

# 创建必要的子目录
mkdir -p data reference results logs qc
```

### 2️⃣ 准备实验数据（推荐方式）

**自动化准备**（最简单）：
```bash
# 运行自动化脚本完成所有准备
bash scripts/download_data.sh

# 或指定准备特定部分
bash scripts/download_data.sh --reference    # 仅下载参考基因组
bash scripts/download_data.sh --sequencing   # 仅准备ChIP-seq数据
bash scripts/download_data.sh --validate     # 仅验证数据完整性
```

**详细信息请参考：** [DATA_SOURCES.md](DATA_SOURCES.md)

### 3️⃣ 检查软件环境

```bash
# 检查核心工具是否已安装
fastqc --version
bwa
samtools --version
macs2 --version
deeptools --version
```

### 4️⃣ 开始分析

```bash
# 运行完整分析流程
bash scripts/chipseq_pipeline.sh

# 或仅运行特定步骤
bash scripts/chipseq_pipeline.sh --qc         # 仅质量控制
bash scripts/chipseq_pipeline.sh --align      # 仅序列比对
bash scripts/chipseq_pipeline.sh --peaks      # 仅Peak calling
bash scripts/chipseq_pipeline.sh --visualize  # 仅生成可视化
```

## 软件要求

### 核心工具

| 软件名称 | 版本要求 | 安装方式 | 用途 |
|---------|---------|---------|------|
| BWA | >= 0.7.17 | `conda install bwa` | 序列比对 |
| SAMtools | >= 1.10 | `conda install samtools` | BAM文件处理 |
| MACS2 | >= 2.2.0 | `pip install MACS2` | Peak calling |
| FastQC | >= 0.11.9 | `conda install fastqc` | 质量控制 |
| deepTools | >= 3.5.0 | `pip install deepTools` | 可视化 |
| R | >= 4.0 | `conda install r-base` | 统计分析 |

### 可选工具

- **ChIPseeker** - R包，用于Peak注释
- **IGV** - 基因组浏览器，用于可视化
- **bedtools** - 基因组工具集

### 安装示例

```bash
# 使用Conda安装主要工具（推荐）
conda install -c bioconda bwa samtools macs2 fastqc deeptools

# 使用pip安装MACS2（备选）
pip install MACS2 deepTools
```

## 数据要求

### 参考数据

| 数据类型 | 版本 | 大小 | 来源 |
|---------|------|------|------|
| 人类基因组 | GRCh38 (hg38) | 完整版~3GB / 演示版~5MB | Ensembl / UCSC |
| 基因注释 | Ensembl v110 | ~50MB | ftp://ftp.ensembl.org |

### ChIP-seq数据

- **样本类型**: ChIP样本 + Input对照
- **测序类型**: 双端配对末端 (PE, 2×50bp)
- **数据量**:
  - 演示数据：每样本约500k reads
  - 完整数据：每样本约10-50M reads
- **总大小**:
  - 演示数据：~100MB
  - 完整数据：~5-10GB

### 数据来源选项

✅ **演示数据** - 快速学习（推荐初学者）
- 数据量：500k reads/样本
- 时间：5-10分钟准备

📊 **模拟数据** - 完整学习（推荐中级）
- 数据量：100万 reads/样本
- 时间：15-20分钟生成

🔬 **真实数据** - 深度学习（推荐高级）
- 来源：ENCODE / GEO数据库
- 时间：30分钟-2小时下载

详见 [DATA_SOURCES.md](DATA_SOURCES.md) 获取完整说明

## 预计时间

- 理论讲授：2小时
- 实践操作：2小时
- 总计：4学时