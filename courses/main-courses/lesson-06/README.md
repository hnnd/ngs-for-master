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
├── README.md                 # 本文件
├── slides/                   # 理论课幻灯片
│   ├── slides.md            # Marp源文件
│   └── slides.html          # 生成的HTML文件
├── manual/                   # 实践操作手册
│   ├── manual.md            # 手册主文件
│   ├── scripts/             # 分析脚本
│   └── data/                # 示例数据
└── images/                   # 课程专用图片
    ├── chip_principle.svg   # ChIP-seq原理图
    ├── peak_calling.svg     # Peak calling流程
    └── epigenome_marks.svg  # 表观遗传标记
```

## 软件要求

- MACS2 (>= 2.2.0)
- deepTools (>= 3.5.0)
- ChIPseeker (R包)
- IGV (Integrative Genomics Viewer)
- Python 3.8+
- R 4.0+

## 数据要求

- ChIP-seq FASTQ文件（H3K4me3, H3K27ac示例）
- 对应的Input对照文件
- 参考基因组文件（hg38）
- 基因注释文件（GTF格式）

## 预计时间

- 理论讲授：2小时
- 实践操作：2小时
- 总计：4学时