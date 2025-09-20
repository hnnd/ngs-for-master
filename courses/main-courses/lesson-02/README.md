# 第2次课：测序数据质量控制与预处理

## 课程信息
- **课程名称**：高通量测序数据分析
- **主讲教师**：王运生
- **联系邮箱**：wangys@hunau.edu.cn
- **办公室**：16教420室
- **上课地点**：105机房
- **课程时长**：4学时（理论2学时 + 实践2学时）

## 课程概述

本次课程重点介绍高通量测序数据的质量控制与预处理技术，包括测序错误的来源、质量评估方法、数据清洗策略等理论知识，以及FastQC、MultiQC等工具的实际操作。

## 学习目标

### 理论目标
- 理解测序数据中各种错误的来源和特征
- 掌握测序数据质量评估的主要指标和方法
- 了解不同数据清洗策略的原理和适用场景
- 熟悉质量控制在测序数据分析流程中的重要性

### 实践目标
- 熟练使用FastQC进行单个样本的质量评估
- 掌握MultiQC进行批量样本质量报告的生成
- 学会使用Trimmomatic等工具进行数据清洗
- 能够解读质量控制报告并制定相应的处理策略

## 课程内容

### 理论部分（2学时，约50张幻灯片）
1. **测序错误来源分析**
   - 测序化学反应中的错误
   - 光学系统引入的噪声
   - 碱基识别算法的局限性
   - 样本制备过程中的问题

2. **质量评估指标体系**
   - Phred质量分数系统
   - 每碱基质量分布
   - 每序列质量分布
   - GC含量分析
   - 序列长度分布
   - 重复序列检测

3. **数据清洗策略**
   - 低质量碱基修剪
   - 接头序列去除
   - 重复序列处理
   - 污染序列过滤
   - 长度过滤策略

4. **质量控制工具比较**
   - FastQC功能特点
   - MultiQC整合分析
   - 其他QC工具介绍
   - 工具选择策略

### 实践部分（2学时）
1. **环境准备和数据下载**
   - 软件安装和配置
   - 测试数据集准备
   - 工作目录组织

2. **FastQC质量评估**
   - 单样本质量分析
   - 报告解读和问题识别
   - 批量处理脚本编写

3. **MultiQC整合分析**
   - 多样本报告生成
   - 质量趋势分析
   - 异常样本识别

4. **数据清洗实战**
   - Trimmomatic参数设置
   - 清洗效果评估
   - 前后质量对比

## 文件结构

```
lesson-02/
├── README.md                    # 本文件
├── slides/                     # 理论幻灯片
│   └── slides.md              # Marp格式幻灯片源文件
├── manual/                    # 实践操作手册
│   ├── manual.md              # 操作手册主文件
│   ├── scripts/               # 相关脚本
│   │   ├── setup.sh          # 环境设置脚本
│   │   ├── download_data.sh   # 数据下载脚本
│   │   ├── run_fastqc.sh     # FastQC批量运行脚本
│   │   ├── run_multiqc.py    # MultiQC分析脚本
│   │   ├── quality_trim.sh   # 数据清洗脚本
│   │   └── compare_quality.py # 质量对比脚本
│   └── data/                  # 示例数据
│       ├── raw_reads/         # 原始测序数据
│       ├── fastqc_results/    # FastQC结果
│       ├── multiqc_results/   # MultiQC结果
│       └── trimmed_reads/     # 清洗后数据
└── images/                    # 课程专用图片
    ├── error_sources.svg      # 错误来源示意图
    ├── quality_metrics.svg    # 质量指标图表
    ├── qc_workflow.svg        # 质量控制流程图
    ├── fastqc_interface.svg   # FastQC界面截图
    ├── multiqc_dashboard.svg  # MultiQC仪表板
    └── before_after_qc.svg    # 质量控制前后对比
```

## 重要概念

### 关键术语
- **Phred Score**：碱基质量分数，表示碱基识别的准确性
- **Per-base Quality**：每个位置的质量分布
- **GC Content**：鸟嘌呤和胞嘧啶的含量百分比
- **Adapter Contamination**：接头序列污染
- **Overrepresented Sequences**：过度表达的序列

### 质量标准
- **高质量数据**：Q30以上碱基占比 > 85%
- **可接受数据**：Q20以上碱基占比 > 90%
- **需要清洗**：存在明显的质量下降趋势或接头污染

## 预备知识

学生应具备以下基础知识：
- Linux基本命令操作
- 测序技术基本原理
- FASTQ文件格式理解
- 基本的统计概念

## 后续课程

本次课程为后续分析奠定基础：
- 第3次课：序列比对算法与工具
- 第4次课：变异检测与基因分型
- 数据质量直接影响后续分析的准确性

## 评估方式

- **课堂参与**：积极参与讨论和提问
- **实践操作**：完成所有实验步骤
- **作业提交**：提交质量控制报告和分析结果
- **理解程度**：能够解释质量控制的重要性和方法选择

## 参考资料

### 主要文献
1. Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data.
2. Ewels, P., et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048.
3. Bolger, A.M., et al. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114-2120.

### 在线资源
- FastQC官方文档：https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- MultiQC官方文档：https://multiqc.info/
- Trimmomatic用户手册：http://www.usadellab.org/cms/?page=trimmomatic

---

**最后更新**：2025年1月
**版本**：v1.0