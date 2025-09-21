# 第5次课：转录组测序数据分析

## 课程概述

本次课程将深入学习转录组测序（RNA-seq）数据分析的完整流程，从原始数据处理到差异表达分析，再到功能富集分析。通过理论学习和实践操作，掌握现代转录组学分析的核心技术和方法。

## 学习目标

### 理论目标
- 理解RNA-seq技术原理和发展历程
- 掌握转录组实验设计的关键考虑因素
- 理解负二项分布模型在RNA-seq分析中的应用
- 掌握多重检验校正的原理和方法

### 实践目标
- 掌握转录组数据预处理和质量控制
- 学会使用HISAT2进行剪接感知序列比对
- 掌握featureCounts进行基因表达定量
- 学会使用DESeq2进行差异表达分析
- 掌握结果可视化和功能富集分析方法

## 课程内容

### 理论部分（约50张幻灯片）
1. **转录组测序技术概述**
   - 技术发展历程
   - RNA-seq vs 传统方法
   - 技术原理和流程

2. **RNA-seq实验设计**
   - 样本设计和重复
   - 测序深度选择
   - 链特异性测序

3. **数据预处理**
   - 质量控制流程
   - 数据清洗方法
   - 污染检测

4. **序列比对与定量**
   - 比对工具比较
   - 基因定量策略
   - 表达量标准化

5. **差异表达分析**
   - 统计学原理
   - DESeq2方法
   - 多重检验校正

6. **结果解读与可视化**
   - MA图和火山图
   - 热图分析
   - 功能富集分析

7. **高级应用**
   - 单细胞RNA-seq
   - 时间序列分析
   - 共表达网络

### 实践部分
1. **环境准备和数据获取**
   - 软件安装和配置
   - 参考数据下载
   - 目录结构创建

2. **数据质量控制**
   - FastQC质量评估
   - MultiQC报告生成
   - 数据清洗（如需要）

3. **序列比对**
   - HISAT2索引构建
   - 序列比对执行
   - BAM文件处理

4. **基因表达定量**
   - featureCounts计数
   - 计数矩阵生成
   - 质量检查

5. **差异表达分析**
   - DESeq2分析流程
   - 结果筛选和注释
   - 统计摘要生成

6. **结果可视化**
   - MA图和火山图
   - PCA和聚类分析
   - 差异基因热图

7. **功能富集分析**
   - GO富集分析
   - KEGG通路分析
   - 结果可视化

## 文件结构

```
lesson-05/
├── README.md                           # 课程说明文档
├── slides/                             # 课程幻灯片
│   └── slides.md                      # 50张理论幻灯片
├── manual/                             # 实践操作手册
│   └── manual.md                      # 详细操作指南
├── scripts/                            # 分析脚本集合
│   ├── rnaseq_pipeline.sh             # 完整分析流水线脚本
│   ├── deseq2_analysis.R              # DESeq2差异分析脚本
│   └── functional_enrichment.R        # 功能富集分析脚本
├── data/                               # 示例数据目录
│   ├── README.md                      # 数据说明文档
│   ├── sample_info.txt                # 样本信息表
│   └── *.fastq.gz                     # RNA-seq测序数据（需下载）
└── images/                             # 课程配图
    ├── rnaseq_timeline.svg            # RNA-seq技术发展时间线
    ├── rnaseq_workflow.svg            # RNA-seq实验流程图
    ├── sequencing_depth_detection.svg  # 测序深度与基因检出关系
    ├── strand_specific_sequencing.svg  # 链特异性测序原理
    ├── alignment_tools_comparison.svg  # 比对工具性能比较
    ├── quantification_methods.svg     # 定量方法分类
    ├── normalization_comparison.svg   # 标准化方法对比
    ├── negative_binomial_model.svg    # 负二项分布模型
    ├── multiple_testing_correction.svg # 多重检验校正效果
    ├── deseq2_results_plots.svg       # DESeq2结果可视化
    ├── go_hierarchy.svg               # GO层次结构示例
    └── scrna_workflow.svg             # 单细胞RNA-seq流程
```

## 软件环境要求

### 核心软件
| 软件名称 | 版本要求 | 用途 | 安装方式 |
|---------|---------|------|---------|
| FastQC | ≥0.11.9 | 质量控制 | `conda install fastqc` |
| MultiQC | ≥1.9 | 报告汇总 | `conda install multiqc` |
| Trimmomatic | ≥0.39 | 序列清洗 | `conda install trimmomatic` |
| HISAT2 | ≥2.2.1 | 序列比对 | `conda install hisat2` |
| SAMtools | ≥1.12 | BAM处理 | `conda install samtools` |
| featureCounts | ≥2.0.1 | 基因定量 | `conda install subread` |
| R | ≥4.1.0 | 统计分析 | `conda install r-base` |

### R包依赖
```r
# Bioconductor核心包
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db"))

# 可视化包
install.packages(c("ggplot2", "pheatmap", "RColorBrewer", "EnhancedVolcano"))
```

### 硬件要求
- **内存**: 至少8GB RAM（推荐16GB）
- **存储**: 至少20GB可用空间
- **CPU**: 4核心以上推荐
- **网络**: 稳定连接（用于下载参考数据）

## 数据准备

### 参考数据
- **人类基因组**: GRCh38 (~3GB)
- **基因注释**: Ensembl GTF v104 (~50MB)
- **下载地址**: ftp://ftp.ensembl.org/pub/release-104/

### 测试数据
- **样本数量**: 6个（3个处理组 + 3个对照组）
- **测序类型**: 双端150bp
- **数据量**: 每样本约100万reads
- **总大小**: 约2GB（压缩后）

## 预期学习成果

### 技能掌握
1. **数据处理能力**
   - 独立完成RNA-seq数据分析全流程
   - 熟练使用命令行工具和R语言
   - 掌握数据质量评估和问题诊断

2. **分析能力**
   - 理解差异表达分析的统计原理
   - 能够解释分析结果的生物学意义
   - 掌握结果可视化的最佳实践

3. **问题解决能力**
   - 能够诊断和解决常见分析问题
   - 掌握参数优化和结果验证方法
   - 具备独立学习新方法的能力

### 输出成果
1. **分析结果**
   - 差异表达基因列表（预期1000-3000个）
   - 标准化表达矩阵
   - 功能富集分析结果

2. **可视化图表**
   - 质量控制报告
   - MA图和火山图
   - PCA和聚类热图
   - 功能富集图表

3. **分析报告**
   - 方法描述
   - 结果解读
   - 生物学意义讨论

## 课后作业

### 必做作业
1. **完整分析流程**
   - 按照操作手册完成所有步骤
   - 生成完整的分析结果
   - 提交关键输出文件

2. **结果解读报告**
   - 分析结果摘要（500字）
   - 主要发现讨论（1000字）
   - 方法学反思（500字）

### 选做练习
1. **参数优化实验**
   - 比较不同比对参数的影响
   - 测试不同标准化方法
   - 分析批次效应校正

2. **扩展分析**
   - 时间序列分析（如有数据）
   - 共表达网络构建
   - 与公共数据库比较

### 思考题
1. 为什么RNA-seq数据适合用负二项分布建模？
2. 如何选择合适的差异表达分析阈值？
3. 功能富集分析结果如何与实验假设联系？
4. 单细胞RNA-seq与bulk RNA-seq有什么区别？

## 评分标准

| 评分项目 | 权重 | 评分标准 |
|---------|------|---------|
| 分析流程完整性 | 40% | 所有步骤正确执行，结果文件完整 |
| 结果解读准确性 | 30% | 对分析结果的理解和解释正确 |
| 报告撰写质量 | 20% | 逻辑清晰，表达准确，格式规范 |
| 创新性和深度 | 10% | 额外分析，深入思考，方法改进 |

## 常见问题解答

### Q1: 内存不足怎么办？
A: 可以减少线程数，使用更小的测试数据集，或者使用云计算资源。

### Q2: 软件安装失败怎么办？
A: 建议使用conda环境管理，按照官方文档逐步安装。

### Q3: 分析结果异常怎么办？
A: 首先检查输入数据质量，然后检查参数设置，最后查看日志文件。

### Q4: 如何解释生物学意义？
A: 结合实验设计、文献调研和功能富集结果进行综合分析。

## 扩展资源

### 在线教程
- [Bioconductor RNA-seq workflow](https://www.bioconductor.org/help/workflows/rnaseqGene/)
- [Galaxy RNA-seq tutorial](https://training.galaxyproject.org/training-material/topics/transcriptomics/)
- [HISAT2 manual](http://daehwankimlab.github.io/hisat2/)

### 参考文献
1. Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15, 550 (2014).
2. Kim, D., Langmead, B. & Salzberg, S.L. HISAT: a fast spliced aligner with low memory requirements. *Nature Methods* 12, 357–360 (2015).
3. Yu, G., Wang, L.G., Han, Y. & He, Q.Y. clusterProfiler: an R package for comparing biological themes among gene clusters. *OMICS* 16, 284-287 (2012).

### 相关工具
- [IGV](http://software.broadinstitute.org/software/igv/): 基因组浏览器
- [UCSC Genome Browser](https://genome.ucsc.edu/): 在线基因组浏览
- [STRING](https://string-db.org/): 蛋白质相互作用网络

## 技术支持

### 联系方式
- **主讲教师**: 王运生
- **邮箱**: wangys@hunau.edu.cn
- **办公室**: 16教420室

### 在线支持
- **课程QQ群**: [群号]
- **课程网站**: [网址]
- **GitHub仓库**: [仓库地址]

---

**课程更新**: 2024年版本  
**文档版本**: v1.0